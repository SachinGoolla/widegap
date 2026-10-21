"""Reverse-bias electrostatics and the field-based breakdown index (stage 4).

Zero-current approximation (docs/benchmarks/stage4_breakdown_index.md): flat
quasi-Fermi potentials phi_n (cathode) and phi_p (anode),

    -psi'' = delta * (exp(phi_p - psi) - exp(psi - phi_n)) + C .

With phi_n = phi_p = 0 the residual is bit-identical to the stage-2 operator
(``widegap.fv.equilibrium.residual``). Breakdown index after Zhu & Andrei,
AIP Adv. 9, 095301 (2019), doi:10.1063/1.5113764, eq. (1):
B = max |E| / E_crit, with E taken on elements. The field-based BV (B = 1) is
a numerical proxy; no impact ionization is modelled.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.linalg import solve_banded
from scipy.optimize import brentq

from widegap.fv.equilibrium import control_volumes, neutral_potential

EPS_MACH = np.finfo(float).eps


def residual(psi, x, c_cv, delta, phi_n=0.0, phi_p=0.0):
    h = np.diff(x)
    vol = control_volumes(x)
    flux = -np.diff(psi) / h
    f = np.zeros_like(psi)
    pi = psi[1:-1]
    f[1:-1] = (flux[1:] - flux[:-1]) - vol[1:-1] * (
        delta * (np.exp(phi_p - pi) - np.exp(pi - phi_n)) + c_cv[1:-1]
    )
    return f


def jacobian_banded(psi, x, delta, phi_n=0.0, phi_p=0.0):
    g = 1.0 / np.diff(x)
    vol = control_volumes(x)
    pi = psi[1:-1]
    ab = np.zeros((3, len(x)))
    ab[1, 1:-1] = g[:-1] + g[1:] + vol[1:-1] * delta * (np.exp(phi_p - pi) + np.exp(pi - phi_n))
    ab[0, 2:] = -g[1:]
    ab[2, :-2] = -g[:-1]
    ab[1, 0] = ab[1, -1] = 1.0
    return ab


def roundoff_tolerance(psi, x, c_cv, delta, phi_n=0.0, phi_p=0.0):
    """1e2 * eps * max_i (sum_j |J_ij psi_j| + V_i |rho_i|) over interior rows."""
    ab = jacobian_banded(psi, x, delta, phi_n, phi_p)
    a = np.abs(psi)
    s = np.abs(ab[1]) * a
    s[:-1] += np.abs(ab[0, 1:]) * a[1:]
    s[1:] += np.abs(ab[2, :-1]) * a[:-1]
    pi = psi[1:-1]
    rho = np.abs(delta * (np.exp(phi_p - pi) - np.exp(pi - phi_n)) + c_cv[1:-1])
    return 1e2 * EPS_MACH * float(np.max(s[1:-1] + control_volumes(x)[1:-1] * rho))


@dataclass
class RBReport:
    converged: bool = False
    iterations: int = 0
    limited_steps: int = 0
    residual_history: list[float] = field(default_factory=list)
    tolerance: float = 0.0


@dataclass
class Diode1D:
    """Abrupt multi-layer diode: anode at x[0], cathode at x[-1]."""

    x: np.ndarray
    c_cv: np.ndarray
    delta: float

    def __post_init__(self):
        self.x = np.asarray(self.x, dtype=float)
        if self.x.ndim != 1 or len(self.x) < 3 or np.any(np.diff(self.x) <= 0):
            raise ValueError("x must be a strictly increasing 1D mesh with >= 3 nodes")
        self.c_cv = np.asarray(self.c_cv, dtype=float)
        self.psi_anode_eq = float(neutral_potential(self.c_cv[0], self.delta))
        self.psi_cathode_eq = float(neutral_potential(self.c_cv[-1], self.delta))

    def solve(self, v_r, psi0=None, max_iter=60, max_step=2.0):
        """Electrostatics at reverse bias v_r (anode at -v_r, phi_p = -v_r, phi_n = 0)."""
        phi_p, phi_n = -float(v_r), 0.0
        if psi0 is None:
            psi = neutral_potential(self.c_cv, self.delta)
        else:
            psi = np.array(psi0, dtype=float)
        psi[0] = self.psi_anode_eq - v_r
        psi[-1] = self.psi_cathode_eq
        args = (self.x, self.c_cv, self.delta, phi_n, phi_p)
        rep = RBReport()
        last_step = np.inf
        for _ in range(max_iter + 1):
            f = residual(psi, *args)
            rn = float(np.max(np.abs(f)))
            tol = roundoff_tolerance(psi, *args)
            rep.residual_history.append(rn)
            rep.tolerance = tol
            if rn <= tol and last_step <= 1e-8:
                rep.converged = True
                break
            if rep.iterations == max_iter:
                break
            dpsi = solve_banded((1, 1), jacobian_banded(psi, self.x, self.delta, phi_n, phi_p), -f)
            step = float(np.max(np.abs(dpsi)))
            if step > max_step:
                dpsi *= max_step / step
                rep.limited_steps += 1
            last_step = float(np.max(np.abs(dpsi)))
            psi = psi + dpsi
            rep.iterations += 1
        return psi, rep

    def element_field(self, psi):
        return -np.diff(psi) / np.diff(self.x)

    def max_field(self, psi):
        return float(np.max(np.abs(self.element_field(psi))))


@dataclass
class SweepResult:
    v: list[float]
    e_max: list[float]
    psi: list[np.ndarray]
    reports: list[RBReport]


def sweep(diode: Diode1D, v_values) -> SweepResult:
    """Bias continuation over increasing reverse biases (warm starts)."""
    out = SweepResult([], [], [], [])
    psi = None
    for v in v_values:
        psi, rep = diode.solve(v, psi0=psi)
        out.v.append(float(v))
        out.e_max.append(diode.max_field(psi))
        out.psi.append(psi)
        out.reports.append(rep)
    return out


def breakdown_voltage(diode: Diode1D, sw: SweepResult, e_crit: float, xtol=1e-10):
    """BV where B = E_max / e_crit = 1, bracketed on the sweep, refined by Brent."""
    b = np.asarray(sw.e_max) / e_crit
    k = int(np.argmax(b >= 1.0))
    if b[k] < 1.0 or k == 0:
        raise ValueError("B = 1 not bracketed by the sweep")
    psi_lo = sw.psi[k - 1]
    reports = []

    def g(v):
        psi, rep = diode.solve(v, psi0=psi_lo)
        reports.append(rep)
        return diode.max_field(psi) / e_crit - 1.0

    bv = brentq(g, sw.v[k - 1], sw.v[k], xtol=xtol, rtol=4 * EPS_MACH)
    psi_bv, rep = diode.solve(bv, psi0=psi_lo)
    reports.append(rep)
    return bv, psi_bv, reports
