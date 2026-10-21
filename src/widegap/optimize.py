"""Constrained doping optimization in 1D (stage 6).

docs/benchmarks/stage6_optimization.md. Maximizes a smoothed field-based
breakdown voltage BV_s of a p+ / n(design) / n+ punch-through diode at an
exactly fixed drift-resistance proxy R = int_0^W dx / b(x) = R0, using a
Legendre log-doping parametrization, the stage-5 adjoint, and a monotone
backtracking BFGS ascent. BV here is a proxy (no impact ionization).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from numpy.polynomial import legendre
from scipy.linalg import solve_banded
from scipy.optimize import brentq

from widegap.adjoint import banded_transpose, polish_reverse
from widegap.fv import reverse_bias as rb
from widegap.fv.equilibrium import control_volumes

_GX, _GW = np.polynomial.legendre.leggauss(8)
_QX, _QW = np.polynomial.legendre.leggauss(96)


@dataclass
class PTDiodeDesign:
    """p+ (a) on [-lp, 0], drift b(x) on [0, W], n+ (bn) on [W, W + ln]."""

    W: float = 12.0
    lp: float = 3.0
    ln: float = 3.0
    a: float = 100.0
    bn: float = 100.0
    R0: float = 12.0
    n_modes: int = 5
    h: float = 0.025
    delta: float = 1e-3
    x: np.ndarray = field(init=False)

    def __post_init__(self):
        n_el = int(round((self.lp + self.W + self.ln) / self.h))
        self.x = np.linspace(-self.lp, self.W + self.ln, n_el + 1)

    # g(x) = sum_k theta_k P_k(2 x / W - 1), k = 1..n_modes
    def _basis(self, xs):
        t = 2.0 * np.asarray(xs) / self.W - 1.0
        return np.stack(
            [legendre.legval(t, np.eye(self.n_modes + 1)[k]) for k in range(1, self.n_modes + 1)]
        )  # (n_modes, len(xs))

    def _norm(self, theta):
        """S = int_0^W e^{-g} dx and dS/dtheta (96-point Gauss)."""
        xq = 0.5 * self.W * (_QX + 1.0)
        wq = 0.5 * self.W * _QW
        P = self._basis(xq)
        eg = np.exp(-(theta @ P))
        return float(wq @ eg), -(P * eg) @ wq

    def b(self, theta, xs):
        """Drift doping b(x) and db/dtheta at points xs in [0, W]."""
        S, dS = self._norm(theta)
        P = self._basis(xs)
        bx = np.exp(theta @ P) * S / self.R0
        return bx, bx * (P + (dS / S)[:, None])

    def resistance(self, theta):
        xq = 0.5 * self.W * (_QX + 1.0)
        wq = 0.5 * self.W * _QW
        return float(wq @ (1.0 / self.b(theta, xq)[0]))

    def c_cv(self, theta):
        """Control-volume-averaged net doping C and dC/dtheta at the nodes."""
        x = self.x
        lo = np.r_[x[0], 0.5 * (x[:-1] + x[1:])]
        hi = np.r_[0.5 * (x[:-1] + x[1:]), x[-1]]
        length = hi - lo
        acc = -self.a * np.clip(np.minimum(hi, 0.0) - lo, 0.0, None)
        acc += self.bn * np.clip(hi - np.maximum(lo, self.W), 0.0, None)
        dacc = np.zeros((self.n_modes, len(x)))
        dlo, dhi = np.maximum(lo, 0.0), np.minimum(hi, self.W)
        idx = np.nonzero(dhi > dlo)[0]
        for i in idx:
            xs = 0.5 * (dhi[i] + dlo[i]) + 0.5 * (dhi[i] - dlo[i]) * _GX
            w = 0.5 * (dhi[i] - dlo[i]) * _GW
            bx, dbx = self.b(theta, xs)
            acc[i] += w @ bx
            dacc[:, i] = dbx @ w
        return acc / length, dacc / length

    def drift_nodes(self):
        return self.x[(self.x >= 0.0) & (self.x <= self.W)]


# ----- smoothed breakdown index and BV_s --------------------------------------
def smoothed_index(d: rb.Diode1D, psi, e_crit, kappa):
    """B_s = log(sum exp(kappa B_e)) / kappa and dB_s/dpsi."""
    e = d.element_field(psi)
    be = np.abs(e) / e_crit
    m = be.max()
    w = np.exp(kappa * (be - m))
    bs = m + np.log(w.sum()) / kappa
    sm = w / w.sum()
    h = np.diff(d.x)
    coef = sm * np.sign(e) / (h * e_crit)
    g = np.zeros_like(psi)
    g[1:] += -coef
    g[:-1] += coef
    return float(bs), g


def bv_smoothed(d: rb.Diode1D, e_crit=15.0, kappa=200.0, vmax=150, xtol=1e-12):
    """BV_s (B_s = 1) by continuation bracketing and Brent; returns (BV_s, psi)."""
    psi, prev = None, None
    for v in range(0, vmax + 1):
        psi, rep = d.solve(v, psi0=psi)
        if not rep.converged:
            raise RuntimeError(f"forward solve did not converge at V_R={v}")
        bs, _ = smoothed_index(d, psi, e_crit, kappa)
        if bs >= 1.0:
            if prev is None:
                raise ValueError("B_s >= 1 at V_R = 0")
            psi_lo = prev[1]

            def f(vr, psi_lo=psi_lo):
                p, _ = d.solve(vr, psi0=psi_lo)
                return smoothed_index(d, p, e_crit, kappa)[0] - 1.0

            bv = brentq(f, v - 1, v, xtol=xtol, rtol=4 * np.finfo(float).eps)
            p, _ = d.solve(bv, psi0=psi_lo)
            return bv, polish_reverse(d, p, bv)
        prev = (v, psi)
    raise ValueError("B_s = 1 not reached by vmax")


def bv_smoothed_grad_c(d: rb.Diode1D, psi, v_r, e_crit=15.0, kappa=200.0):
    """dBV_s/dC (interior nodes) via the adjoint and the implicit function theorem."""
    phi_n, phi_p = 0.0, -float(v_r)
    _, gpsi = smoothed_index(d, psi, e_crit, kappa)
    ab = rb.jacobian_banded(psi, d.x, d.delta, phi_n, phi_p)
    lam = solve_banded((1, 1), banded_transpose(ab), gpsi)
    vol = control_volumes(d.x)
    db_dc = lam * vol
    db_dc[0] = db_dc[-1] = 0.0
    f_v = np.zeros_like(psi)
    f_v[1:-1] = vol[1:-1] * d.delta * np.exp(phi_p - psi[1:-1])
    f_v[0] = 1.0
    return -db_dc / float(-lam @ f_v)


def evaluate(des: PTDiodeDesign, theta, e_crit=15.0, kappa=200.0, grad=True):
    c, dc = des.c_cv(theta)
    d = rb.Diode1D(des.x, c, des.delta)
    bv, psi = bv_smoothed(d, e_crit, kappa)
    if not grad:
        return bv, None, psi
    g_c = bv_smoothed_grad_c(d, psi, bv, e_crit, kappa)
    return bv, dc @ g_c, psi


# ----- optimizer -------------------------------------------------------------
@dataclass
class OptRecord:
    thetas: list = field(default_factory=list)
    bv_s: list = field(default_factory=list)
    grad_norm: list = field(default_factory=list)
    r_rel_err: list = field(default_factory=list)
    b_min: list = field(default_factory=list)
    b_max: list = field(default_factory=list)
    evaluations: int = 0
    stop_reason: str = ""


def feasible(des: PTDiodeDesign, theta, lo=0.05, hi=20.0):
    bx = des.b(theta, des.drift_nodes())[0]
    return bool(np.all(bx >= lo) and np.all(bx <= hi)), float(bx.min()), float(bx.max())


def optimize(des: PTDiodeDesign, theta0=None, max_iter=60, e_crit=15.0, kappa=200.0):
    """Monotone BFGS ascent with Armijo backtracking and bound rejection."""
    rec = OptRecord()
    th = np.zeros(des.n_modes) if theta0 is None else np.array(theta0, float)
    f, g, _ = evaluate(des, th, e_crit, kappa)
    rec.evaluations += 1
    H = np.eye(des.n_modes)

    def log(th, f, g):
        ok, bmin, bmax = feasible(des, th)
        rec.thetas.append(th.tolist())
        rec.bv_s.append(f)
        rec.grad_norm.append(float(np.linalg.norm(g)))
        rec.r_rel_err.append(abs(des.resistance(th) / des.R0 - 1.0))
        rec.b_min.append(bmin)
        rec.b_max.append(bmax)

    log(th, f, g)
    for _ in range(max_iter):
        if np.linalg.norm(g) <= 1e-6 * f:
            rec.stop_reason = "gradient"
            break
        if len(rec.bv_s) >= 4 and (rec.bv_s[-1] - rec.bv_s[-4]) <= 1e-8 * rec.bv_s[-1]:
            rec.stop_reason = "stagnation"
            break
        p = H @ g
        if p @ g <= 0:  # safeguard: fall back to steepest ascent
            H = np.eye(des.n_modes)
            p = g.copy()
        step, accepted = 1.0, False
        for _ in range(30):
            tn = th + step * p
            if feasible(des, tn)[0]:
                fn, gn, _ = evaluate(des, tn, e_crit, kappa)
                rec.evaluations += 1
                if fn >= f + 1e-4 * step * (p @ g) and fn > f:
                    accepted = True
                    break
            step *= 0.5
        if not accepted:
            rec.stop_reason = "line search"
            break
        s, y = tn - th, gn - g
        sy = s @ y
        if sy < 0:  # ascent: curvature condition is s.(-y) > 0 for -f
            rho = -1.0 / sy
            V = np.eye(des.n_modes) - rho * np.outer(s, -y)
            H = V @ H @ V.T + rho * np.outer(s, s)
        th, f, g = tn, fn, gn
        log(th, f, g)
    else:
        rec.stop_reason = "max_iter"
    return th, rec
