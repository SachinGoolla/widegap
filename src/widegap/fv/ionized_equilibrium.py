"""Equilibrium Poisson with incomplete ionization (stage 7).

    -psi'' = delta (e^{-psi} - e^{psi}) + N_D^+(psi) - N_A^-(psi) + s

Box-method FV as in stage 2; nodal (control-volume-averaged) N_D and N_A;
Newton with analytic Jacobian, the stage-2 step limiter and the stage-3b
round-off-relative tolerance. Imports nothing from widegap.pinn.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.linalg import solve_banded

from widegap.fv.equilibrium import control_volumes
from widegap.materials.ionization import Level, charge, dcharge_dpsi

EPS_MACH = np.finfo(float).eps


@dataclass
class IonizedReport:
    converged: bool = False
    iterations: int = 0
    limited_steps: int = 0
    residual_history: list[float] = field(default_factory=list)
    tolerance: float = 0.0


@dataclass
class IonizedJunction:
    x: np.ndarray
    nd: np.ndarray  # scaled donor density at nodes (CV-averaged)
    na: np.ndarray  # scaled acceptor density at nodes (CV-averaged)
    delta: float
    donors: tuple[Level, ...]
    acceptors: tuple[Level, ...]
    source: np.ndarray | None = None

    def __post_init__(self):
        self.x = np.asarray(self.x, dtype=float)
        if self.x.ndim != 1 or len(self.x) < 3 or np.any(np.diff(self.x) <= 0):
            raise ValueError("x must be a strictly increasing 1D mesh with >= 3 nodes")
        self.vol = control_volumes(self.x)
        self.g = 1.0 / np.diff(self.x)
        if self.source is None:
            self.source = np.zeros_like(self.x)

    def rho(self, psi):
        return charge(psi, self.nd, self.na, self.delta, self.donors, self.acceptors) + self.source

    def residual(self, psi):
        flux = -np.diff(psi) * self.g
        f = np.zeros_like(psi)
        f[1:-1] = (flux[1:] - flux[:-1]) - self.vol[1:-1] * self.rho(psi)[1:-1]
        return f

    def jacobian_banded(self, psi):
        dr = dcharge_dpsi(psi, self.nd, self.na, self.delta, self.donors, self.acceptors)
        ab = np.zeros((3, len(psi)))
        ab[1, 1:-1] = self.g[:-1] + self.g[1:] - self.vol[1:-1] * dr[1:-1]
        ab[0, 2:] = -self.g[1:]
        ab[2, :-2] = -self.g[:-1]
        ab[1, 0] = ab[1, -1] = 1.0
        return ab

    def tolerance(self, psi):
        ab = self.jacobian_banded(psi)
        a = np.abs(psi)
        s = np.abs(ab[1]) * a
        s[:-1] += np.abs(ab[0, 1:]) * a[1:]
        s[1:] += np.abs(ab[2, :-1]) * a[:-1]
        mag = (
            self.delta * (np.exp(-psi) + np.exp(psi))
            + np.abs(self.nd)
            + np.abs(self.na)
            + np.abs(self.source)
        )
        return 1e2 * EPS_MACH * float(np.max((s + self.vol * mag)[1:-1]))

    def solve(self, psi0, psi_left, psi_right, max_iter=300, max_step=2.0):
        psi = np.array(psi0, dtype=float)
        psi[0], psi[-1] = psi_left, psi_right
        rep = IonizedReport()
        last = np.inf
        for _ in range(max_iter + 1):
            f = self.residual(psi)
            rn = float(np.max(np.abs(f)))
            tol = self.tolerance(psi)
            rep.residual_history.append(rn)
            rep.tolerance = tol
            if rn <= tol and last <= 1e-8:
                rep.converged = True
                break
            if rep.iterations == max_iter:
                break
            d = solve_banded((1, 1), self.jacobian_banded(psi), -f)
            step = float(np.max(np.abs(d)))
            if step > max_step:
                d *= max_step / step
                rep.limited_steps += 1
            last = float(np.max(np.abs(d)))
            psi = psi + d
            rep.iterations += 1
        return psi, rep
