"""Steady-state 1D drift-diffusion with Scharfetter-Gummel fluxes (stage 3).

Scaled model (docs/benchmarks/stage3_drift_diffusion.md, docs/equations.md):

    -psi'' = p - n + C + s_psi
     j_n'  =  R + s_n,   j_n =  mu_n (n' - n psi')
     j_p'  = -R + s_p,   j_p = -mu_p (p' + p psi')
     R = (n p - delta^2) / (tau_p (n + delta) + tau_n (p + delta))   (or R = 0)

Box-method finite volume; Scharfetter-Gummel fluxes (Scharfetter & Gummel,
IEEE Trans. Electron Devices 16, 64 (1969), doi:10.1109/T-ED.1969.16566);
node-lumped R and sources; fully coupled Newton with analytic Jacobian.
Unknowns are interleaved per node: X = [psi_0, n_0, p_0, psi_1, ...].
Imports nothing from ``widegap.pinn``.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve

from widegap.fv.equilibrium import control_volumes

EPS_MACH = np.finfo(float).eps


def bernoulli(t: np.ndarray) -> np.ndarray:
    """B(t) = t / (exp(t) - 1), with a series near 0."""
    t = np.asarray(t, dtype=float)
    out = np.empty_like(t)
    small = np.abs(t) < 1e-4
    ts = t[small]
    out[small] = 1.0 - ts / 2.0 + ts**2 / 12.0
    tl = t[~small]
    out[~small] = tl / np.expm1(tl)
    return out


def bernoulli_prime(t: np.ndarray) -> np.ndarray:
    """dB/dt = (expm1(t) - t e^t) / expm1(t)^2, with a series near 0."""
    t = np.asarray(t, dtype=float)
    out = np.empty_like(t)
    small = np.abs(t) < 1e-4
    ts = t[small]
    out[small] = -0.5 + ts / 6.0
    tl = t[~small]
    em = np.expm1(tl)
    out[~small] = (em - tl * np.exp(tl)) / em**2
    return out


@dataclass(frozen=True)
class DDParams:
    delta: float
    mu_n: float = 1.0
    mu_p: float = 1.0
    tau_n: float = 1.0
    tau_p: float = 1.0
    recombination: bool = True


@dataclass
class DDReport:
    converged: bool = False
    iterations: int = 0
    limited_steps: int = 0
    residual_history: list[tuple[float, float, float]] = field(default_factory=list)
    tolerance: tuple[float, float, float] = (0.0, 0.0, 0.0)
    last_dpsi: float = float("nan")


@dataclass
class DDSolution:
    x: np.ndarray
    psi: np.ndarray
    n: np.ndarray
    p: np.ndarray
    jn: np.ndarray  # element electron current
    jp: np.ndarray  # element hole current
    report: DDReport

    @property
    def j(self) -> np.ndarray:
        return self.jn + self.jp


class DriftDiffusion1D:
    """Residual, Jacobian and Newton solver for a fixed mesh and doping."""

    def __init__(self, x, c_cv, params: DDParams, sources=None):
        self.x = np.asarray(x, dtype=float)
        if self.x.ndim != 1 or len(self.x) < 3 or np.any(np.diff(self.x) <= 0):
            raise ValueError("x must be a strictly increasing 1D mesh with >= 3 nodes")
        self.c = np.asarray(c_cv, dtype=float)
        self.prm = params
        n = len(self.x)
        zero = np.zeros(n)
        s = sources or {}
        self.s_psi = np.asarray(s.get("psi", zero), dtype=float)
        self.s_n = np.asarray(s.get("n", zero), dtype=float)
        self.s_p = np.asarray(s.get("p", zero), dtype=float)
        self.h = np.diff(self.x)
        self.vol = control_volumes(self.x)

    # ----- physics pieces -------------------------------------------------
    def recombination(self, n, p):
        prm = self.prm
        if not prm.recombination:
            z = np.zeros_like(n)
            return z, z, z
        d = prm.tau_p * (n + prm.delta) + prm.tau_n * (p + prm.delta)
        num = n * p - prm.delta**2
        r = num / d
        dr_dn = (p * d - num * prm.tau_p) / d**2
        dr_dp = (n * d - num * prm.tau_n) / d**2
        return r, dr_dn, dr_dp

    def fluxes(self, psi, n, p):
        dl = np.diff(psi)
        bp, bm = bernoulli(dl), bernoulli(-dl)
        jn = self.prm.mu_n / self.h * (bp * n[1:] - bm * n[:-1])
        jp = -self.prm.mu_p / self.h * (bm * p[1:] - bp * p[:-1])
        return jn, jp

    # ----- residual and Jacobian -----------------------------------------
    @staticmethod
    def split(X):
        return X[0::3], X[1::3], X[2::3]

    def residual(self, X, bc):
        psi, n, p = self.split(X)
        fpsi_e = -np.diff(psi) / self.h
        jn, jp = self.fluxes(psi, n, p)
        r, _, _ = self.recombination(n, p)
        F = np.zeros_like(X)
        Fpsi, Fn, Fp = F[0::3], F[1::3], F[2::3]
        for Fb, fe in ((Fpsi, fpsi_e), (Fn, jn), (Fp, jp)):
            Fb[:-1] += fe
            Fb[1:] -= fe
        Fpsi -= self.vol * (p - n + self.c + self.s_psi)
        Fn -= self.vol * (r + self.s_n)
        Fp -= self.vol * (-r + self.s_p)
        for k in (0, -1):  # Dirichlet rows
            Fpsi[k], Fn[k], Fp[k] = psi[k] - bc[k][0], n[k] - bc[k][1], p[k] - bc[k][2]
        return F

    def term_scale(self, X):
        """Per-block S_b: max over rows of the sum of |terms| entering the row."""
        psi, n, p = self.split(X)
        dl = np.diff(psi)
        bp, bm = bernoulli(dl), bernoulli(-dl)
        tpsi = np.abs(np.diff(psi) / self.h)
        tn = self.prm.mu_n / self.h * (np.abs(bp * n[1:]) + np.abs(bm * n[:-1]))
        tp = self.prm.mu_p / self.h * (np.abs(bm * p[1:]) + np.abs(bp * p[:-1]))
        r, _, _ = self.recombination(n, p)
        out = []
        for te, vol_term in (
            (tpsi, np.abs(p) + np.abs(n) + np.abs(self.c) + np.abs(self.s_psi)),
            (tn, np.abs(r) + np.abs(self.s_n)),
            (tp, np.abs(r) + np.abs(self.s_p)),
        ):
            s = self.vol * vol_term
            s[:-1] += te
            s[1:] += te
            out.append(float(np.max(s[1:-1])))
        return tuple(out)

    def jacobian(self, X):
        psi, n, p = self.split(X)
        N = len(psi)
        mu_n, mu_p, h = self.prm.mu_n, self.prm.mu_p, self.h
        dl = np.diff(psi)
        bp, bm = bernoulli(dl), bernoulli(-dl)
        dbp, dbm = bernoulli_prime(dl), bernoulli_prime(-dl)
        a = np.arange(N - 1)  # left node of each element
        b = a + 1
        rows, cols, vals = [], [], []

        def add(r_node, r_var, c_node, c_var, v):
            rows.append(3 * r_node + r_var)
            cols.append(3 * c_node + c_var)
            vals.append(np.broadcast_to(v, r_node.shape))

        # Element flux contributions: row a gets +flux, row b gets -flux.
        g = 1.0 / h
        jn_dn_b, jn_dn_a = mu_n * g * bp, -mu_n * g * bm
        jn_dD = mu_n * g * (dbp * n[1:] + dbm * n[:-1])
        jp_dp_b, jp_dp_a = -mu_p * g * bm, mu_p * g * bp
        jp_dD = mu_p * g * (dbm * p[1:] + dbp * p[:-1])
        for sign, row in ((1.0, a), (-1.0, b)):
            # psi block: flux = -(psi_b - psi_a)/h
            add(row, 0, a, 0, sign * g)
            add(row, 0, b, 0, -sign * g)
            # electron
            add(row, 1, b, 1, sign * jn_dn_b)
            add(row, 1, a, 1, sign * jn_dn_a)
            add(row, 1, b, 0, sign * jn_dD)
            add(row, 1, a, 0, -sign * jn_dD)
            # hole
            add(row, 2, b, 2, sign * jp_dp_b)
            add(row, 2, a, 2, sign * jp_dp_a)
            add(row, 2, b, 0, sign * jp_dD)
            add(row, 2, a, 0, -sign * jp_dD)

        # Node (volume) terms.
        i = np.arange(N)
        _, dr_dn, dr_dp = self.recombination(n, p)
        vol = self.vol
        add(i, 0, i, 1, vol)  # -vol*(p - n): d/dn = +vol
        add(i, 0, i, 2, -vol)
        add(i, 1, i, 1, -vol * dr_dn)
        add(i, 1, i, 2, -vol * dr_dp)
        add(i, 2, i, 1, vol * dr_dn)
        add(i, 2, i, 2, vol * dr_dp)

        R = np.concatenate(rows)
        C = np.concatenate(cols)
        V = np.concatenate(vals)
        bc_rows = np.r_[0, 1, 2, 3 * N - 3, 3 * N - 2, 3 * N - 1]
        keep = ~np.isin(R, bc_rows)
        R = np.r_[R[keep], bc_rows]
        C = np.r_[C[keep], bc_rows]
        V = np.r_[V[keep], np.ones(6)]
        return sp.csc_matrix((V, (R, C)), shape=(3 * N, 3 * N))

    # ----- Newton ---------------------------------------------------------
    def solve(self, X0, bc, max_iter=60, max_step=2.0):
        X = np.array(X0, dtype=float)
        for k in (0, -1):
            j = 0 if k == 0 else len(X) // 3 - 1
            X[3 * j : 3 * j + 3] = bc[k]
        rep = DDReport()
        for _ in range(max_iter + 1):
            F = self.residual(X, bc)
            norms = tuple(float(np.max(np.abs(F[v::3]))) for v in range(3))
            tol = tuple(1e2 * EPS_MACH * s for s in self.term_scale(X))
            rep.residual_history.append(norms)
            rep.tolerance = tol
            if all(nb <= tb for nb, tb in zip(norms, tol, strict=True)) and rep.last_dpsi <= 1e-8:
                rep.converged = True
                break
            if rep.iterations == max_iter:
                break
            dX = spsolve(self.jacobian(X), -F)
            step = float(np.max(np.abs(dX[0::3])))
            if step > max_step:
                dX *= max_step / step
                rep.limited_steps += 1
            rep.last_dpsi = float(np.max(np.abs(dX[0::3])))
            X += dX
            rep.iterations += 1
        psi, n, p = self.split(X)
        jn, jp = self.fluxes(psi, n, p)
        return DDSolution(self.x, psi.copy(), n.copy(), p.copy(), jn, jp, rep)
