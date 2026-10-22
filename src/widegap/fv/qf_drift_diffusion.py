"""Quasi-Fermi-variable drift-diffusion, 1D (stage 7.1).

Same physics and box/Scharfetter-Gummel discretization as stage 3
(``widegap.fv.drift_diffusion``), but with unknowns X = (psi, phi_n, phi_p)
per node (interleaved) and

    n = exp(log_delta + psi - phi_n),   p = exp(log_delta + phi_p - psi).

The fluxes are algebraically identical to stage 3 but are written without
cancellation, using expm1 on quasi-Fermi differences, so that minority
quantities of order delta^2 ~ 1e-48 remain resolved
(docs/benchmarks/stage7_1_quasi_fermi_dd.md). ``log_delta`` is passed
explicitly so that delta = 1e-24 never has to be formed in a product with
exponentials that over/underflow; ``params.delta`` must equal
``exp(log_delta)`` and is used only where delta appears additively (SRH).
Imports nothing from ``widegap.pinn``.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
from scipy.special import logsumexp

from widegap.fv.drift_diffusion import DDParams, bernoulli, bernoulli_prime
from widegap.fv.equilibrium import control_volumes

EPS_MACH = np.finfo(float).eps


@dataclass
class QFReport:
    converged: bool = False
    iterations: int = 0
    limited_steps: int = 0
    residual_history: list[float] = field(default_factory=list)  # max_i |F_i| / S_i
    last_dX: float = float("inf")


@dataclass
class QFSolution:
    x: np.ndarray
    X: np.ndarray
    psi: np.ndarray
    phi_n: np.ndarray
    phi_p: np.ndarray
    n: np.ndarray
    p: np.ndarray
    jn: np.ndarray
    jp: np.ndarray
    report: QFReport

    @property
    def j(self) -> np.ndarray:
        return self.jn + self.jp


class QFDriftDiffusion1D:
    """Residual, Jacobian and Newton solver for a fixed mesh and doping."""

    def __init__(self, x, c_cv, params: DDParams, log_delta: float, sources=None):
        self.x = np.asarray(x, dtype=float)
        if self.x.ndim != 1 or len(self.x) < 3 or np.any(np.diff(self.x) <= 0):
            raise ValueError("x must be a strictly increasing 1D mesh with >= 3 nodes")
        if not np.isclose(params.delta, np.exp(log_delta), rtol=1e-12, atol=0.0):
            raise ValueError("params.delta must equal exp(log_delta)")
        self.c = np.asarray(c_cv, dtype=float)
        self.prm = params
        self.ld = float(log_delta)
        n = len(self.x)
        zero = np.zeros(n)
        s = sources or {}
        self.s_psi = np.asarray(s.get("psi", zero), dtype=float)
        self.s_n = np.asarray(s.get("n", zero), dtype=float)
        self.s_p = np.asarray(s.get("p", zero), dtype=float)
        self.h = np.diff(self.x)
        self.vol = control_volumes(self.x)

    @staticmethod
    def split(X):
        return X[0::3], X[1::3], X[2::3]

    def densities(self, X):
        psi, pn, pp = self.split(X)
        return np.exp(self.ld + psi - pn), np.exp(self.ld + pp - psi)

    # ----- SRH -------------------------------------------------------------
    def _srh(self, X):
        """R and (dR/dpsi, dR/dphin, dR/dphip) at nodes."""
        psi, pn, pp = self.split(X)
        z = np.zeros_like(psi)
        if not self.prm.recombination:
            return z, (z, z, z)
        prm = self.prm
        n, p = self.densities(X)
        d2 = np.exp(2.0 * self.ld)
        ex = np.exp(pp - pn)
        num = d2 * np.expm1(pp - pn)
        dnum_p = d2 * ex
        D = prm.tau_p * (n + prm.delta) + prm.tau_n * (p + prm.delta)
        R = num / D
        dD_psi = prm.tau_p * n - prm.tau_n * p
        dD_pn = -prm.tau_p * n
        dD_pp = prm.tau_n * p
        dR_psi = -num * dD_psi / D**2
        dR_pn = (-dnum_p * D - num * dD_pn) / D**2
        dR_pp = (dnum_p * D - num * dD_pp) / D**2
        return R, (dR_psi, dR_pn, dR_pp)

    # ----- element fluxes ---------------------------------------------------
    def _flux_parts(self, X):
        psi, pn, pp = self.split(X)
        D = np.diff(psi)
        g = 1.0 / self.h
        B, Bm = bernoulli(D), bernoulli(-D)
        A = np.exp(self.ld + psi[1:] - pn[:-1])
        E = np.expm1(-(pn[1:] - pn[:-1]))
        P = np.exp(self.ld + pp[:-1] - psi[1:])
        Fp = np.expm1(pp[1:] - pp[:-1])
        return D, g, B, Bm, A, E, P, Fp

    def fluxes(self, X):
        _, g, B, Bm, A, E, P, Fp = self._flux_parts(X)
        return self.prm.mu_n * g * B * A * E, -self.prm.mu_p * g * Bm * P * Fp

    # ----- residual / Jacobian ---------------------------------------------
    def residual(self, X, bc):
        psi, pn, pp = self.split(X)
        n, p = self.densities(X)
        fpsi = -np.diff(psi) / self.h
        jn, jp = self.fluxes(X)
        R, _ = self._srh(X)
        F = np.zeros_like(X)
        Fpsi, Fn, Fp = F[0::3], F[1::3], F[2::3]
        for Fb, fe in ((Fpsi, fpsi), (Fn, jn), (Fp, jp)):
            Fb[:-1] += fe
            Fb[1:] -= fe
        Fpsi -= self.vol * (p - n + self.c + self.s_psi)
        Fn -= self.vol * (R + self.s_n)
        Fp -= self.vol * (-R + self.s_p)
        for k in (0, -1):
            Fpsi[k], Fn[k], Fp[k] = psi[k] - bc[k][0], pn[k] - bc[k][1], pp[k] - bc[k][2]
        return F

    def jacobian(self, X):
        psi, pn, pp = self.split(X)
        N = len(psi)
        mu_n, mu_p = self.prm.mu_n, self.prm.mu_p
        D, g, B, Bm, A, E, P, Fp = self._flux_parts(X)
        Bp_, Bmp = bernoulli_prime(D), bernoulli_prime(-D)
        a = np.arange(N - 1)
        b = a + 1
        rows, cols, vals = [], [], []

        def add(r_node, r_var, c_node, c_var, v):
            rows.append(3 * r_node + r_var)
            cols.append(3 * c_node + c_var)
            vals.append(np.broadcast_to(v, r_node.shape))

        jn_psib = mu_n * g * (Bp_ + B) * A * E
        jn_psia = -mu_n * g * Bp_ * A * E
        jn_pnb = -mu_n * g * B * A * (E + 1.0)
        jn_pna = mu_n * g * B * A
        jp_psib = mu_p * g * (Bmp + Bm) * P * Fp
        jp_psia = -mu_p * g * Bmp * P * Fp
        jp_ppb = -mu_p * g * Bm * P * (Fp + 1.0)
        jp_ppa = mu_p * g * Bm * P
        for sign, row in ((1.0, a), (-1.0, b)):
            add(row, 0, a, 0, sign * g)
            add(row, 0, b, 0, -sign * g)
            add(row, 1, b, 0, sign * jn_psib)
            add(row, 1, a, 0, sign * jn_psia)
            add(row, 1, b, 1, sign * jn_pnb)
            add(row, 1, a, 1, sign * jn_pna)
            add(row, 2, b, 0, sign * jp_psib)
            add(row, 2, a, 0, sign * jp_psia)
            add(row, 2, b, 2, sign * jp_ppb)
            add(row, 2, a, 2, sign * jp_ppa)

        i = np.arange(N)
        n, p = self.densities(X)
        vol = self.vol
        add(i, 0, i, 0, vol * (p + n))
        add(i, 0, i, 1, -vol * n)
        add(i, 0, i, 2, -vol * p)
        if self.prm.recombination:
            _, (r_psi, r_pn, r_pp) = self._srh(X)
            for cv, dr in ((0, r_psi), (1, r_pn), (2, r_pp)):
                add(i, 1, i, cv, -vol * dr)
                add(i, 2, i, cv, vol * dr)

        R_ = np.concatenate(rows)
        C_ = np.concatenate(cols)
        V_ = np.concatenate(vals)
        bc_rows = np.r_[0, 1, 2, 3 * N - 3, 3 * N - 2, 3 * N - 1]
        keep = ~np.isin(R_, bc_rows)
        R_ = np.r_[R_[keep], bc_rows]
        C_ = np.r_[C_[keep], bc_rows]
        V_ = np.r_[V_[keep], np.ones(6)]
        return sp.csc_matrix((V_, (R_, C_)), shape=(3 * N, 3 * N))

    def _volume_terms(self, X):
        n, p = self.densities(X)
        R, _ = self._srh(X)
        v = self.vol
        out = np.empty_like(X)
        out[0::3] = v * (np.abs(p) + np.abs(n) + np.abs(self.c) + np.abs(self.s_psi))
        out[1::3] = v * (np.abs(R) + np.abs(self.s_n))
        out[2::3] = v * (np.abs(R) + np.abs(self.s_p))
        return out

    def row_scale(self, X, J=None):
        """S_i = sum_j |J_ij| max(|X_j|, 1) + V_i |volume terms of row i|."""
        J = self.jacobian(X) if J is None else J
        return abs(J) @ np.maximum(np.abs(X), 1.0) + self._volume_terms(X)

    # ----- Newton ------------------------------------------------------------
    def solve(self, X0, bc, max_iter=60, max_step=2.0):
        X = np.array(X0, dtype=float)
        N = len(X) // 3
        X[:3] = bc[0]
        X[3 * N - 3 :] = bc[-1]
        rep = QFReport()
        for _ in range(max_iter + 1):
            F = self.residual(X, bc)
            J = self.jacobian(X)
            S = self.row_scale(X, J)
            ratio = np.abs(F[3:-3]) / S[3:-3]
            rep.residual_history.append(float(np.max(ratio)))
            if np.all(ratio <= 1e2 * EPS_MACH) and rep.last_dX <= 1e-8:
                rep.converged = True
                break
            if rep.iterations == max_iter:
                break
            rs = np.asarray(abs(J).max(axis=1).todense()).ravel()
            Js = sp.diags(1.0 / rs) @ J
            dX = spsolve(sp.csc_matrix(Js), -F / rs)
            step = float(np.max(np.abs(dX[0::3])))
            if step > max_step:
                dX *= max_step / step
                rep.limited_steps += 1
            rep.last_dX = float(np.max(np.abs(dX)))
            X += dX
            rep.iterations += 1
        psi, pn, pp = self.split(X)
        n, p = self.densities(X)
        jn, jp = self.fluxes(X)
        return QFSolution(self.x, X.copy(), psi.copy(), pn.copy(), pp.copy(), n, p, jn, jp, rep)


def pack(psi, phi_n, phi_p):
    X = np.empty(3 * len(psi))
    X[0::3], X[1::3], X[2::3] = psi, phi_n, phi_p
    return X


def exact_sg_identity_currents(
    x, psi, phi_n_left, phi_n_right, phi_p_left, phi_p_right, log_delta, mu_n, mu_p
):
    """Exact discrete SG currents (jn, jp) for R == 0, given psi.

    Electrons: on element e=(a,b), D=psi_b-psi_a,
      jn_e = mu_n/h B(D) delta e^{psi_b} (e^{-phi_n,b} - e^{-phi_n,a}),
    so (e^{-phi_n,b}-e^{-phi_n,a}) = jn h/(mu_n delta B(D) e^{psi_b}); with jn
    constant (R=0) the sum telescopes:
      e^{-phi_n,R} - e^{-phi_n,L} = jn * Sn,  Sn = sum_e h/(mu_n delta B(D) e^{psi_b}).
    Holes: jp_e = -mu_p/h B(-D) delta e^{-psi_b} (e^{phi_p,b} - e^{phi_p,a}), so
      e^{phi_p,R} - e^{phi_p,L} = -jp * Sp,  Sp = sum_e h/(mu_p delta B(-D) e^{-psi_b}).
    Sums use logsumexp; numerators use expm1 of the quasi-Fermi difference:
      jn = e^{-phi_nL} expm1(-(phi_nR-phi_nL)) / Sn
      jp = -e^{phi_pL} expm1(phi_pR-phi_pL) / Sp.
    """
    x = np.asarray(x, dtype=float)
    psi = np.asarray(psi, dtype=float)
    h = np.diff(x)
    D = np.diff(psi)
    psib = psi[1:]
    lsn = logsumexp(np.log(h) - np.log(mu_n) - log_delta - np.log(bernoulli(D)) - psib)
    lsp = logsumexp(np.log(h) - np.log(mu_p) - log_delta - np.log(bernoulli(-D)) + psib)
    jn = np.exp(-phi_n_left - lsn) * np.expm1(-(phi_n_right - phi_n_left))
    jp = -np.exp(phi_p_left - lsp) * np.expm1(phi_p_right - phi_p_left)
    return float(jn), float(jp)
