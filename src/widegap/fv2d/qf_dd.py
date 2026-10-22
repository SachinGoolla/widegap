"""Quasi-Fermi drift-diffusion on the edge graph (stage-7.1 physics, 1D and 2D).

Identical to ``widegap.fv.qf_drift_diffusion`` edge by edge (SG fluxes in
expm1 form, node-lumped SRH, analytic Jacobian, per-row tolerance,
row-equilibrated Newton), except that each edge flux is multiplied by the
edge's dual-face weight ``mesh.w`` (1D: 1) and the Poisson coupling is the
edge conductance ``G_e``.  Unknowns X = (psi, phi_n, phi_p) interleaved per
node.  Dirichlet (ohmic) rows are identity rows on chosen node sets; all other
boundaries carry zero flux.  Semiconductor everywhere, constant mobility.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve

from widegap.fv.drift_diffusion import DDParams, bernoulli, bernoulli_prime
from widegap.fv2d.core import edge_conductance
from widegap.fv2d.mesh import TensorMesh

EPS_MACH = np.finfo(float).eps


@dataclass
class QFReport:
    converged: bool = False
    iterations: int = 0
    limited_steps: int = 0
    residual_history: list[float] = field(default_factory=list)
    last_dX: float = float("inf")


@dataclass
class QFSolution:
    X: np.ndarray
    psi: np.ndarray
    phi_n: np.ndarray
    phi_p: np.ndarray
    n: np.ndarray
    p: np.ndarray
    jn: np.ndarray  # per-edge electron current through the dual face (a -> b)
    jp: np.ndarray
    report: QFReport

    @property
    def j(self) -> np.ndarray:
        return self.jn + self.jp


def pack(psi, phi_n, phi_p):
    X = np.empty(3 * len(psi))
    X[0::3], X[1::3], X[2::3] = psi, phi_n, phi_p
    return X


class QFDriftDiffusion:
    def __init__(
        self, mesh: TensorMesh, c_nodes, params: DDParams, log_delta, eps_cell=1.0, sources=None
    ):
        if not np.isclose(params.delta, np.exp(log_delta), rtol=1e-12, atol=0.0):
            raise ValueError("params.delta must equal exp(log_delta)")
        self.mesh = mesh
        self.c = np.asarray(c_nodes, dtype=float)
        self.prm = params
        self.ld = float(log_delta)
        n = mesh.n_nodes
        zero = np.zeros(n)
        s = sources or {}
        self.s_psi = np.asarray(s.get("psi", zero), dtype=float)
        self.s_n = np.asarray(s.get("n", zero), dtype=float)
        self.s_p = np.asarray(s.get("p", zero), dtype=float)
        self.a, self.b = mesh.ea, mesh.eb
        self.g = 1.0 / mesh.h
        self.gw = self.g * mesh.w  # SG coupling x dual-face weight
        self.G = edge_conductance(mesh, eps_cell)  # Poisson coupling
        self.vol = mesh.vol

    @staticmethod
    def split(X):
        return X[0::3], X[1::3], X[2::3]

    def densities(self, X):
        psi, pn, pp = self.split(X)
        return np.exp(self.ld + psi - pn), np.exp(self.ld + pp - psi)

    def _srh(self, X):
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

    def _flux_parts(self, X):
        psi, pn, pp = self.split(X)
        a, b = self.a, self.b
        D = psi[b] - psi[a]
        B, Bm = bernoulli(D), bernoulli(-D)
        A = np.exp(self.ld + psi[b] - pn[a])
        E = np.expm1(-(pn[b] - pn[a]))
        P = np.exp(self.ld + pp[a] - psi[b])
        Fp = np.expm1(pp[b] - pp[a])
        return D, B, Bm, A, E, P, Fp

    def fluxes(self, X):
        _, B, Bm, A, E, P, Fp = self._flux_parts(X)
        gw = self.gw
        return self.prm.mu_n * gw * B * A * E, -self.prm.mu_p * gw * Bm * P * Fp

    # ------------------------------------------------------------------
    def residual(self, X, dir_nodes, dir_vals):
        psi, pn, pp = self.split(X)
        n, p = self.densities(X)
        a, b = self.a, self.b
        fpsi = -(psi[b] - psi[a]) * self.G
        jn, jp = self.fluxes(X)
        R, _ = self._srh(X)
        F = np.zeros_like(X)
        Fpsi, Fn, Fp = F[0::3], F[1::3], F[2::3]
        for Fb, fe in ((Fpsi, fpsi), (Fn, jn), (Fp, jp)):
            Fb += np.bincount(a, fe, minlength=len(psi))
            Fb -= np.bincount(b, fe, minlength=len(psi))
        Fpsi -= self.vol * (p - n + self.c + self.s_psi)
        Fn -= self.vol * (R + self.s_n)
        Fp -= self.vol * (-R + self.s_p)
        dn = np.asarray(dir_nodes, dtype=int)
        dv = np.asarray(dir_vals, dtype=float).reshape(len(dn), 3)
        Fpsi[dn] = psi[dn] - dv[:, 0]
        Fn[dn] = pn[dn] - dv[:, 1]
        Fp[dn] = pp[dn] - dv[:, 2]
        return F

    def jacobian(self, X, dir_nodes):
        psi, pn, pp = self.split(X)
        N = len(psi)
        mu_n, mu_p = self.prm.mu_n, self.prm.mu_p
        D, B, Bm, A, E, P, Fp = self._flux_parts(X)
        Bp_, Bmp = bernoulli_prime(D), bernoulli_prime(-D)
        a, b, g, G = self.a, self.b, self.gw, self.G
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
            add(row, 0, a, 0, sign * G)
            add(row, 0, b, 0, -sign * G)
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
        dn = np.asarray(dir_nodes, dtype=int)
        bc_rows = np.sort(np.concatenate([3 * dn, 3 * dn + 1, 3 * dn + 2]))
        keep = ~np.isin(R_, bc_rows)
        R_ = np.r_[R_[keep], bc_rows]
        C_ = np.r_[C_[keep], bc_rows]
        V_ = np.r_[V_[keep], np.ones(len(bc_rows))]
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

    def row_scale(self, X, J):
        """S_i = sum_j |J_ij| max(|X_j|, 1) + V_i |volume terms of row i|."""
        return abs(J) @ np.maximum(np.abs(X), 1.0) + self._volume_terms(X)

    def solve(self, X0, dir_nodes, dir_vals, max_iter=60, max_step=2.0):
        dn = np.asarray(dir_nodes, dtype=int)
        dv = np.asarray(dir_vals, dtype=float).reshape(len(dn), 3)
        X = np.array(X0, dtype=float)
        X[3 * dn], X[3 * dn + 1], X[3 * dn + 2] = dv[:, 0], dv[:, 1], dv[:, 2]
        free = np.ones(len(X), dtype=bool)
        free[np.concatenate([3 * dn, 3 * dn + 1, 3 * dn + 2])] = False
        rep = QFReport()
        for _ in range(max_iter + 1):
            F = self.residual(X, dn, dv)
            J = self.jacobian(X, dn)
            S = self.row_scale(X, J)
            ratio = np.abs(F[free]) / S[free]
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
        return QFSolution(X.copy(), psi.copy(), pn.copy(), pp.copy(), n, p, jn, jp, rep)

    # ----- terminal currents -------------------------------------------
    def contact_current(self, sol: QFSolution, nodes) -> float:
        """Total current (jn + jp) leaving the node set into the rest of the
        device, summed over edges crossing the set boundary (unit depth)."""
        inset = np.zeros(self.mesh.n_nodes, dtype=bool)
        inset[np.asarray(nodes, dtype=int)] = True
        f = sol.j
        out = inset[self.a] & ~inset[self.b]
        inn = ~inset[self.a] & inset[self.b]
        return float(f[out].sum() - f[inn].sum())

    def cut_current(self, sol: QFSolution, axis: int, k: int) -> float:
        """Current in +axis direction through the cut between grid lines k and k+1."""
        m = self.mesh
        coord = m.x if axis == 0 else m.y
        xa, xb = m.node_xy[self.a, axis], m.node_xy[self.b, axis]
        sel = (m.edge_dir == axis) & (xa == coord[k]) & (xb == coord[k + 1])
        return float(sol.j[sel].sum())
