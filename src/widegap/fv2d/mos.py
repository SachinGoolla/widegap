"""Region-aware quasi-Fermi drift-diffusion for MOS structures (stage 8.1).

Unknowns X = (psi, phi_n, phi_p) interleaved per node (stage 7.1/8):

    n = delta exp(psi - phi_n),  p = delta exp(phi_p - psi),  rho = p - n + C.

* Poisson holds on *all* nodes with cell-wise eps (semiconductor eps_s = 1,
  oxide eps_ox/eps_s = ``eps_ratio``); the space charge acts only over the
  semiconductor part of each control volume (``Regions.vol_s``).
* Carrier equations hold on semiconductor nodes with semiconductor-only
  dual-face weights; nodes touching no semiconductor cell carry identity rows
  phi_n = phi_p = 0.  Interface/oxide edges have zero weight, hence zero flux.
* Dirichlet data: ohmic contacts (psi = psi_eq + V, phi_n = phi_p = V),
  gate (psi = V_G - V_FB), and, for equilibrium, phi = 0 on every
  semiconductor node (the special case phi = 0).
* Newton: analytic Jacobian (parent class assembly with masked weights),
  per-row round-off tolerance |F_i| <= 1e2 eps S_i, row equilibration,
  psi step limiter, bias continuation with step halving.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import scipy.sparse as sp

from widegap.fv.drift_diffusion import DDParams, bernoulli
from widegap.fv2d.qf_dd import EPS_MACH, QFDriftDiffusion, QFReport, QFSolution
from widegap.fv2d.regions import Regions
from widegap.linalg import safe_spsolve


class Dirichlet:
    """Row-wise Dirichlet data on X (row = 3*node + var, var 0/1/2 = psi/phi_n/phi_p)."""

    def __init__(self):
        self._d: dict[int, float] = {}

    def _set(self, rows, vals):
        rows = np.atleast_1d(rows)
        vals = np.broadcast_to(np.asarray(vals, dtype=float), rows.shape)
        for r, v in zip(rows, vals, strict=True):
            self._d[int(r)] = float(v)

    def psi(self, nodes, vals):
        nodes = np.asarray(nodes, dtype=int)
        self._set(3 * nodes, vals)
        return self

    def phi(self, nodes, vn, vp=None):
        nodes = np.asarray(nodes, dtype=int)
        self._set(3 * nodes + 1, vn)
        self._set(3 * nodes + 2, vn if vp is None else vp)
        return self

    def ohmic(self, c_nodes, delta, nodes, v=0.0):
        """psi = asinh(C/2delta) + V, phi_n = phi_p = V."""
        nodes = np.asarray(nodes, dtype=int)
        psi = np.arcsinh(np.asarray(c_nodes)[nodes] / (2.0 * delta)) + v
        return self.psi(nodes, psi).phi(nodes, v)

    def gate(self, nodes, vg, vfb):
        return self.psi(nodes, vg - vfb)

    def arrays(self):
        rows = np.array(sorted(self._d), dtype=int)
        return rows, np.array([self._d[r] for r in rows])


class MOSDriftDiffusion(QFDriftDiffusion):
    def __init__(
        self,
        regions: Regions,
        c_nodes,
        params: DDParams,
        log_delta,
        eps_ratio=0.38,
        sources=None,
    ):
        mesh = regions.mesh
        eps_cell = regions.cell_values(1.0, eps_ratio)
        if params.recombination:
            raise ValueError("stage 8.1 uses R = 0 (recombination=False)")
        super().__init__(mesh, c_nodes, params, log_delta, eps_cell=eps_cell, sources=sources)
        self.reg = regions
        self.gw = self.g * regions.w_s
        self.vol = regions.vol_s
        self.c = np.where(regions.semi_node, self.c, 0.0)
        ox = np.flatnonzero(regions.oxide_node)
        self.oxide_rows = np.r_[3 * ox + 1, 3 * ox + 2]

    # ----- exponentials only on semiconductor nodes (oxide psi can be huge) ---
    def _psi_c(self, psi):
        return np.where(self.reg.semi_node, psi, 0.0)

    def densities(self, X):
        psi, pn, pp = self.split(X)
        psi = self._psi_c(psi)
        return np.exp(self.ld + psi - pn), np.exp(self.ld + pp - psi)

    def _flux_parts(self, X):
        psi, pn, pp = self.split(X)
        psi = self._psi_c(psi)
        a, b = self.a, self.b
        D = psi[b] - psi[a]
        B, Bm = bernoulli(D), bernoulli(-D)
        A = np.exp(self.ld + psi[b] - pn[a])
        E = np.expm1(-(pn[b] - pn[a]))
        P = np.exp(self.ld + pp[a] - psi[b])
        Fp = np.expm1(pp[b] - pp[a])
        return D, B, Bm, A, E, P, Fp

    # ------------------------------------------------------------------
    def _all_bc(self, bc: Dirichlet):
        rows, vals = bc.arrays()
        rows = np.r_[rows, self.oxide_rows].astype(int)
        vals = np.r_[vals, np.zeros(len(self.oxide_rows))]
        order = np.argsort(rows, kind="stable")
        rows, vals = rows[order], vals[order]
        keep = np.r_[rows[1:] != rows[:-1], True]  # identity rows (last) win
        return rows[keep], vals[keep]

    def equilibrium_bc(self, bc: Dirichlet) -> Dirichlet:
        """Copy of ``bc`` with phi = 0 on all semiconductor nodes (equilibrium)."""
        out = Dirichlet()
        out._d = dict(bc._d)
        out.phi(np.flatnonzero(self.reg.semi_node), 0.0)
        return out

    def initial_guess(self, bc: Dirichlet | None = None):
        n = self.mesh.n_nodes
        psi = np.where(self.reg.semi_node, np.arcsinh(self.c / (2.0 * np.exp(self.ld))), 0.0)
        X = np.zeros(3 * n)
        X[0::3] = psi
        if bc is not None:
            rows, vals = self._all_bc(bc)
            X[rows] = vals
        return X

    def residual_bc(self, X, rows, vals):
        F = super().residual(X, [], [])
        F[rows] = X[rows] - vals
        return F

    def jacobian_bc(self, X, rows):
        J = sp.csr_matrix(super().jacobian(X, []))
        mask = np.ones(J.shape[0])
        mask[rows] = 0.0
        ident = sp.csr_matrix((np.ones(len(rows)), (rows, rows)), shape=J.shape)
        return sp.csc_matrix(sp.diags(mask) @ J + ident)

    def solve_bc(self, X0, bc: Dirichlet, max_iter=200, max_step=2.0) -> QFSolution:
        rows, vals = self._all_bc(bc)
        X = np.array(X0, dtype=float)
        X[rows] = vals
        free = np.ones(len(X), dtype=bool)
        free[rows] = False
        rep = QFReport()
        for _ in range(max_iter + 1):
            F = self.residual_bc(X, rows, vals)
            J = self.jacobian_bc(X, rows)
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
            try:
                dX = safe_spsolve(sp.csc_matrix(Js), -F / rs)
            except RuntimeError:
                break
            if not np.all(np.isfinite(dX)):
                break
            step = float(np.max(np.abs(dX[0::3][self.reg.semi_node])))  # limiter: semiconductor psi
            if step > max_step:
                dX *= max_step / step
                rep.limited_steps += 1
            rep.last_dX = float(np.max(np.abs(dX)))
            X += dX
            rep.iterations += 1
        psi, pn, pp = self.split(X)
        n, p = self.densities(X)
        semi = self.reg.semi_node
        n, p = np.where(semi, n, 0.0), np.where(semi, p, 0.0)
        jn, jp = self.fluxes(X)
        return QFSolution(X.copy(), psi.copy(), pn.copy(), pp.copy(), n, p, jn, jp, rep)

    # ------------------------------------------------------------------
    def fd_jacobian(self, X, bc: Dirichlet, h=1e-5):
        """Central-difference Jacobian of the BC-modified residual."""
        rows, vals = self._all_bc(bc)
        cols = []
        for k in range(len(X)):
            Xp, Xm = X.copy(), X.copy()
            Xp[k] += h
            Xm[k] -= h
            cols.append(
                (self.residual_bc(Xp, rows, vals) - self.residual_bc(Xm, rows, vals)) / (2 * h)
            )
        return sp.csc_matrix(np.array(cols).T)

    def space_charge(self, sol: QFSolution) -> float:
        """Integral of rho over the semiconductor part of the control volumes."""
        return float(np.sum(self.vol * (sol.p - sol.n + self.c)))


def continue_bias(
    model: MOSDriftDiffusion,
    X0,
    bc_of: Callable[[np.ndarray], Dirichlet],
    p_start,
    targets,
    equilibrium=False,
    min_frac=1 / 256,
    **kw,
):
    """March the bias parameter vector from ``p_start`` (X0 consistent with it) through
    ``targets`` with automatic step halving.

    Returns a list of (target, QFSolution, newton_iterations, solves); a target whose
    solve failed down to ``min_frac`` is returned with a non-converged solution and the
    march stops there.
    """
    p_cur = np.atleast_1d(np.asarray(p_start, dtype=float))
    X = np.array(X0, dtype=float)
    out = []

    def build(p):
        bc = bc_of(p)
        return model.equilibrium_bc(bc) if equilibrium else bc

    for tgt in targets:
        tgt = np.atleast_1d(np.asarray(tgt, dtype=float))
        done, step, its, solves = 0.0, 1.0, 0, 0
        p0 = p_cur.copy()
        sol = None
        while done < 1.0 - 1e-12:
            step = min(step, 1.0 - done)
            p_try = p0 + (done + step) * (tgt - p0)
            s = model.solve_bc(X, build(p_try), **kw)
            solves += 1
            its += s.report.iterations
            if s.report.converged:
                done += step
                X, sol, p_cur = s.X.copy(), s, p_try
                if s.report.iterations < 8:
                    step = min(1.0, 2 * step)
            else:
                step *= 0.5
                if step < min_frac:
                    out.append((tgt, s, its, solves))
                    return out
        if sol is None:
            sol = model.solve_bc(X, build(tgt), **kw)
            solves += 1
            its += sol.report.iterations
        out.append((tgt, sol, its, solves))
    return out
