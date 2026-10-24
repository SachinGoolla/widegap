r"""Stage 24: smooth n+ -> body barrier at a fixed blocking voltage, with its adjoint.

Spec: docs/benchmarks/stage24_guard_ablation.md (pre-registered; AUTHORITATIVE).

Per stage-20 barrier column c (mesh column under the n+ source, frozen at theta = 0) with window
nodes y_0 = surface, ..., the hard barrier is ``psi_c,0 - min_y psi_c,y`` (:func:`widegap.bodyredesign
.barrier_min`).  The smooth barrier nests two log-sum-exp soft-minima with temperature ``tau``
(scaled potential units = V_T):

    m_c = -tau log sum_y exp(-psi_c,y / tau),      b_c = psi_c,0 - m_c,
    B   = -tau log sum_c exp(-b_c / tau).

``B`` depends on the electrostatic potential only, so at fixed drain voltage the gradient is the
off-state adjoint of :class:`widegap.adjoint_sic_oxide.OxideAdjoint` (same structure as ``E_ox`` at
400 V): ``F_X^T lam = dB/dX`` (BC-modified Jacobian, Dirichlet rows do not depend on theta),
``dB/dnd = -lam^T dF/dnd``, ``dB/dtheta = (dnd/dtheta)^T dB/dnd``.  Additive module: nothing existing
is modified.
"""

from __future__ import annotations

import numpy as np
from scipy.special import logsumexp

from widegap.adjoint2d import solve_transposed
from widegap.adjoint_sic import EvalError

TAU_VT = 0.5
LEVEL_VT = 10.0
V_STAR_V = 650.0


class SoftBarrier:
    """Nested soft-min barrier of a (columns x window) node-index matrix."""

    def __init__(self, idx, tau=TAU_VT):
        self.idx = np.asarray(idx)
        self.tau = float(tau)

    def value_grad(self, psi):
        """(B, dB/dpsi over all nodes, info) with info = hard barrier and per-column b_c."""
        psi = np.asarray(psi, dtype=float)
        p = psi[self.idx]  # (n_cols, n_win)
        tau = self.tau
        lse_y = logsumexp(-p / tau, axis=1)
        m = -tau * lse_y
        w = np.exp(-p / tau - lse_y[:, None])  # d m / d p
        b = p[:, 0] - m
        lse_c = logsumexp(-b / tau)
        big = -tau * float(lse_c)
        u = np.exp(-b / tau - lse_c)  # dB/db_c
        dp = -u[:, None] * w  # dB/dp_cy through m_c: db_c/dp_cy = delta_y0 - w_cy
        dp[:, 0] += u
        g = np.zeros(len(psi))
        np.add.at(g, self.idx.ravel(), dp.ravel())
        hard = float(np.min(p[:, 0] - p.min(axis=1)))
        return float(big), g, {"hard": hard, "b_cols": b, "m_cols": m}

    def value(self, psi):
        return self.value_grad(psi)[0]


class BarrierAdjoint:
    """B_tau(theta; V*) and its gradient on an :class:`OxideAdjoint`-like evaluator ``ev``.

    ``ev`` must be dedicated to this measure (its cached off-state is at V*).  ``idx`` is the
    window node-index matrix (:func:`widegap.bodyredesign.column_index_matrix`, theta = 0)."""

    def __init__(self, ev, idx, v_star_v=V_STAR_V, tau=TAU_VT):
        self.ev, self.v_star, self.measure = ev, float(v_star_v), SoftBarrier(idx, tau)

    def evaluate(self, theta, grad=True):
        ev = self.ev
        try:
            X = ev.state_blk(theta, self.v_star)
        except EvalError:
            ev.stats["failed_evals"] += 1
            raise
        psi = X[0::3]
        b, gpsi, info = self.measure.value_grad(psi)
        out = {"B": b, "hard": info["hard"], "b_cols": info["b_cols"], "X": X}
        if grad:
            rows, _ = ev._off_bc(float(ev.sc.volt_s(self.v_star)))
            jac = ev.mod.jacobian_bc(X, rows)
            g = np.zeros(len(X))
            g[0::3] = gpsi
            lam = solve_transposed(jac, g)
            lam[rows] = 0.0
            g_nd = ev.dF_dnd(-lam, X)
            out["g_nd"] = g_nd
            out["grad"] = ev.dnd_dtheta.T @ g_nd
        return out


__all__ = ["BarrierAdjoint", "SoftBarrier", "LEVEL_VT", "TAU_VT", "V_STAR_V"]
