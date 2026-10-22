r"""Stage 13d: gate-oxide field and local SiC field measures at a fixed blocking voltage, with adjoints.

Spec: docs/benchmarks/stage13d_oxide_constrained.md (pre-registered; AUTHORITATIVE).  Built on the
stage-13c model (:class:`widegap.offstate.CorrectedCornerAdjoint`: rounded device, corrected
zero-current off state ``phi_n = V_D w``, ``phi_p = 0``); nothing existing is modified.

Field measures (all at the fixed drain voltage ``V_BLK_V`` = 400 V, gate and source at 0 V)
---------------------------------------------------------------------------------------------
A *smooth maximum* of edge field magnitudes (log-sum-exp with fixed ``kappa``)

    S_kappa(E) = (1/kappa) log sum_i exp(kappa E_i),      max E <= S_kappa <= max E + ln(N)/kappa,

of the physical edge fields ``E_i = c |psi_b - psi_a| / h_i`` [MV/cm] (``c`` converts the scaled
field, V_T / L_D, to MV/cm; the oxide field is the same potential gradient, in the oxide).

* ``E_ox``: the vertical edges of the first oxide layer, i.e. the oxide edges that touch the
  SiC/SiO2 interface (``y = 0``) at every mesh column (under the gate and over the contacts).
* ``W1``, ``W2``, ``W3``: the semiconductor edges (dual-face weight > 0, both directions) whose
  midpoint lies in a fixed physical box.  Recorded BEFORE any stage-13d run (from the stage-13c
  geometry and locations only: JFET/body junction corner (x_jfet, d_body) = (0.8, -0.8) um;
  global field peak at (1.69, -0.86) um under the source edge; JFET-centre surface):

  ``WINDOWS_UM = {W1: x in [0.5, 1.1], y in [-1.1, -0.5];  W2: x in [1.3, 2.0], y in [-1.2, -0.6];
                  W3: x in [0.0, 0.4], y in [-0.2, 0.0]}``

  (W2 stays below the rounded n+ source junction, y <= -0.6 um, whose built-in field is
  independent of V_D.)  ``KAPPA_PER_MV_CM = 20`` (smoothing length 0.05 MV/cm; bias at most
  ln(N)/20 ~ 0.15-0.2 MV/cm upward, the same for all designs).

Adjoint at fixed ``V_blk``
--------------------------
With the electrostatic residual ``F(X; nd) = 0`` (BC-modified, zero-current ``phi_n = V_D w`` with
the theta = 0 classification, so the Dirichlet rows do not depend on theta), a field measure
``S(X)`` that depends on theta only through ``X`` has ``dS/dnd = -lam^T dF/dnd`` with
``F_X^T lam = dS/dX`` (exactly the :meth:`SiCAdjoint.bv` machinery without the voltage root-find),
and ``dS/dtheta = (dnd/dtheta)^T dS/dnd``.  ``dS/dpsi`` is analytic (soft-max weights times
``c sign(dpsi)/h``).
"""

from __future__ import annotations

import numpy as np
from scipy.special import logsumexp

from widegap.adjoint2d import solve_transposed
from widegap.adjoint_sic import EvalError
from widegap.adjoint_sic_corner import solve_at_voltage
from widegap.offstate import CorrectedCornerAdjoint

V_BLK_V = 400.0
KAPPA_PER_MV_CM = 20.0
WINDOWS_UM = {
    "W1": ((0.5, 1.1), (-1.1, -0.5)),  # JFET/body corner (x_jfet, d_body) = (0.8, -0.8)
    "W2": ((1.3, 2.0), (-1.2, -0.6)),  # body junction under the source edge x_src0 = 1.6
    "W3": ((0.0, 0.4), (-0.2, 0.0)),  # JFET-centre surface under the gate oxide
}
MEASURES = ("E_ox", "W1", "W2", "W3")
_TOL = 1e-9


def smooth_max(e, kappa):
    """(S_kappa(e), softmax weights) with S = log sum exp(kappa e) / kappa."""
    z = kappa * np.asarray(e, dtype=float)
    lse = float(logsumexp(z))
    return lse / kappa, np.exp(z - lse)


class FieldMeasure:
    """Smooth maximum of ``c |psi_b - psi_a| / h`` over a fixed edge set."""

    def __init__(self, a, b, h, c_mv_cm, kappa, mid):
        self.a, self.b, self.h = np.asarray(a), np.asarray(b), np.asarray(h, dtype=float)
        self.c, self.kappa = float(c_mv_cm), float(kappa)
        self.mid = np.asarray(mid, dtype=float)  # (n_edges, 2) midpoints (scaled units)
        if len(self.a) == 0:
            raise ValueError("empty edge set")

    def edge_field(self, psi):
        return self.c * np.abs(psi[self.b] - psi[self.a]) / self.h

    def value(self, psi):
        return smooth_max(self.edge_field(psi), self.kappa)[0]

    def value_grad(self, psi):
        """(S, dS/dpsi over all nodes)."""
        d = psi[self.b] - psi[self.a]
        e = self.c * np.abs(d) / self.h
        s, w = smooth_max(e, self.kappa)
        coef = w * self.c * np.sign(d) / self.h
        g = np.bincount(self.b, coef, minlength=len(psi)) - np.bincount(
            self.a, coef, minlength=len(psi)
        )
        return s, g

    def hard_max(self, psi):
        return float(np.max(self.edge_field(psi)))


class OxideAdjoint(CorrectedCornerAdjoint):
    """:class:`CorrectedCornerAdjoint` plus the fixed-voltage field measures of stage 13d."""

    def __init__(self, sc, design, level=0, r_um=0.1, cuts_um=(), tau_v=5.0, v_blk_v=V_BLK_V,
                 kappa=KAPPA_PER_MV_CM, windows_um=None, **kw):  # fmt: skip
        super().__init__(sc, design, level, r_um=r_um, cuts_um=cuts_um, tau_v=tau_v, **kw)
        self.v_blk_v, self.kappa = float(v_blk_v), float(kappa)
        self.windows_um = dict(WINDOWS_UM if windows_um is None else windows_um)
        mesh, mod = self.mesh, self.mod
        a, b, h = mod.a, mod.b, mesh.h
        xy = mesh.node_xy
        mid = 0.5 * (xy[a] + xy[b])
        c = float(sc.e_unit) * 1e-6  # scaled field -> MV/cm
        self._c = c
        t = self.prm.t_ox
        ya, yb = xy[a, 1], xy[b, 1]
        ox = (mesh.edge_dir == 1) & (np.minimum(ya, yb) > -_TOL) & (np.minimum(ya, yb) < _TOL) & (
            np.maximum(ya, yb) > _TOL)  # fmt: skip
        if t <= 0 or not ox.any():
            raise ValueError("no oxide edges adjacent to the interface")
        self.measures: dict[str, FieldMeasure] = {
            "E_ox": FieldMeasure(a[ox], b[ox], h[ox], c, kappa, mid[ox])
        }
        semi = self.dev.reg.w_s > 0
        um = 1e-4 / float(sc.l_d)  # um -> scaled length
        self.window_masks = {}
        for name, ((x0, x1), (y0, y1)) in self.windows_um.items():
            sel = (semi & (mid[:, 0] >= x0 * um - _TOL) & (mid[:, 0] <= x1 * um + _TOL)
                   & (mid[:, 1] >= y0 * um - _TOL) & (mid[:, 1] <= y1 * um + _TOL))  # fmt: skip
            self.window_masks[name] = sel
            self.measures[name] = FieldMeasure(a[sel], b[sel], h[sel], c, kappa, mid[sel])
        self._x_blk = None

    # ----------------------------------------------------------------------------- state
    def state_blk(self, theta, v_volt=None):
        """Converged electrostatic state X at the blocking voltage (warm start from the previous
        state of this evaluator, else the 20 V march from the equilibrium)."""
        v_volt = self.v_blk_v if v_volt is None else float(v_volt)
        self.set_theta(theta)
        v_s = float(self.sc.volt_s(v_volt))
        if self._x_blk is not None and abs(self._x_blk[0] - v_s) < 1e-12:
            s = self._off_solve(v_s, self._x_blk[1])
            if s is not None:
                self._x_blk = (v_s, s.X.copy())
                return s.X
        X = solve_at_voltage(self, theta, v_volt)
        self._x_blk = (v_s, X.copy())
        return X

    def fields(self, theta, grad=True, which=MEASURES, v_volt=None):
        """Smooth-max field measures [MV/cm] at ``V_blk`` (and dS/dtheta, dS/dnd) for ``which``."""
        v_volt = self.v_blk_v if v_volt is None else float(v_volt)
        try:
            X = self.state_blk(theta, v_volt)
        except EvalError:
            self.stats["failed_evals"] += 1
            raise
        psi = X[0::3]
        out = {"X": X, "v_volt": v_volt}
        rows = None
        if grad:
            rows, _ = self._off_bc(float(self.sc.volt_s(v_volt)))
            jac = self.mod.jacobian_bc(X, rows)
        for name in which:
            m = self.measures[name]
            if not grad:
                out[name] = {"S": m.value(psi), "hard_max": m.hard_max(psi)}
                continue
            s, gpsi = m.value_grad(psi)
            g = np.zeros(len(X))
            g[0::3] = gpsi
            lam = solve_transposed(jac, g)
            lam[rows] = 0.0
            g_nd = self.dF_dnd(-lam, X)
            out[name] = {"S": s, "hard_max": m.hard_max(psi), "g_nd": g_nd,
                         "grad": self.dnd_dtheta.T @ g_nd}  # fmt: skip
        return out

    def edge_fields(self, X):
        """Per-edge physical field [MV/cm] of every semiconductor edge (dual weight > 0) and of
        the interface-adjacent oxide edges, with midpoints in um (diagnostics only)."""
        psi = X[0::3]
        l_um = float(self.sc.l_d) * 1e4
        mod, mesh = self.mod, self.mesh
        e_all = self._c * np.abs(psi[mod.b] - psi[mod.a]) / mesh.h
        mid = 0.5 * (mesh.node_xy[mod.a] + mesh.node_xy[mod.b]) * l_um
        semi = np.flatnonzero(self.dev.reg.w_s > 0)
        mo = self.measures["E_ox"]
        return {"semi_E": e_all[semi], "semi_x_um": mid[semi, 0], "semi_y_um": mid[semi, 1],
                "semi_dir": mesh.edge_dir[semi], "ox_E": mo.edge_field(psi),
                "ox_x_um": mo.mid[:, 0] * l_um}  # fmt: skip

    def window_info(self):
        return {k: {"box_um": [list(map(float, v[0])), list(map(float, v[1]))],
                    "n_edges": int(self.window_masks[k].sum())} for k, v in self.windows_um.items()} | {
            "E_ox": {"n_edges": len(self.measures["E_ox"].a)}}  # fmt: skip


__all__ = [
    "KAPPA_PER_MV_CM",
    "MEASURES",
    "V_BLK_V",
    "WINDOWS_UM",
    "FieldMeasure",
    "OxideAdjoint",
    "smooth_max",
]
