r"""Stage 19: 1D drift-diffusion punch-through check of the study DMOSFET.

Spec: docs/benchmarks/stage19_punchthrough_check.md (pre-registered; AUTHORITATIVE).

The n+ source / p-body / n- drift / n+ drain stack along one vertical line of the stage-13d
device (rounded junctions, theta = 0) is solved with the stage-9a edge core
(:class:`widegap.sic.AnisoMobilityDD`, a stage-8 one-column mesh = the stage-7.1 quasi-Fermi
scheme in 1D) with incomplete ionization, Caughey-Thomas velocity saturation and the device's
doping-dependent vertical electron mobility, no impact ionization, no recombination.  Ohmic
contacts at the top (source, grounded) and bottom (drain, V_D); the p-body has no contact
(``body="floating"``) or a hole-only contact ``phi_p = 0`` at one body node (``body="tied"``).

Doping.  The cell-centre doping of the 2D device is evaluated pointwise at the column's x
(the exact functions of :class:`widegap.adjoint_sic_corner.RoundedDMOSFET`, the same cell
membership masks) and node-averaged over the 1D control volume, so both mesh levels discretize
the same continuous profile.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from widegap.adjoint_sic_corner import smooth_step
from widegap.fv2d.mesh import TensorMesh
from widegap.fv2d.mos import Dirichlet, continue_bias
from widegap.fv2d.physics import Ionization, neutral_psi
from widegap.fv2d.regions import Regions
from widegap.linalg import safe_spsolve
from widegap.sic import AnisoMobilityDD

BARRIER_BELOW_BODY_UM = 0.5  # barrier window: surface to d_body + 0.5 um (stage-17c scoping)
Q_E = 1.602176634e-19


def column_x(ev, x_um):
    """Mesh column nearest to ``x_um`` of the device ``ev`` (a stage-12c/13d evaluator)."""
    sc, mesh = ev.sc, ev.dev.mesh
    xs = mesh.x
    i = int(np.argmin(np.abs(xs - float(sc.length_s(x_um * 1e-4)))))
    return i, float(xs[i])


def column_doping(ev, x_s, theta):
    """(t, nd, na, cell_nd, cell_na) along the vertical line x = ``x_s`` (scaled units).

    ``t = -y`` ascending from the surface; node values are 1D control-volume averages of the
    cell-centre values.  Cell membership follows ``DMOSFET`` (drain / drift / upper layer)."""
    dev = ev.dev
    p = dev.prm
    yv = dev.mesh.y
    ys = np.sort(yv[yv <= 1e-12])[::-1]  # 0, -dy, ... (top first)
    t = -ys
    yc = 0.5 * (ys[:-1] + ys[1:])
    st = lambda s: smooth_step(s, dev.r)  # noqa: E731
    y_drift0 = -(p.d_body + p.w_drift)
    th = np.asarray(theta, dtype=float)
    nd = np.zeros_like(yc)
    na = np.zeros_like(yc)
    m_drain = yc < y_drift0
    m_drift = (yc >= y_drift0) & (yc < -p.d_body)
    m_upper = yc >= -p.d_body
    nd[m_drain] = p.b_drain
    nd[m_drift] = dev.drift_profile(th, yc[m_drift])
    nd += np.where(m_upper, p.n_jfet * np.exp(th[4]) * st(p.x_jfet - x_s), 0.0)
    nd += np.where(
        m_upper, p.b_src * st(x_s - p.x_src0) * st(p.x_src1 - x_s) * st(yc + p.d_src), 0.0
    )
    na += np.where(m_upper, p.a_body * st(x_s - p.x_jfet) * st(yc + p.d_body), 0.0)
    hc = np.diff(t)
    vol = np.zeros(len(t))
    vol[:-1] += 0.5 * hc
    vol[1:] += 0.5 * hc

    def cv(q):
        out = np.zeros(len(t))
        out[:-1] += 0.5 * hc * q
        out[1:] += 0.5 * hc * q
        return out / vol

    return t, cv(nd), cv(na), nd, na


class Column1D:
    """One vertical line of the study device as a 1D drift-diffusion problem."""

    def __init__(self, ev, x_um, theta, body="floating", body_depth_um=None):
        if body not in ("floating", "tied"):
            raise ValueError("body must be 'floating' or 'tied'")
        self.ev, self.sc = ev, ev.sc
        sc, base, p = self.sc, ev.base, ev.dev.prm
        self.i_col, self.x_s = column_x(ev, x_um)
        self.x_um = float(sc.length_cm(self.x_s)) * 1e4
        self.t, self.nd, self.na, self.cell_nd, self.cell_na = column_doping(ev, self.x_s, theta)
        self.body = body
        n = len(self.t)
        mesh = TensorMesh(self.t)
        reg = Regions.from_function(mesh, lambda xc: np.ones_like(xc, dtype=bool))
        nd_cm = sc.n0 * 0.5 * (self.nd[:-1] + self.nd[1:] + self.na[:-1] + self.na[1:])
        mu_edge = sc.mu_s(sc.mobility_cm2(nd_cm, "n", "par"))
        self.mu_edge = np.asarray(mu_edge, dtype=float)
        self.donors, self.acceptors = base.donors, base.acceptors
        self.model = AnisoMobilityDD(
            reg, self.nd - self.na, base.dp, float(np.log(sc.delta)), eps_ratio=1.0,
            ionization=Ionization(self.nd, self.na, self.donors, self.acceptors),
            mobility=base.mob, mu_edge=self.mu_edge,
        )  # fmt: skip
        self.mesh, self.n = mesh, n
        self.psi_c = np.array(
            [
                neutral_psi(self.nd[k], self.na[k], sc.delta, self.donors, self.acceptors)
                for k in (0, n - 1)
            ]
        )
        d_um = lambda v: float(sc.length_cm(v)) * 1e4  # noqa: E731
        self.d_body_um, self.d_src_um = d_um(p.d_body), d_um(p.d_src)
        self.t_um = np.array([d_um(v) for v in self.t])
        self.win = self.t_um <= self.d_body_um + BARRIER_BELOW_BODY_UM + 1e-9
        bd = body_depth_um if body_depth_um is not None else 0.5 * (self.d_src_um + self.d_body_um)
        self.i_body = int(np.argmin(np.abs(self.t_um - bd)))
        self.body_node_um = float(self.t_um[self.i_body])
        self.stats = {"solves": 0, "failed": 0}

    # ----- boundary conditions ----------------------------------------------------
    def bc(self, v_d_s):
        bc = Dirichlet()
        bc.psi([0], self.psi_c[0]).phi([0], 0.0)
        bc.psi([self.n - 1], self.psi_c[1] + v_d_s).phi([self.n - 1], v_d_s)
        if self.body == "tied":
            bc._set(3 * self.i_body + 2, 0.0)  # hole-only contact: phi_p = 0
        return bc

    def equilibrium(self):
        mod = self.model
        bc0 = self.bc(0.0)
        eq = mod.equilibrium_bc(bc0)
        return mod.solve_bc(mod.initial_guess(bc0), eq, max_step=20.0)

    def march(self, X0, v_start_volt, targets_volt):
        """Continuation through ``targets_volt`` from ``X0`` (consistent with ``v_start_volt``)."""
        sc = self.sc
        tg = [[float(sc.volt_s(v))] for v in targets_volt]
        res = continue_bias(
            self.model, X0, lambda q: self.bc(q[0]), [float(sc.volt_s(v_start_volt))], tg,
            max_step=20.0,
        )  # fmt: skip
        out = []
        for v, (_, sol, its, solves) in zip(targets_volt, res, strict=False):
            self.stats["solves"] += solves
            self.stats["failed"] += int(not sol.report.converged)
            out.append((v, sol, its, solves))
        return out

    # ----- observables -------------------------------------------------------------
    def j_scale(self):
        """Scaled 1D current density -> A/cm^2: q N_0 mu_0 V_T / L_D."""
        sc = self.sc
        return Q_E * sc.n0 * sc.mu0 * sc.vt / sc.l_d

    def observe(self, sol):
        """Observables of one converged state.

        ``Jn_A_cm2`` is the electron flux on the edge just below the potential minimum of the
        barrier window.  In neutral n material phi_n differs from V_D by less than round-off
        for the currents of interest, so edge fluxes there are not resolved (their spread is
        recorded, not used); the barrier edges carry the current to relative 1e-14 (the
        per-row tolerance is relative to the row's own tiny coefficients)."""
        sc = self.sc
        psi = sol.psi
        js = self.j_scale()
        jn = np.asarray(sol.jn, dtype=float) * js
        jp = np.asarray(sol.jp, dtype=float) * js
        iw = np.flatnonzero(self.win)
        k = int(iw[np.argmin(psi[iw])])
        bar = float(psi[0] - psi[k])
        jb = float(jn[min(k, len(jn) - 1)])
        ew = jn[: iw[-1]]  # edges from the surface to the end of the barrier window
        return {
            "Jn_A_cm2": jb,
            "Jn_edge_spread_rel": float((ew.max() - ew.min()) / max(abs(jb), 1e-300)),
            "Jp_A_cm2": float(jp[min(k, len(jp) - 1)]),
            "barrier_VT": bar,
            "psi_min_t_um": float(self.t_um[k]),
            "phi_p_body_V": float(sol.phi_p[self.i_body] * sc.vt),
            "phi_n_body_V": float(sol.phi_n[self.i_body] * sc.vt),
            "psi_body_V": float(sol.psi[self.i_body] * sc.vt),
        }


# ------------------------------------------------------------------------- analysis helpers
def first_crossing(v, y, level, above=True):
    """Smallest V at which ``y(V)`` first crosses ``level`` (above: y > level; else y <= level),
    by linear interpolation of ``y`` between the two bracketing grid points; None if never."""
    v, y = np.asarray(v, dtype=float), np.asarray(y, dtype=float)
    hit = np.flatnonzero((y > level) if above else (y <= level))
    if len(hit) == 0:
        return None
    k = int(hit[0])
    if k == 0:
        return float(v[0])
    f = (level - y[k - 1]) / (y[k] - y[k - 1])
    return float(v[k - 1] + f * (v[k] - v[k - 1]))


def interp_at(v, y, v0):
    """Linear interpolation of ``y(V)`` at ``v0`` (None outside the grid)."""
    v = np.asarray(v, dtype=float)
    if v0 < v[0] or v0 > v[-1]:
        return None
    return float(np.interp(v0, v, np.asarray(y, dtype=float)))


def punch_voltage(volts, jn_abs, v_ref=10.0, factor=1.0e3):
    """V_pt: smallest V_D at which |J_n| exceeds ``factor`` |J_n(v_ref)| (log-linear interpolation
    between grid points).  Returns dict(V_pt, threshold, J_ref); V_pt None if not reached."""
    volts, jn_abs = np.asarray(volts, dtype=float), np.asarray(jn_abs, dtype=float)
    j_ref = interp_at(volts, jn_abs, v_ref)
    if j_ref is None or j_ref <= 0:
        return {"V_pt": None, "threshold": None, "J_ref": j_ref}
    ln = np.log(np.maximum(jn_abs, 1e-300))
    thr = float(np.log(factor * j_ref))
    mask = volts > v_ref
    vp = first_crossing(volts[mask], ln[mask], thr, above=True)
    return {"V_pt": vp, "threshold": float(np.exp(thr)), "J_ref": float(j_ref)}


def zero_current_barrier(col: Column1D, theta, volts):
    """Stage-12c zero-current barrier (psi_surface - min psi over the same window) on the same
    column of the 2D device, from the existing off-state march ``solve_at_voltage``."""
    from widegap.adjoint_sic_corner import solve_at_voltage

    ev = col.ev
    mesh = ev.dev.mesh
    xy = mesh.node_xy
    nodes = np.flatnonzero((xy[:, 0] == mesh.x[col.i_col]) & (xy[:, 1] <= 0.0))
    nodes = nodes[np.argsort(-xy[nodes, 1])]
    if len(nodes) != col.n:
        raise RuntimeError("1D and 2D columns differ")
    out = {}
    for v in volts:
        psi = solve_at_voltage(ev, theta, float(v))[0::3][nodes]
        out[float(v)] = float(psi[0] - psi[col.win].min())
    return out


def failure_record(col: Column1D, sol, v_target_volt):
    """Where and how a non-converged solve fails: which of the two Newton tests (residual
    ``|F|/S <= 1e2 eps``, last step ``<= 1e-8``) is violated, the worst residual row, and the
    variable/node of the largest Newton step.  ``sol`` is the last (sub-)step returned by
    ``continue_bias``; free-row residuals do not depend on the boundary values."""
    rep = sol.report
    eps_row = 1e2 * np.finfo(float).eps
    res = float(rep.residual_history[-1]) if rep.residual_history else float("nan")
    rec = {
        "target_V": float(v_target_volt),
        "iterations": int(rep.iterations),
        "last_residual_ratio": res,
        "last_dX": float(rep.last_dX),
        "residual_test_ok": bool(res <= eps_row),
        "step_test_ok": bool(rep.last_dX <= 1e-8),
        "state_finite": bool(np.all(np.isfinite(sol.X))),
    }
    if rec["state_finite"]:
        mod = col.model
        bc = col.bc(float(col.sc.volt_s(v_target_volt)))
        rows, vals = mod._all_bc(bc)
        free = np.ones(len(sol.X), dtype=bool)
        free[rows] = False
        F = mod.residual_bc(sol.X, rows, vals)
        J = mod.jacobian_bc(sol.X, rows)
        S = mod.row_scale(sol.X, J)
        ratio = np.where(free, np.abs(F) / S, 0.0)
        i = int(np.argmax(ratio))
        n, p = mod.densities(sol.X)
        rs = np.asarray(abs(J).max(axis=1).todense()).ravel()
        dX = safe_spsolve(sp.csc_matrix(sp.diags(1.0 / rs) @ J), -F / rs)
        if np.all(np.isfinite(dX)):
            comp = {}
            for k, nm in enumerate(("psi", "phi_n", "phi_p")):
                a = np.abs(np.where(free, dX, 0.0)[k::3])  # Dirichlet rows excluded
                j = int(np.argmax(a))
                comp[nm] = {"max_abs_step": float(a[j]), "node": j, "depth_um": float(col.t_um[j]),
                            "p_at_node": float(p[j]), "n_at_node": float(n[j])}  # fmt: skip
            rec["newton_step_by_variable"] = comp
        rec["worst_row"] = {
            "variable": ("psi", "phi_n", "phi_p")[i % 3],
            "node": i // 3,
            "depth_um": float(col.t_um[i // 3]),
            "ratio": float(ratio[i]),
            "F": float(F[i]),
            "S": float(S[i]),
            "p_at_node": float(p[i // 3]),
            "n_at_node": float(n[i // 3]),
        }
    return rec
