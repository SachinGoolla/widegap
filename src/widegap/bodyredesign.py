r"""Stage 20: p-body depth family and the electrostatic punch-through voltage.

Spec: docs/benchmarks/stage20_body_redesign.md (pre-registered; AUTHORITATIVE).

* :func:`body_design` changes ONLY ``d_body`` of a :class:`widegap.sic.DMOSDesign`; the y-mesh
  segment of the body (``ny_body`` cells between ``-d_body`` and ``-d_src``) gets a cell count
  proportional to its thickness (``ceil(11 (d_body - d_src) / 0.5 um)``, equal to the stage-12
  value 11 at 0.8 um) so that the vertical cell size stays at the 0.8 um design's.  Everything
  else (surface and channel doping, N_A, d_src, drift thickness ``w_drift`` and doping, JFET
  width/doping, oxide) is untouched: the drift layer starts below the body junction
  (``y = -(d_body + w_drift)``), so the device gets deeper by the same amount.
* :func:`source_columns` / :func:`column_index_matrix`: mesh columns under the n+ source (surface
  node net n-type, ``x_src0 <= x <= x_bc``) and their window nodes (surface to ``d_body + 0.5 um``).
* :func:`barrier_min`: ``min over columns of psi_surface - min_window psi``.
* :class:`OffMarch`: zero-current off-state march (stage-12c model) with step halving and
  :meth:`vpt_es`, the smallest V_D at which the minimum barrier falls below a level.
"""

from __future__ import annotations

import dataclasses
import math

import numpy as np

from widegap.adjoint_sic import EvalError
from widegap.adjoint_sic_corner import _AffineBC

BARRIER_BELOW_BODY_UM = 0.5
STAGE12_NY_BODY = 11
STAGE12_BODY_THICKNESS_UM = 0.5  # d_body - d_src of the 0.8 um design


def body_design(design, d_body_um: float):
    """``design`` with only the body depth changed (and the body y-cell count scaled)."""
    n = math.ceil(STAGE12_NY_BODY * (d_body_um - design.d_src) / STAGE12_BODY_THICKNESS_UM - 1e-9)
    kw = dict(design.mesh_kw)
    kw["ny_body"] = max(int(n), 1)
    return dataclasses.replace(design, d_body=float(d_body_um), mesh_kw=kw)


def source_columns(ev):
    """Indices of mesh columns under the n+ source: surface node net n-type and
    x_src0 <= x <= x_bc (the p-body columns whose surface is p-type carry no n+ -> p barrier)."""
    mesh, p = ev.dev.mesh, ev.dev.prm
    xy = mesh.node_xy
    cols = []
    for i, x in enumerate(mesh.x):
        if x < p.x_src0 - 1e-9 or x > p.x_bc + 1e-9:
            continue
        top = np.flatnonzero((xy[:, 0] == x) & (xy[:, 1] == 0.0))[0]
        if ev.nd[top] > ev.na[top]:
            cols.append(i)
    return cols


def column_index_matrix(ev, cols):
    """(len(cols), n_window) node indices, surface first, window = y >= -(d_body + 0.5 um)."""
    mesh, p, sc = ev.dev.mesh, ev.dev.prm, ev.sc
    xy = mesh.node_xy
    y_lo = -(p.d_body + float(sc.length_s(BARRIER_BELOW_BODY_UM * 1e-4)))
    rows = []
    for i in cols:
        nodes = np.flatnonzero(
            (xy[:, 0] == mesh.x[i]) & (xy[:, 1] <= 0.0) & (xy[:, 1] >= y_lo - 1e-12)
        )
        rows.append(nodes[np.argsort(-xy[nodes, 1])])
    return np.array(rows)


def barrier_min(psi, idx):
    """(min barrier, column row of the minimum, all barriers) for ``psi`` over index matrix."""
    p = np.asarray(psi)[idx]
    bar = p[:, 0] - p.min(axis=1)
    k = int(np.argmin(bar))
    return float(bar[k]), k, bar


class OffMarch:
    """Stage-12c zero-current off-state marching with step halving (scaled voltages inside)."""

    def __init__(self, ev, theta):
        self.ev, self.sc = ev, ev.sc
        ev.set_theta(theta)
        eq = ev._eq()
        rows0, vals0 = ev._off_bc(0.0)
        s0 = ev.mod.solve_bc(eq.X, _AffineBC(rows0, vals0), max_step=20.0)
        if not s0.report.converged:
            raise EvalError("off-state zero-bias solve did not converge")
        self.v, self.x = 0.0, s0.X.copy()
        self.n_solves = 0

    def advance(self, state, dv_volt):
        """March from ``state = (v_s, X)`` by ``dv_volt`` [V]; returns the new (v_s, X)."""
        ev, sc = self.ev, self.sc
        v, x = state
        target = v + float(sc.volt_s(dv_volt))
        h = target - v
        while v < target - 1e-12:
            h_try = min(h, target - v)
            rows, dvals = ev._dvals(v)
            t = ev._tangent(x, rows, dvals)
            s = ev._off_solve(v + h_try, x + h_try * t)
            self.n_solves += 1
            if s is None:
                h = h_try * 0.5
                if h < 1e-3:
                    raise EvalError(f"off-state march failed at V_D = {sc.volt_phys(v):.1f} V")
                continue
            v, x = v + h_try, s.X.copy()
        return v, x

    def vpt_es(self, idx, level_vt=10.0, step_v=5.0, coarse_after_v=200.0, coarse_step_v=10.0,
               v_max_v=1200.0, tol_v=0.5):  # fmt: skip
        """Smallest V_D [V] with min-barrier < ``level_vt`` (bisection to ``tol_v``); None if not
        reached below ``v_max_v``.  Returns dict with the barrier curve."""
        sc = self.sc
        state = (self.v, self.x)
        b0, k0, _ = barrier_min(state[1][0::3], idx)
        curve = [(0.0, b0, k0)]
        prev = state
        while True:
            v_phys = float(sc.volt_phys(state[0]))
            if v_phys >= v_max_v:
                return {"V_pt_es": None, "curve": curve, "n_solves": self.n_solves}
            step = step_v if v_phys < coarse_after_v else coarse_step_v
            prev = state
            state = self.advance(state, step)
            b, k, _ = barrier_min(state[1][0::3], idx)
            curve.append((float(sc.volt_phys(state[0])), b, k))
            if b < level_vt:
                break
        lo, hi = prev, state
        b_lo = barrier_min(lo[1][0::3], idx)[0]
        while float(sc.volt_phys(hi[0] - lo[0])) > tol_v:
            dv = 0.5 * float(sc.volt_phys(hi[0] - lo[0]))
            mid = self.advance(lo, dv)
            bm = barrier_min(mid[1][0::3], idx)[0]
            if bm < level_vt:
                hi = mid
            else:
                lo, b_lo = mid, bm
        v_lo, v_hi = float(sc.volt_phys(lo[0])), float(sc.volt_phys(hi[0]))
        b_hi, k_hi, bars = barrier_min(hi[1][0::3], idx)
        # linear interpolation of the barrier between the final bracket ends
        f = (b_lo - level_vt) / (b_lo - b_hi) if b_lo != b_hi else 0.5
        return {
            "V_pt_es": float(v_lo + f * (v_hi - v_lo)),
            "V_bracket": [v_lo, v_hi],
            "worst_column_row": int(k_hi),
            "barriers_at_hi": bars.tolist(),
            "curve": curve,
            "n_solves": self.n_solves,
        }
