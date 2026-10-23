r"""Stage 15d: local-field mesh adequacy (diagnosis and the replacement area-mean measure).

Spec: docs/benchmarks/stage15_caveat_closure.md, section 15d (pre-registered; AUTHORITATIVE).
Everything here is post-processing of a converged 400 V electrostatic state ``psi`` of the stage-13d
evaluator (:class:`widegap.adjoint_sic_oxide.OxideAdjoint`); nothing existing is modified.

Definitions (frozen before the first 15d run; none is adjusted after seeing numbers)
------------------------------------------------------------------------------------
Hard-max position
    The 13d measures E_ox, W1, W2, W3 are smooth maxima of per-edge fields
    ``E_e = c |psi_b - psi_a| / h_e``; their hard maximum is the largest ``E_e`` over the
    measure's fixed edge set.  The *position* of a hard maximum is the midpoint of that edge
    (um, physical), together with the edge direction.

Fixed features (um; geometry of the device with ``theta = 0``, mesh independent)
    * ``gate_edge``         : the point (x_gate_end, t_ox), where the gate electrode ends on top
      of the oxide (``x_gate_end = x_src0``, the gate edge of the 13d device);
    * ``interface_corner``  : the planar SiC/SiO2 interface has no geometric corner.  It is
      defined here as the points of the interface ``y = 0`` where the boundary condition at the
      surface changes: (x_gate_end, 0) (end of the gate over the interface), the start of the
      source contact (x_src0 + gap, 0) and the start of the body contact (x_bc, 0).  FLAGGED:
      this is a definition made for the diagnosis, the spec names the feature but does not
      define it;
    * ``n_source_junction`` / ``p_body_junction`` : the ``N_D = N_A`` contour of the analytic
      (erfc-rounded, ``theta = 0``) doping of :class:`widegap.adjoint_sic_corner.RoundedDMOSFET`,
      sampled on a fixed 5 nm grid.  A contour point is assigned to the n+ source junction if
      the source donors exceed all other donors there, else to the p-body junction;
    * ``window_boundary``   : the boundary of the measure's own box (``WINDOWS_UM``); E_ox has no
      window.

Nearness
    A hard-max position is *within 2 cells* of a point/contour feature if its Euclidean distance
    to the feature is at most ``NEAR_CELLS`` times the larger side of the primal cell containing
    it (cell of that mesh).  It is *on the window boundary* if its distance to the box boundary is at
    most ``ON_BOUNDARY_TOL_UM``.  A measure is **pinned** if on all three meshes some feature
    satisfies this and the set of such features common to the three meshes is non-empty; the common
    features are the pinning features (R3).

Replacement measure (area-mean of |E|)
    ``Ebar`` is the area-weighted mean of ``|E|`` over the square of side ``SIDE_UM`` centred at a
    fixed point, intersected with the SiC (``y < 0``) or the oxide (``y > 0``).  Edge fields are
    assigned to their dual-face rectangles; each primal cell splits into four quarters, each of
    which lies in exactly one x-edge dual rectangle (field ``Ex``) and one y-edge dual rectangle
    (field ``Ey``), so ``|E| = sqrt(Ex^2 + Ey^2)`` is constant on a quarter.  Weights are the exact
    areas of the intersection of each quarter with the square (no midpoint-in-square test), and
    the mean divides by the intersected area actually present (reported as a fraction of L^2).
    ``y = 0`` is a mesh line, so no quarter straddles the interface.

L1 assessment (:func:`l1_assess`)
    Three values on h, h/2, h/4: ``e1 = v1 - v0``, ``e2 = v2 - v1``.  Monotone convergence is
    required (``e1 e2 > 0``); observed order ``p = log2(|e1/e2|)``; Richardson value
    ``v2 + e2/(2^p - 1)``.  Pass iff ``p >= ORDER_MIN`` (no upper bound) and the finest value is
    within ``RICH_TOL`` of the Richardson value.  A non-monotone sequence fails, with the reason.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree
from scipy.special import erfc

SIDE_UM = 0.1  # L, the rounding length r
NEAR_CELLS = 2.0
ON_BOUNDARY_TOL_UM = 1e-6
ORDER_MIN = 0.8
RICH_TOL = 0.02
CONTOUR_STEP_UM = 0.005
POINT_FEATURES = ("gate_edge", "interface_corner", "n_source_junction", "p_body_junction")
FEATURES = (*POINT_FEATURES, "window_boundary")


def smooth_step(s, r):
    """F(s) = erfc(-s / (sqrt(2) r)) / 2 (as :func:`widegap.adjoint_sic_corner.smooth_step`)."""
    return 0.5 * erfc(-np.asarray(s, dtype=float) / (np.sqrt(2.0) * r))


# ------------------------------------------------------------------------------- geometry
@dataclass(frozen=True)
class Geometry:
    """Device layout in um and cm^-3 (``theta = 0``)."""

    x_jfet: float
    x_src0: float
    x_src1: float
    x_bc: float
    x_max: float
    gate_x1: float
    src_contact_x0: float
    d_src: float
    d_body: float
    t_ox: float
    r_um: float
    n_drift: float
    n_jfet: float
    n_body: float
    n_src: float

    @classmethod
    def from_device(cls, sc, design, prm, r_um):
        """Geometry from the scaled :class:`widegap.devices.dmosfet.DMOSParams` ``prm`` of a built
        device (physical lengths via ``sc.l_d``) and the physical :class:`DMOSDesign` dopings."""
        um = float(sc.l_d) * 1e4
        return cls(
            x_jfet=prm.x_jfet * um, x_src0=prm.x_src0 * um, x_src1=prm.x_src1 * um,
            x_bc=prm.x_bc * um, x_max=prm.x_max * um, gate_x1=prm.gate_x1 * um,
            src_contact_x0=(prm.x_src0 + prm.src_contact_gap) * um,
            d_src=prm.d_src * um, d_body=prm.d_body * um, t_ox=prm.t_ox * um, r_um=float(r_um),
            n_drift=design.n_drift, n_jfet=design.n_jfet, n_body=design.n_body,
            n_src=design.n_src,
        )  # fmt: skip

    def doping(self, x, y):
        """(N_src, N_other_donors, N_A) at (x, y) [um], the ``theta = 0`` rounded doping."""
        r = self.r_um
        x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
        na = self.n_body * smooth_step(x - self.x_jfet, r) * smooth_step(y + self.d_body, r)
        src = (
            self.n_src
            * smooth_step(x - self.x_src0, r)
            * smooth_step(self.x_src1 - x, r)
            * smooth_step(y + self.d_src, r)
        )
        upper = y >= -self.d_body
        other = np.where(upper, self.n_jfet * smooth_step(self.x_jfet - x, r), self.n_drift)
        return src, other, na

    def contours(self, step=CONTOUR_STEP_UM):
        """Points (um) of the N_D = N_A contour, split into (n_source, p_body) arrays (n, 2).

        Sign changes of ``N_D,total - N_A`` between neighbouring points of a fixed grid are
        located by linear interpolation along the grid edge."""
        xs = np.arange(0.0, self.x_max + 0.5 * step, step)
        ys = np.arange(-(self.d_body + 0.6), 0.5 * step, step)
        xx, yy = np.meshgrid(xs, ys, indexing="ij")
        src, other, na = self.doping(xx, yy)
        s = src + other - na
        pts, from_src = [], []
        for ax in (0, 1):
            a = s[:-1, :] if ax == 0 else s[:, :-1]
            b = s[1:, :] if ax == 0 else s[:, 1:]
            cross = (a * b) < 0
            i, j = np.nonzero(cross)
            t = a[i, j] / (a[i, j] - b[i, j])
            x0, y0 = xs[i], ys[j]
            x1 = x0 + t * step if ax == 0 else x0
            y1 = y0 if ax == 0 else y0 + t * step
            pts.append(np.column_stack([x1, y1]))
            sr, ot, _ = self.doping(x1, y1)
            from_src.append(sr > ot)
        pts, from_src = np.vstack(pts), np.concatenate(from_src)
        return pts[from_src], pts[~from_src]


class Features:
    """Distances from a point to the fixed features of a :class:`Geometry`."""

    def __init__(self, geo: Geometry, step=CONTOUR_STEP_UM):
        self.geo = geo
        self.src_pts, self.body_pts = geo.contours(step)
        if len(self.src_pts) == 0 or len(self.body_pts) == 0:
            raise ValueError("junction contour not found")
        self._tree = {"n_source_junction": cKDTree(self.src_pts),
                      "p_body_junction": cKDTree(self.body_pts)}  # fmt: skip
        g = geo
        self.corners = np.array([[g.gate_x1, 0.0], [g.src_contact_x0, 0.0], [g.x_bc, 0.0]])
        self.gate = np.array([g.gate_x1, g.t_ox])

    def distances(self, x, y, window=None):
        """{feature: distance [um]} (window_boundary only if ``window`` = ((x0,x1),(y0,y1)))."""
        p = np.array([x, y], dtype=float)
        d = {"gate_edge": float(np.hypot(*(p - self.gate))),
             "interface_corner": float(np.min(np.hypot(*(self.corners - p).T)))}  # fmt: skip
        for k, t in self._tree.items():
            d[k] = float(t.query(p)[0])
        if window is not None:
            d["window_boundary"] = window_boundary_distance(x, y, window)
        return d


def window_boundary_distance(x, y, window):
    """Distance from (x, y) to the boundary of the rectangle ``((x0,x1),(y0,y1))`` (inside or out)."""
    (x0, x1), (y0, y1) = window
    dx = max(x0 - x, 0.0, x - x1)
    dy = max(y0 - y, 0.0, y - y1)
    if dx > 0.0 or dy > 0.0:  # outside: distance to the rectangle
        return float(np.hypot(dx, dy))
    return float(min(x - x0, x1 - x, y - y0, y1 - y))


# ------------------------------------------------------------------------- hard-max position
def cell_size_um(mesh, x_s, y_s, um_per_s):
    """Larger side [um] of the primal cell of the tensor mesh containing the point (scaled)."""
    i = int(np.clip(np.searchsorted(mesh.x, x_s, side="right") - 1, 0, mesh.nx - 1))
    j = int(np.clip(np.searchsorted(mesh.y, y_s, side="right") - 1, 0, mesh.ny - 1))
    return float(max(mesh.x[i + 1] - mesh.x[i], mesh.y[j + 1] - mesh.y[j]) * um_per_s)


def hard_max_record(measure, psi, mesh, um_per_s, features: Features, window=None):
    """Hard maximum of a :class:`widegap.adjoint_sic_oxide.FieldMeasure`, its position
    (edge midpoint, um), distances to the fixed features, and the features it lies near."""
    e = measure.edge_field(psi)
    k = int(np.argmax(e))
    xs, ys = float(measure.mid[k, 0]), float(measure.mid[k, 1])
    x_um, y_um = xs * um_per_s, ys * um_per_s
    cell = cell_size_um(mesh, xs, ys, um_per_s)
    dist = features.distances(x_um, y_um, window)
    near = [f for f in POINT_FEATURES if dist[f] <= NEAR_CELLS * cell]
    if "window_boundary" in dist and dist["window_boundary"] <= ON_BOUNDARY_TOL_UM:
        near.append("window_boundary")
    nearest = min(dist, key=dist.get)
    return {
        "hard_max_MV_cm": float(e[k]), "x_um": x_um, "y_um": y_um,
        "edge_dir": int(abs(mesh.node_xy[measure.b[k], 0] - mesh.node_xy[measure.a[k], 0]) == 0.0),
        "cell_um": cell, "distance_um": dist, "nearest_feature": nearest,
        "nearest_distance_um": dist[nearest], "nearest_distance_cells": dist[nearest] / cell,
        "near_features": near,
    }  # fmt: skip


def classify_pinned(records):
    """Pinned classification from the per-mesh :func:`hard_max_record` dicts (h, h/2, h/4).

    Pinned iff the features near the argmax on every mesh have a non-empty common set; returns
    ``{"pinned": bool, "pinning_features": [...]}``."""
    common = set(records[0]["near_features"])
    for r in records[1:]:
        common &= set(r["near_features"])
    return {"pinned": bool(common), "pinning_features": [f for f in FEATURES if f in common]}


# --------------------------------------------------------------------------- area-mean |E|
def mean_field(mesh, psi, c_mv_cm, um_per_s, centre_um, region, side_um=SIDE_UM):
    """Area-weighted mean of ``|E|`` [MV/cm] over the square of side ``side_um`` centred at
    ``centre_um`` (x, y in um), intersected with the SiC (``region='sic'``, y < 0) or the oxide
    (``region='oxide'``, y > 0).  See the module docstring for the construction.

    Returns {"mean_MV_cm", "area_um2", "area_fraction"} (fraction of ``side_um**2`` present)."""
    if region not in ("sic", "oxide"):
        raise ValueError("region must be 'sic' or 'oxide'")
    nx, ny = mesh.nx, mesh.ny
    x, y = mesh.x, mesh.y
    hx, hy = np.diff(x), np.diff(y)
    p2 = np.asarray(psi, dtype=float).reshape(nx + 1, ny + 1)  # node (i, j) = i (ny+1) + j
    ex = c_mv_cm * np.abs(np.diff(p2, axis=0)) / hx[:, None]  # (nx, ny+1): x-edge (i,j)-(i+1,j)
    ey = c_mv_cm * np.abs(np.diff(p2, axis=1)) / hy[None, :]  # (nx+1, ny): y-edge (i,j)-(i,j+1)
    h = 0.5 * side_um / um_per_s
    sx0, sx1 = centre_um[0] / um_per_s - h, centre_um[0] / um_per_s + h
    sy0, sy1 = centre_um[1] / um_per_s - h, centre_um[1] / um_per_s + h
    yc = 0.5 * (y[:-1] + y[1:])
    if np.any((y[:-1] < 0.0) & (y[1:] > 0.0)):
        raise ValueError("y = 0 must be a mesh line")
    ok_y = (yc < 0.0) if region == "sic" else (yc > 0.0)
    num = area = 0.0
    for qx in (0, 1):
        # quarter column: x in [x_i + qx hx/2, x_i + (qx+1) hx/2]
        xa = x[:-1] + qx * 0.5 * hx
        xb = xa + 0.5 * hx
        wx = np.clip(np.minimum(xb, sx1) - np.maximum(xa, sx0), 0.0, None)  # (nx,)
        for qy in (0, 1):
            ya = y[:-1] + qy * 0.5 * hy
            yb = ya + 0.5 * hy
            wy = np.clip(np.minimum(yb, sy1) - np.maximum(ya, sy0), 0.0, None) * ok_y  # (ny,)
            a = np.outer(wx, wy)  # scaled area of the intersection
            exq = ex[:, qy : qy + ny]  # x-edge on node row j + qy of column i
            eyq = ey[qx : qx + nx, :]  # y-edge on node column i + qx of row j
            emag = np.sqrt(exq**2 + eyq**2)
            num += float(np.sum(a * emag))
            area += float(np.sum(a))
    if area <= 0.0:
        raise ValueError("the clipped square has zero area")
    a_um2 = area * um_per_s**2
    return {"mean_MV_cm": num / area, "area_um2": a_um2, "area_fraction": a_um2 / side_um**2}


# ------------------------------------------------------------------------------- L1 / L2
def l1_assess(vals, order_min=ORDER_MIN, rich_tol=RICH_TOL):
    """L1 on three successive values (h, h/2, h/4); see the module docstring."""
    v0, v1, v2 = (float(v) for v in vals)
    e1, e2 = v1 - v0, v2 - v1
    out = {"values": [v0, v1, v2], "order": None, "richardson": None, "rel_dist": None,
           "monotone": False, "pass": False, "reason": None}  # fmt: skip
    if e1 == 0.0 and e2 == 0.0:
        return {**out, "monotone": True, "richardson": v2, "rel_dist": 0.0, "pass": True,
                "reason": "constant"}  # fmt: skip
    if e2 == 0.0 and e1 != 0.0:  # exact convergence on the finest step: order infinite
        return {**out, "monotone": True, "order": float("inf"), "richardson": v2,
                "rel_dist": 0.0, "pass": True, "reason": "e2 = 0"}  # fmt: skip
    if e1 * e2 <= 0.0:
        out["reason"] = "not monotone"
        return out
    p = float(np.log2(abs(e1 / e2)))
    rich = v2 + e2 / (2.0**p - 1.0)
    rel = abs(v2 - rich) / abs(rich) if rich != 0.0 else float("inf")
    ok = bool(p >= order_min and rel <= rich_tol)
    reason = None if ok else (f"order {p:.3g} < {order_min}" if p < order_min else
                              f"finest {100 * rel:.2g}% from Richardson")  # fmt: skip
    return {**out, "order": p, "richardson": float(rich), "rel_dist": float(rel),
            "monotone": True, "pass": ok, "reason": reason}  # fmt: skip


__all__ = [
    "FEATURES",
    "NEAR_CELLS",
    "ON_BOUNDARY_TOL_UM",
    "ORDER_MIN",
    "POINT_FEATURES",
    "RICH_TOL",
    "SIDE_UM",
    "Features",
    "Geometry",
    "classify_pinned",
    "hard_max_record",
    "l1_assess",
    "mean_field",
    "window_boundary_distance",
]
