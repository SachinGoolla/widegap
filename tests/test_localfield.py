"""Stage 15d: widegap.localfield (geometry features, area-mean |E|, pinned classification, L1)."""

from __future__ import annotations

import math

import numpy as np
import pytest

from widegap import localfield as lf
from widegap.adjoint_sic_oxide import FieldMeasure
from widegap.fv2d.mesh import TensorMesh

GEO = lf.Geometry(
    x_jfet=0.8, x_src0=1.6, x_src1=2.6, x_bc=3.0, x_max=3.4, gate_x1=1.6, src_contact_x0=1.64,
    d_src=0.3, d_body=0.8, t_ox=0.05, r_um=0.1, n_drift=5e16, n_jfet=1e17, n_body=3e17, n_src=1e19,
)  # fmt: skip


def _mesh():
    x = np.array([0.0, 0.3, 0.7, 1.0, 1.6, 2.0])
    y = np.array([-1.0, -0.55, -0.2, -0.05, 0.0, 0.03, 0.08])
    return TensorMesh(x, y)


# ------------------------------------------------------------------------- mean |E|
@pytest.mark.parametrize("region", ["sic", "oxide"])
@pytest.mark.parametrize("centre", [(1.1, -0.1), (0.0, 0.0), (1.7, 0.02), (0.05, -0.04)])
def test_mean_field_linear_potential_is_exact(region, centre):
    """psi = a x + b y has |E| = c sqrt(a^2 + b^2) everywhere: the area mean is that value for any
    clip.  (Pooling x- and y-edge components as separate samples would not give it.)"""
    m = _mesh()
    a, b, c, um = 0.7, -1.9, 1.5, 2.0
    psi = a * m.node_xy[:, 0] + b * m.node_xy[:, 1]
    cen_um = (centre[0] * um, centre[1] * um)
    r = lf.mean_field(m, psi, c, um, cen_um, region, side_um=0.6)
    assert r["mean_MV_cm"] == pytest.approx(c * math.hypot(a, b), rel=1e-13)


def _brute(m, psi, c, um, cen_um, region, side_um):
    """Independent quarter-cell integration through the mesh edge list."""
    lookup = {(int(a), int(b)): (float(h)) for a, b, h in zip(m.ea, m.eb, m.h, strict=True)}
    h2 = 0.5 * side_um / um
    sx0, sx1 = cen_um[0] / um - h2, cen_um[0] / um + h2
    sy0, sy1 = cen_um[1] / um - h2, cen_um[1] / um + h2
    num = area = 0.0
    for i in range(m.nx):
        for j in range(m.ny):
            yc = 0.5 * (m.y[j] + m.y[j + 1])
            if (region == "sic") != (yc < 0):
                continue
            hx, hy = m.x[i + 1] - m.x[i], m.y[j + 1] - m.y[j]
            for qx in (0, 1):
                for qy in (0, 1):
                    xa, ya = m.x[i] + qx * hx / 2, m.y[j] + qy * hy / 2
                    w = max(0.0, min(xa + hx / 2, sx1) - max(xa, sx0))
                    v = max(0.0, min(ya + hy / 2, sy1) - max(ya, sy0))
                    na, nb = m.node(i, j + qy), m.node(i + 1, j + qy)
                    ex = c * abs(psi[nb] - psi[na]) / lookup[(na, nb)]
                    na, nb = m.node(i + qx, j), m.node(i + qx, j + 1)
                    ey = c * abs(psi[nb] - psi[na]) / lookup[(na, nb)]
                    num += w * v * math.hypot(ex, ey)
                    area += w * v
    return num / area, area * um**2


@pytest.mark.parametrize("region", ["sic", "oxide"])
def test_mean_field_matches_edge_list_brute_force(region):
    m = _mesh()
    rng = np.random.default_rng(3)
    psi = rng.normal(size=m.n_nodes) + 0.4 * m.node_xy[:, 0] ** 2
    for cen in [(1.05, -0.12), (0.1, 0.01), (1.9, 0.04)]:
        um = 1.7
        got = lf.mean_field(m, psi, 2.0, um, (cen[0] * um, cen[1] * um), region, side_um=0.5)
        ref, area = _brute(m, psi, 2.0, um, (cen[0] * um, cen[1] * um), region, 0.5)
        assert got["mean_MV_cm"] == pytest.approx(ref, rel=1e-12)
        assert got["area_um2"] == pytest.approx(area, rel=1e-12)


def test_mean_field_clip_area_and_errors():
    m = _mesh()
    psi = np.zeros(m.n_nodes)
    psi[:] = m.node_xy[:, 1]
    # square centred on the interface at the symmetry plane, side 0.04 (scaled = um): quarter of
    # the square lies in the oxide (x >= 0, y > 0), a quarter in the SiC
    r = lf.mean_field(m, psi, 1.0, 1.0, (0.0, 0.0), "oxide", side_um=0.04)
    assert r["area_fraction"] == pytest.approx(0.25)
    r = lf.mean_field(m, psi, 1.0, 1.0, (0.0, 0.0), "sic", side_um=0.04)
    assert r["area_fraction"] == pytest.approx(0.25)
    with pytest.raises(ValueError):
        lf.mean_field(m, psi, 1.0, 1.0, (0.0, 0.0), "gold")
    with pytest.raises(ValueError):  # square entirely outside the region
        lf.mean_field(m, psi, 1.0, 1.0, (1.0, 0.5), "sic", side_um=0.04)
    bad = TensorMesh(np.array([0.0, 1.0]), np.array([-1.0, 0.5]))
    with pytest.raises(ValueError):  # y = 0 not a mesh line
        lf.mean_field(bad, np.zeros(bad.n_nodes), 1.0, 1.0, (0.5, 0.0), "sic")


# --------------------------------------------------------------------------- geometry
def test_contours_and_feature_distances():
    ft = lf.Features(GEO)
    assert len(ft.src_pts) > 50 and len(ft.body_pts) > 50
    # every contour point is a zero of the net doping (linear interpolation on a 5 nm grid)
    for pts in (ft.src_pts, ft.body_pts):
        s, o, na = GEO.doping(pts[:, 0], pts[:, 1])
        smooth = np.abs(pts[:, 1] + GEO.d_body) > 0.01  # the donor level jumps at y = -d_body
        assert np.max(np.abs((s + o - na) / na)[smooth]) < 0.05
    # n+ source junction lies around the source box; p-body junction near the body edges
    assert 1.4 < ft.src_pts[:, 0].min() < 1.6 and ft.src_pts[:, 1].min() > -0.7
    assert ft.body_pts[:, 1].min() < -0.8
    p = ft.body_pts[7]
    assert ft.distances(p[0], p[1])["p_body_junction"] < 1e-12
    d = ft.distances(GEO.gate_x1, GEO.t_ox)
    assert d["gate_edge"] == 0.0 and d["interface_corner"] == pytest.approx(GEO.t_ox)
    d = ft.distances(GEO.x_bc, 0.0, ((2.9, 3.1), (-0.1, 0.1)))
    assert d["interface_corner"] == 0.0 and d["window_boundary"] == pytest.approx(0.1)


def test_window_boundary_distance():
    w = ((0.0, 1.0), (-1.0, 0.0))
    assert lf.window_boundary_distance(0.0, -0.5, w) == 0.0  # on the boundary
    assert lf.window_boundary_distance(0.25, -0.5, w) == pytest.approx(0.25)  # inside
    assert lf.window_boundary_distance(1.3, -0.5, w) == pytest.approx(0.3)  # outside
    assert lf.window_boundary_distance(1.3, 0.4, w) == pytest.approx(0.5)  # outside corner


def test_cell_size_and_hard_max_record():
    m = _mesh()
    um = 2.0
    assert lf.cell_size_um(m, 0.5, -0.1, um) == pytest.approx(
        2.0 * 0.4
    )  # cell [0.3,0.7]x[-0.2,-0.05]
    assert lf.cell_size_um(m, 99.0, 99.0, um) == pytest.approx(
        2.0 * 0.4
    )  # clipped to the last cell
    # measure on three vertical edges; psi chosen so that the middle one is the hard maximum
    n = m.n_nodes
    psi = np.zeros(n)
    a, b = m.node(2, 1), m.node(2, 2)
    psi[b] = 3.0
    sel = np.array(
        [
            k
            for k in range(m.n_edges)
            if (m.ea[k], m.eb[k]) in {(a, b), (m.node(1, 1), m.node(1, 2))}
        ]
    )
    mid = 0.5 * (m.node_xy[m.ea[sel]] + m.node_xy[m.eb[sel]])
    fm = FieldMeasure(m.ea[sel], m.eb[sel], m.h[sel], 1.0, 20.0, mid)
    ft = lf.Features(GEO)
    rec = lf.hard_max_record(fm, psi, m, um, ft, window=((0.0, 2.0), (-1.0, 0.0)))
    assert rec["hard_max_MV_cm"] == pytest.approx(3.0 / (m.y[2] - m.y[1]))
    assert (rec["x_um"], rec["y_um"]) == pytest.approx((um * m.x[2], um * 0.5 * (m.y[1] + m.y[2])))
    assert rec["edge_dir"] == 1 and rec["cell_um"] > 0
    assert set(rec["distance_um"]) == set(lf.FEATURES)


def test_classify_pinned():
    r = lambda near: {"near_features": near}  # noqa: E731
    assert lf.classify_pinned([r(["p_body_junction"]), r(["p_body_junction", "gate_edge"]), r(["p_body_junction"])]) == {
        "pinned": True, "pinning_features": ["p_body_junction"]}  # fmt: skip
    # near a feature on only two meshes, or different features per mesh: not pinned
    assert not lf.classify_pinned([r(["gate_edge"]), r([]), r(["gate_edge"])])["pinned"]
    assert not lf.classify_pinned([r(["gate_edge"]), r(["interface_corner"]), r(["gate_edge"])])[
        "pinned"
    ]
    assert lf.classify_pinned([r(["window_boundary"])] * 3)["pinning_features"] == [
        "window_boundary"
    ]


# ------------------------------------------------------------------------------- L1
def test_l1_assess():
    # exact second-order sequence v(h) = 5 + h^2: order 2, finest within 2 % of Richardson
    v = [5 + 1.0, 5 + 0.25, 5 + 0.0625]
    r = lf.l1_assess(v)
    assert r["pass"] and r["order"] == pytest.approx(2.0) and r["richardson"] == pytest.approx(5.0)
    # first order passes the order test (>= 0.8) but a large remaining distance fails the 2 % test
    r = lf.l1_assess([1.0, 0.5, 0.25])  # p = 1, Richardson value 0: relative distance infinite
    assert r["order"] == pytest.approx(1.0) and not r["pass"] and r["rel_dist"] == float("inf")
    r = lf.l1_assess([10.0 + 4.0, 10.0 + 2.0, 10.0 + 1.0])
    assert r["pass"] is False and r["rel_dist"] == pytest.approx(1.0 / 10.0, rel=1e-9)
    assert "Richardson" in r["reason"]
    # slow convergence (order 0.5) fails on the order
    r = lf.l1_assess([0.0, 1.0, 1.0 + 1.0 / np.sqrt(2.0)])
    assert r["order"] == pytest.approx(0.5) and not r["pass"] and "order" in r["reason"]
    # non-monotone, and a diverging (negative-order) sequence
    assert lf.l1_assess([1.0, 2.0, 1.5])["reason"] == "not monotone"
    assert lf.l1_assess([1.0, 2.0, 4.0])["order"] < 0 and not lf.l1_assess([1.0, 2.0, 4.0])["pass"]
    # constant and exactly converged on the finest step
    assert lf.l1_assess([2.0, 2.0, 2.0])["pass"]
    assert lf.l1_assess([3.0, 2.0, 2.0])["pass"]
