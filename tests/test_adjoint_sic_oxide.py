"""Fast tests of the stage-13d fixed-voltage field measures (widegap.adjoint_sic_oxide) on a small
4H-SiC DMOSFET mesh (187 nodes): smooth maximum, oxide/window edge sets, and the adjoint gradients
of E_ox and W1-W3 at V_blk against Richardson-extrapolated central differences (eps = 1e-3, 5e-4).
< 60 s."""

import numpy as np
import pytest

from widegap.adjoint_sic import richardson_fd
from widegap.adjoint_sic_oxide import (
    KAPPA_PER_MV_CM,
    MEASURES,
    WINDOWS_UM,
    FieldMeasure,
    OxideAdjoint,
    smooth_max,
)
from widegap.sic import DMOSDesign, SiC4H, TracedRegistry

pytestmark = pytest.mark.filterwarnings(
    "ignore::RuntimeWarning", "ignore::scipy.sparse.linalg.MatrixRankWarning"
)

DESIGN = DMOSDesign(
    x_jfet=0.8, ch_len=0.8, src_w=1.0, bc_w=0.8, pad_w=0.4,
    d_src=0.3, d_body=0.8, w_drift=3.0, t_drain=0.3, t_ox=0.05,
    n_drift=5e16, n_jfet=1e17, n_body=3e17, n_src=1e19, n_drain=1e19,
    mesh_kw=dict(nx=(2, 1, 2, 2, 2, 1), ny_drain=1, ny_drift=5, ny_body=4, ny_surf=5, ny_ox=1,
                 q_surf=0.6),
)  # fmt: skip
TH = np.array([0.1, -0.05, 0.03, 0.02, 0.1])
V_TEST = 100.0  # fixed blocking voltage of the fast tests (spec value 400 V is used by the run)


@pytest.fixture(scope="module")
def sc():
    return SiC4H(TracedRegistry.from_json())


@pytest.fixture(scope="module")
def ev(sc):
    return OxideAdjoint(sc, DESIGN, 0, v_blk_v=V_TEST)


def test_smooth_max_bounds_and_weights():
    e = np.array([1.0, 2.5, 2.4, 0.1])
    s, w = smooth_max(e, KAPPA_PER_MV_CM)
    assert e.max() <= s <= e.max() + np.log(len(e)) / KAPPA_PER_MV_CM
    assert w.sum() == pytest.approx(1.0) and np.argmax(w) == 1
    # translation covariance, and no overflow for large arguments
    assert smooth_max(e + 1e3, KAPPA_PER_MV_CM)[0] == pytest.approx(s + 1e3)


def test_field_measure_psi_gradient():
    rng = np.random.default_rng(1)
    a, b = np.arange(6), np.arange(1, 7)
    m = FieldMeasure(a, b, np.full(6, 0.3), 0.7, 5.0, np.zeros((6, 2)))
    psi = rng.standard_normal(7)
    s, g = m.value_grad(psi)
    assert s == pytest.approx(m.value(psi))
    for i in range(7):
        d = np.zeros(7)
        d[i] = 1e-6
        fd = (m.value(psi + d) - m.value(psi - d)) / 2e-6
        assert g[i] == pytest.approx(fd, rel=1e-6, abs=1e-9)
    with pytest.raises(ValueError):
        FieldMeasure([], [], [], 1.0, 1.0, np.zeros((0, 2)))


def test_edge_sets(ev, sc):
    mesh, p = ev.mesh, ev.prm
    mo = ev.measures["E_ox"]
    ya, yb = mesh.node_xy[mo.a, 1], mesh.node_xy[mo.b, 1]
    assert np.all(np.minimum(ya, yb) == 0.0) and np.all(np.maximum(ya, yb) > 0.0)
    assert len(mo.a) == len(mesh.x)  # one interface-adjacent oxide edge per mesh column
    # windows: non-empty, semiconductor only, pairwise disjoint boxes, recorded values
    assert set(WINDOWS_UM) == {"W1", "W2", "W3"}
    masks = list(ev.window_masks.values())
    for m in masks:
        assert m.sum() > 0 and np.all(ev.dev.reg.w_s[m] > 0)
    assert not (masks[0] & masks[1]).any() and not (masks[0] & masks[2]).any()
    assert not (masks[1] & masks[2]).any()
    assert WINDOWS_UM["W1"] == ((0.5, 1.1), (-1.1, -0.5))
    assert WINDOWS_UM["W2"] == ((1.3, 2.0), (-1.2, -0.6)) and KAPPA_PER_MV_CM == 20.0
    # the corner (x_jfet, d_body) lies in W1, the source edge at the body junction in W2
    ux = float(sc.l_d) * 1e4
    assert 0.5 < p.x_jfet * ux < 1.1 and 1.3 < p.x_src0 * ux < 2.0


def test_fields_consistent_with_hard_max(ev):
    o = ev.fields(TH, grad=False)
    for k in MEASURES:
        s, h = o[k]["S"], o[k]["hard_max"]
        n = len(ev.measures[k].a)
        assert h <= s <= h + np.log(n) / KAPPA_PER_MV_CM + 1e-12
    assert 0.0 < o["E_ox"]["S"] < 50.0


def test_state_warm_equals_cold(sc):
    ev1 = OxideAdjoint(sc, DESIGN, 0, v_blk_v=V_TEST)
    x_cold = ev1.state_blk(TH).copy()
    th2 = TH + 0.01
    x_warm = ev1.state_blk(th2).copy()  # warm start from the previous state
    ev2 = OxideAdjoint(sc, DESIGN, 0, v_blk_v=V_TEST)
    x_ref = ev2.state_blk(th2)  # march
    assert np.max(np.abs(x_warm[0::3] - x_ref[0::3])) < 1e-8 * max(1.0, np.max(np.abs(x_ref[0::3])))
    assert not np.allclose(x_cold[0::3], x_ref[0::3])


@pytest.mark.parametrize("seed", [0, 1])
def test_adjoint_matches_richardson_fd(ev, seed):
    u = np.random.default_rng(seed).standard_normal(5)
    u /= np.linalg.norm(u)
    o = ev.fields(TH)
    # nodal sensitivity consistent with the theta gradient (before FD calls change the theta)
    assert np.allclose(o["E_ox"]["grad"], ev.dnd_dtheta.T @ o["E_ox"]["g_nd"])
    for k in MEASURES:
        rich, f1, f2 = richardson_fd(lambda t, k=k: ev.fields(t, grad=False, which=(k,))[k]["S"],
                                     TH, u, 1e-3)  # fmt: skip
        a = float(o[k]["grad"] @ u)
        assert a == pytest.approx(rich, rel=1e-8, abs=1e-9), k
        assert abs(f1 - f2) < 1e-4 * max(abs(rich), 1e-3)
