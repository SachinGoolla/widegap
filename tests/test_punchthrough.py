"""Fast tests of the stage-19 1D punch-through column (widegap.punchthrough) on the small
4H-SiC DMOSFET mesh of test_offstate.py.  < 60 s."""

import numpy as np
import pytest

from widegap.offstate import CorrectedCornerAdjoint
from widegap.punchthrough import (
    Column1D,
    column_doping,
    first_crossing,
    interp_at,
    punch_voltage,
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
TH0 = np.zeros(5)


@pytest.fixture(scope="module")
def ev():
    sc = SiC4H(TracedRegistry.from_json())
    return CorrectedCornerAdjoint(sc, DESIGN, 0, r_um=0.1, cuts_um=(), tau_v=5.0)


def test_column_doping_matches_2d_node_doping_deep_in_the_source(ev):
    """At a column whose lateral neighbours carry the same doping (mid source, away from the
    erfc edges and the junction steps) the 1D control-volume average equals the 2D node value."""
    dev, mesh = ev.dev, ev.dev.mesh
    p = dev.prm
    x_mid = 0.5 * (p.x_src0 + p.x_src1)
    i = int(np.argmin(np.abs(mesh.x - x_mid)))
    t, nd, na, _, _ = column_doping(ev, float(mesh.x[i]), TH0)
    xy = mesh.node_xy
    nodes = np.flatnonzero((xy[:, 0] == mesh.x[i]) & (xy[:, 1] <= 0.0))
    nodes = nodes[np.argsort(-xy[nodes, 1])]
    assert np.allclose(t, -xy[nodes, 1])
    # compare the n+ source layer (interior nodes strictly above the junction) and the drain
    n2d, a2d = ev.nd[nodes], ev.na[nodes]
    src = (t > 0.0) & (t < 0.5 * p.d_src)
    assert np.allclose(nd[src], n2d[src], rtol=2e-2)
    assert np.allclose(nd[-1], n2d[-1], rtol=1e-12)
    # the acceptors are exactly the body profile; check the plateau value deep in the body
    k = int(np.argmin(np.abs(t - 0.5 * (p.d_src + p.d_body))))
    assert 0.95 < na[k] / ev.dev.prm.a_body <= 1.0
    assert a2d[k] > 0


def test_tied_bc_adds_exactly_one_hole_row(ev):
    fl = Column1D(ev, 2.0, TH0, "floating")
    td = Column1D(ev, 2.0, TH0, "tied")
    rf, _ = fl.model._all_bc(fl.bc(0.5))
    rt, vt = td.model._all_bc(td.bc(0.5))
    assert len(rt) == len(rf) + 1
    extra = sorted(set(rt.tolist()) - set(rf.tolist()))
    assert extra == [3 * td.i_body + 2]
    assert vt[list(rt).index(extra[0])] == 0.0
    assert 0 < td.body_node_um < 1.0


def test_floating_march_converges_and_current_is_constant_over_the_barrier(ev):
    col = Column1D(ev, 2.0, TH0, "floating")
    eq = col.equilibrium()
    assert eq.report.converged
    o0 = col.observe(eq)
    assert o0["Jn_A_cm2"] == 0.0 and o0["barrier_VT"] > 100.0
    res = col.march(eq.X, 0.0, [1.0, 2.0, 3.0])
    assert all(r[1].report.converged for r in res)
    s = res[-1][1]
    o = col.observe(s)
    assert 0.0 < abs(o["Jn_A_cm2"]) < 1e-20
    # flux is constant over the edges that resolve it (barrier and depletion edges): the spread
    # over the edges between the barrier minimum and the end of the window is round-off
    js = col.j_scale()
    jn = np.asarray(s.jn) * js
    k = int(np.argmin(s.psi[col.win]))
    sel = jn[k : int(np.flatnonzero(col.win)[-1])]
    assert np.max(np.abs(sel - sel.mean())) <= 1e-9 * abs(sel.mean())
    assert abs(o["phi_p_body_V"]) < 0.1
    # the barrier shrinks with bias
    assert o["barrier_VT"] < o0["barrier_VT"]


def test_crossing_and_punch_voltage_helpers():
    v = np.arange(0.0, 101.0, 1.0)
    j = np.exp(0.1 * v)  # exp growth: |J(10)| = e, threshold 1e3 e at V = 10 + ln(1e3)/0.1
    out = punch_voltage(v, j, 10.0, 1e3)
    assert out["V_pt"] == pytest.approx(10.0 + np.log(1e3) / 0.1, rel=1e-12)
    assert first_crossing(v, 100.0 - v, 50.0, above=False) == pytest.approx(50.0)
    assert first_crossing(v, v, 1e9) is None
    assert interp_at(v, 2.0 * v, 3.5) == pytest.approx(7.0)
    assert interp_at(v, v, 200.0) is None
