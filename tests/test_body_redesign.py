"""Fast tests of the stage-20 body-depth family (widegap.bodyredesign).  < 60 s."""

import numpy as np
import pytest

from widegap.adjoint_sic_oxide import OxideAdjoint
from widegap.bodyredesign import (
    OffMarch,
    barrier_min,
    body_design,
    column_index_matrix,
    source_columns,
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


@pytest.fixture(scope="module")
def sc():
    return SiC4H(TracedRegistry.from_json())


def make(sc, d):
    return OxideAdjoint(sc, body_design(DESIGN, d), 0, r_um=0.1, cuts_um=(), tau_v=5.0)


def test_body_design_changes_only_the_body_depth():
    base = DESIGN
    d2 = body_design(base, 1.2)
    assert d2.d_body == 1.2
    for f in ("x_jfet", "ch_len", "src_w", "d_src", "w_drift", "t_drain", "t_ox", "n_drift",
              "n_jfet", "n_body", "n_src", "n_drain"):  # fmt: skip
        assert getattr(d2, f) == getattr(base, f)
    assert {k: v for k, v in d2.mesh_kw.items() if k != "ny_body"} == {
        k: v for k, v in base.mesh_kw.items() if k != "ny_body"
    }
    # the stage-12 design (0.5 um of body below the source) keeps its cell count
    from dataclasses import replace

    s12 = replace(base, mesh_kw={**base.mesh_kw, "ny_body": 11})
    assert body_design(s12, 0.8).mesh_kw["ny_body"] == 11
    assert body_design(s12, 1.6).mesh_kw["ny_body"] == 29


def test_surface_and_channel_doping_unchanged_and_drift_follows_the_body(sc):
    a, b = make(sc, 0.8), make(sc, 1.2)
    ma, mb = a.dev.mesh, b.dev.mesh
    p = a.dev.prm
    ya = {round(float(y), 9) for y in ma.y if y >= -p.d_src - 1e-12}
    yb = {round(float(y), 9) for y in mb.y if y >= -p.d_src - 1e-12}
    assert ya == yb
    xa = np.round(ma.node_xy[:, 0], 9)
    ya_n = np.round(ma.node_xy[:, 1], 9)
    xb = np.round(mb.node_xy[:, 0], 9)
    yb_n = np.round(mb.node_xy[:, 1], 9)
    top_a = np.flatnonzero(ya_n > round(-p.d_src, 9) + 1e-6)
    top_b = np.flatnonzero(yb_n > round(-p.d_src, 9) + 1e-6)
    ka = sorted(zip(xa[top_a], ya_n[top_a], top_a, strict=True))
    kb = sorted(zip(xb[top_b], yb_n[top_b], top_b, strict=True))
    assert [k[:2] for k in ka] == [k[:2] for k in kb]
    ia, ib = [k[2] for k in ka], [k[2] for k in kb]
    # the acceptor erfc tail at 0.3-0.5 um above the junction is the only difference allowed
    assert np.allclose(a.nd[ia], b.nd[ib], rtol=1e-6, atol=1e-6)
    assert np.allclose(a.na[ia], b.na[ib], rtol=1e-6, atol=1e-6)
    # total depth grows by exactly the body-depth change (drift thickness unchanged)
    assert (ma.y.min() - mb.y.min()) == pytest.approx(float(sc.length_s(0.4e-4)), rel=1e-9)


def test_barrier_min_on_a_toy_profile():
    # two columns x 5 nodes: surface first
    psi = np.array([10.0, 8.0, 3.0, 5.0, 9.0, 7.0, 6.0, 5.0, 6.0, 8.0])
    idx = np.array([[0, 1, 2, 3, 4], [5, 6, 7, 8, 9]])
    m, k, bars = barrier_min(psi, idx)
    assert bars.tolist() == [7.0, 2.0]
    assert (m, k) == (2.0, 1)


def test_source_columns_and_offmarch_barrier_decreases_with_bias(sc):
    ev = make(sc, 1.0)
    ev.set_theta(np.zeros(5))
    cols = source_columns(ev)
    p = ev.dev.prm
    assert cols and all(p.x_src0 - 1e-9 <= ev.dev.mesh.x[i] <= p.x_bc + 1e-9 for i in cols)
    idx = column_index_matrix(ev, cols)
    assert idx.shape[0] == len(cols)
    m = OffMarch(ev, np.zeros(5))
    st = (m.v, m.x)
    b0 = barrier_min(st[1][0::3], idx)[0]
    st2 = m.advance(st, 30.0)
    b1 = barrier_min(st2[1][0::3], idx)[0]
    assert b0 > 100.0 and b1 < b0
    r = OffMarch(ev, np.zeros(5)).vpt_es(idx, level_vt=b1, step_v=15.0, tol_v=0.5)
    assert r["V_pt_es"] == pytest.approx(30.0, abs=1.0)
