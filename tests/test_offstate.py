"""Fast tests of the stage-12c corrected off-state model (widegap.offstate) on a small 4H-SiC
DMOSFET mesh (187 nodes) and a small laterally uniform p+n diode.  < 60 s."""

import numpy as np
import pytest

from widegap.adjoint_sic import richardson_fd
from widegap.adjoint_sic_corner import solve_at_voltage
from widegap.fv.avalanche import avalanche_bv
from widegap.offstate import (
    CorrectedCornerAdjoint,
    UniformPN2D,
    corrected_class,
    jacobian_fd_rows,
    n_consistency,
)
from widegap.sic import DMOSDesign, IonizedReverseDiode, SiC4H, SiCDMOS, TracedRegistry

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
CUTS_UM = (0.8, 1.6, 3.0)


@pytest.fixture(scope="module")
def sc():
    return SiC4H(TracedRegistry.from_json())


@pytest.fixture(scope="module")
def adj(sc):
    return CorrectedCornerAdjoint(sc, DESIGN, 0, r_um=0.1, cuts_um=CUTS_UM, tau_v=5.0)


@pytest.fixture(scope="module")
def adj_dist(sc):
    return CorrectedCornerAdjoint(sc, DESIGN, 0, r_um=0.1, cuts_um=CUTS_UM, tau_v=5.0,
                                  w_mode="distance")  # fmt: skip


def _unit(seed=0):
    u = np.random.default_rng(seed).standard_normal(5)
    return u / np.linalg.norm(u)


def test_w_field_dirichlet_values_laplace_and_maximum_principle(adj, adj_dist):
    for ev in (adj, adj_dist):
        base = ev.base
        w = base.offstate_w
        zero, one = base._w_sets
        semi = base.dev.reg.semi_node
        assert len(zero) > 0 and len(one) > 0
        assert np.all(w[zero] == 0.0) and np.all(w[one] == 1.0)
        assert np.all(np.isfinite(w[semi])) and np.all(np.isnan(w[~semi]))
        assert w[semi].min() >= 0.0 and w[semi].max() <= 1.0
        assert np.all(base.nd[zero] > base.na[zero]) and np.all(base.nd[one] > base.na[one])
        # the n+ source is pinned to the source; the drain and the drift/JFET to the drain
        drain = base.dev.nodes["drain"]
        assert np.all(w[drain] == 1.0)
        mid = semi & (w > 0.01) & (w < 0.99)
        assert mid.sum() > 10
    # harmonic: sum_j gw_ij (w_i - w_j) = 0 on every non-Dirichlet semiconductor node
    base, mod = adj.base, adj.base.model
    w = np.nan_to_num(base.offstate_w)
    dw = w[mod.b] - w[mod.a]
    res = np.bincount(mod.a, mod.gw * dw, minlength=len(w)) - np.bincount(
        mod.b, mod.gw * dw, minlength=len(w)
    )
    zero, one = base._w_sets
    pinned = np.zeros(len(w), dtype=bool)
    pinned[np.r_[zero, one]] = True
    sel = base.dev.reg.semi_node & ~pinned
    scale = np.bincount(mod.a, mod.gw, minlength=len(w)) + np.bincount(
        mod.b, mod.gw, minlength=len(w)
    )
    assert np.max(np.abs(res[sel]) / scale[sel]) < 1e-12
    # the two interpolations agree on the pinned sets and differ in between
    wd = adj_dist.base.offstate_w
    assert np.array_equal(np.isnan(wd), np.isnan(adj.base.offstate_w))
    assert np.nanmax(np.abs(wd - adj.base.offstate_w)) > 1e-3


def test_offstate_bc_structure_and_dvals(adj):
    base = adj.base
    mod = adj.mod
    w = base.offstate_w
    rows0, v0 = mod._all_bc(base.offstate_bc(0.0))
    rows1, v1 = mod._all_bc(base.offstate_bc(7.0))
    assert np.array_equal(rows0, rows1)
    pn_rows = rows1[rows1 % 3 == 1]
    pp_rows = rows1[rows1 % 3 == 2]
    sel = np.isin(rows1, pn_rows)
    free = base._w_free
    d = dict(zip(rows1.tolist(), (v1 - v0).tolist(), strict=True))
    # the V_D-derivative of the Dirichlet data: w on the free phi_n rows, 0 on phi_p rows, 1 on
    # the drain contact, 0 on the source/body contact
    assert np.allclose([d[3 * i + 1] for i in free], 7.0 * w[free])
    assert all(d[int(r)] == 0.0 for r in pp_rows if (r // 3) in set(free.tolist()))
    assert np.allclose([d[3 * i + 1] for i in base.dev.nodes["drain"]], 7.0)
    assert np.allclose([d[3 * i + 1] for i in base.dev.nodes["source_body"]], 0.0)
    assert sel.sum() > 0
    # zero bias: phi_n = phi_p = 0 on every node (equilibrium)
    assert np.allclose(v0[rows0 % 3 != 0], 0.0)
    # the adjoint tangent machinery sees exactly w
    rows, dv = adj._dvals(3.0)
    assert np.allclose(dv[np.isin(rows, 3 * free + 1)], w[free])


def test_n_negligible_where_w_is_intermediate_and_state_is_smooth(sc, adj, adj_dist):
    for ev in (adj, adj_dist):
        for v in (100.0, 300.0):
            r = n_consistency(ev.base, solve_at_voltage(ev, np.zeros(5), v))
            assert r["ok"] and r["n_nodes_checked"] > 10
            assert r["free_nodes"]["worst_ratio"] <= r["tol"]


def test_jacobian_matches_central_difference_per_row(adj):
    ev = adj
    X = solve_at_voltage(ev, np.zeros(5), 200.0)
    rows, vals = ev._off_bc(200.0 / ev.sc.vt)
    out = jacobian_fd_rows(ev.mod, X, rows, vals, step=1e-2)
    assert out["err"].max() <= 1e-6
    # the round-off floor is far below the gate, so the gate is meaningful
    assert out["floor"].max() < 1e-7


def test_bv_adjoint_gradient_and_voltage_derivative(adj):
    adj.tau_v = 80.0
    try:
        obj = adj.objective(TH)
        assert obj["w"].min() > 0.02
        for seed in (0, 1):
            u = _unit(seed)
            rich, _, _ = richardson_fd(lambda t: adj.objective(t, grad=False)["BV"], TH, u, 1e-3)
            assert float(obj["grad"] @ u) == pytest.approx(rich, rel=1e-4, abs=1e-6)
        out = adj.bv(TH, grad=True)  # dI/dV: adjoint (lam^T dvals) vs tangent
        assert out["di_dv_adjoint"] == pytest.approx(out["di_dv_tangent"], rel=1e-6)
    finally:
        adj.tau_v = 5.0


def test_uniform_pn2d_reproduces_the_1d_diode(sc):
    """Laterally uniform p+n (N_D = 1e17, level 0): with w = 1 on the p+ side the 2D corrected
    model is the 12a 1D model (phi_n = V_D everywhere, phi_p = 0)."""
    ch = sc.niwa()
    pn = UniformPN2D(sc, 1e17, level=0)
    assert np.all(pn.w("harmonic")[pn.reg.semi_node] == pytest.approx(1.0, abs=1e-12))
    assert np.allclose(pn.w("distance")[pn.reg.semi_node], 1.0)
    ref = avalanche_bv(
        pn.diode, ch, "n", v_step=pn.v_est / 3.0, v_max=6.0 * pn.v_est, xtol=1e-12
    ).bv
    assert abs(pn.bv(ch, "new", "harmonic") - ref) / ref <= 1e-6
    assert abs(pn.bv(ch, "new", "distance") - ref) / ref <= 1e-6


def test_sharp_device_planar_cut_bv_old_model_is_biased_low_vs_1d_column(sc):
    """The planar cut is NOT unaffected by the model change: the old rule (phi_n = 0 on the
    p-type junction node) puts an electron reservoir (n ~ N_A) at the body/drift junction node
    and suppresses the p-side depletion, so its BV is low.  The corrected 2D BV lies between the
    old one and the flat-quasi-Fermi 1D column diode (no lateral fringing)."""
    ch = sc.niwa()
    out = {}
    for name, cls in (("old", SiCDMOS), ("new", corrected_class(SiCDMOS))):
        d = cls(sc, DESIGN, 0)
        eq = d.equilibrium()
        r = d.offstate_bv(ch, sc.volt_s(20.0), sc.volt_s(1500.0), eq.X, "n", x_min=d.prm.x_bc)
        out[name] = sc.volt_phys(r["bv"])
        if name == "old":  # spurious junction-node electrons in the old state
            nodes = d.column_nodes(d.prm.x_bc)[1]
            psi, pn = r["X"][0::3][nodes], r["X"][1::3][nodes]
            n_ratio = np.exp(d.model.ld + psi - pn) / np.maximum(d.nd[nodes], d.na[nodes])
            assert n_ratio[:-2].max() > 1.0  # n >> N at the junction node (n ~ 10 N_A)
    col, nodes = d.column_nodes(d.prm.x_bc)
    dio = IonizedReverseDiode(-d.dev.mesh.node_xy[nodes, 1], d.nd[nodes], d.na[nodes], sc.delta,
                              d.donors, d.acceptors)  # fmt: skip
    b1d = sc.volt_phys(avalanche_bv(dio, ch, "n", v_step=float(sc.volt_s(20.0)),
                                    v_max=float(sc.volt_s(1500.0)), xtol=1e-12).bv)  # fmt: skip
    assert out["old"] < out["new"] < b1d
    assert abs(out["new"] - b1d) < 0.5 * abs(out["old"] - b1d)
