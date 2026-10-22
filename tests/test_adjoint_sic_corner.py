"""Fast tests of the stage-13b corner-aware objective (widegap.adjoint_sic_corner) on a small
4H-SiC DMOSFET mesh (187 nodes): rounded doping, soft-min, and adjoint gradients of the
multi-cut soft-min BV and of R_on,sp against Richardson-extrapolated central differences
(eps = 1e-3, 5e-4).  < 60 s."""

import numpy as np
import pytest

from widegap.adjoint_sic import richardson_fd
from widegap.adjoint_sic_corner import (
    CornerAdjoint,
    DenseCuts,
    RoundedDMOSFET,
    dense_bv,
    smooth_step,
    softmin,
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
R_UM = 0.1
CUTS_UM = (0.8, 1.6, 3.0)  # BV 469 / 450 / 482 V at TH on this mesh (all within 7%)
TAU = 5.0


@pytest.fixture(scope="module")
def sc():
    return SiC4H(TracedRegistry.from_json())


@pytest.fixture(scope="module")
def adj(sc):
    return CornerAdjoint(sc, DESIGN, 0, r_um=R_UM, cuts_um=CUTS_UM, tau_v=TAU)


def _unit(seed=0):
    u = np.random.default_rng(seed).standard_normal(5)
    return u / np.linalg.norm(u)


def test_smooth_step_and_rounded_doping(adj):
    r = 0.7
    s = np.linspace(-6, 6, 13)
    f = smooth_step(s, r)
    assert f[0] < 1e-12 and f[-1] > 1 - 1e-12 and f[6] == pytest.approx(0.5)
    assert np.allclose(f + smooth_step(-s, r), 1.0)
    # variance of the smoothed edge = r^2
    ds = np.linspace(-8, 8, 4001)
    pdf = np.gradient(smooth_step(ds, r), ds)
    assert np.sqrt(np.trapezoid(pdf * ds**2, ds) / np.trapezoid(pdf, ds)) == pytest.approx(r, rel=1e-3)
    dev = adj.dev
    assert isinstance(dev, RoundedDMOSFET)
    p = dev.prm
    a, d = dev.cell_acceptors(), dev.cell_donors(np.zeros(5))
    assert a.max() == pytest.approx(p.a_body, rel=1e-6)  # deep body reaches the plateau
    assert d.max() == pytest.approx(p.b_src, rel=1e-3) or d.max() == pytest.approx(p.b_drain)
    # rounded edge: acceptors vary smoothly (no cell jump larger than the plateau fraction)
    cc = dev.mesh.cell_centres
    row = (np.abs(cc[:, 1] + 0.5 * p.d_body) == np.abs(cc[:, 1] + 0.5 * p.d_body).min()) & dev.reg.semi_cell
    assert 0.0 < a[row][(cc[row, 0] > p.x_jfet - 3 * dev.r) & (cc[row, 0] < p.x_jfet + 3 * dev.r)].min()
    # compensated node fields: donors and acceptors both present near the junction
    assert np.any((adj.nd > 0) & (adj.na > 0))
    assert np.all(adj.dev.cell_doping(np.zeros(5)) == d - a)


def test_cell_doping_jac_matches_central_difference(adj):
    dev = adj.dev
    jac = dev.cell_doping_jac(TH)
    for k in range(5):
        e = np.zeros(5)
        e[k] = 1e-6
        fd = (dev.cell_donors(TH + e) - dev.cell_donors(TH - e)) / 2e-6
        assert np.allclose(jac[:, k], fd, rtol=1e-7, atol=1e-9)
    # acceptors do not depend on theta: the net-doping Jacobian equals the donor Jacobian
    assert np.allclose(dev.M @ jac, adj.dev.M @ jac)


def test_softmin_bounds_weights_and_gradient():
    bv = np.array([412.0, 415.0, 503.0, 518.0, 400.5])
    tau = 5.0
    obj, w = softmin(bv, tau)
    assert bv.min() - tau * np.log(len(bv)) <= obj <= bv.min()
    assert w.sum() == pytest.approx(1.0) and np.all(w > 0)
    for k in range(len(bv)):
        e = np.zeros(len(bv))
        e[k] = 1e-5
        fd = (softmin(bv + e, tau)[0] - softmin(bv - e, tau)[0]) / 2e-5
        assert fd == pytest.approx(w[k], rel=1e-6, abs=2e-8)  # abs: round-off floor of the FD
    # tau -> 0 recovers the hard minimum
    assert softmin(bv, 1e-3)[0] == pytest.approx(bv.min(), abs=1e-9)


def test_cut_positions_are_mesh_columns_and_cold_vs_warm(adj):
    assert len(adj.cuts) == len(CUTS_UM)
    assert np.allclose(adj.cut_x_um, [adj.sc.length_cm(x) * 1e4 for x in adj.cut_x])
    adj.reset_caches()
    cold = adj.objective(TH, grad=False, warm=False)
    warm = adj.objective(TH, grad=False)  # warm start from the cold caches
    assert warm["BV"] == pytest.approx(cold["BV"], abs=1e-6)
    assert np.all(cold["bv_k"] > 100.0)
    assert cold["BV"] <= cold["bv_k"].min() + 1e-12


def test_objective_gradient_vs_richardson_fd(adj):
    adj.tau_v = 30.0  # mixed soft-min weights (all three cuts contribute)
    try:
        obj = adj.objective(TH)
        assert obj["w"].min() > 0.05
        g = obj["grad"]
        assert np.allclose(g, sum(w * gk for w, gk in zip(obj["w"], obj["grad_k"], strict=True)))
        for seed in (0, 1):
            u = _unit(seed)
            rich, f1, f2 = richardson_fd(lambda t: adj.objective(t, grad=False)["BV"], TH, u, 1e-3)
            assert float(g @ u) == pytest.approx(rich, rel=1e-4, abs=1e-6)
    finally:
        adj.tau_v = TAU


def test_ron_gradient_on_rounded_device(adj):
    r = adj.ron(TH)
    u = _unit(2)
    rich, _, _ = richardson_fd(lambda t: adj.ron(t, grad=False)["R"], TH, u, 1e-3)
    assert float(r["grad"] @ u) == pytest.approx(rich, rel=1e-5)


def test_dense_cut_bv_matches_objective_cut_on_a_mesh_column(sc, adj):
    """The independent dense machinery (bilinear lines, global march + bracket refinement)
    reproduces the adjoint-model BV of a cut that sits on a mesh column."""
    ob = adj.objective(TH, grad=False)
    xs = np.array(adj.cut_x_um)
    dense = DenseCuts(adj, xs)
    res = dense_bv(adj, TH, dense, v_step_v=20.0, tol_v=1e-3, cap=1.15)
    assert np.allclose(res["bv"], ob["bv_k"], atol=0.02)
    assert res["bv_min"] == pytest.approx(ob["bv_k"].min(), abs=0.02)


def test_cached_contact_bc_equals_model_ohmic(adj):
    """The cached contact potentials equal the parent's ``model.ohmic`` (root find per node)."""
    from widegap.fv2d.mos import Dirichlet

    b = adj.base
    rows, vals = b.model._all_bc(b.bc(1.3, 0.0, 2.7))
    bc = Dirichlet()
    b.model.ohmic(bc, b.dev.nodes["source_body"], 0.0)
    b.model.ohmic(bc, b.dev.nodes["drain"], 2.7)
    bc.gate(b.dev.nodes["gate"], 1.3, b.vfb_code)
    rows0, vals0 = b.model._all_bc(bc)
    assert np.array_equal(rows, rows0) and np.allclose(vals, vals0, rtol=0, atol=1e-12)
