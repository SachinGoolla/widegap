"""Fast tests of the stage-17 exact adjoint through the channel mobility (widegap.adjoint_channel)
on a small 4H-SiC DMOSFET mesh (264 nodes, 7 nodes inside d_ch = 20 nm of each channel column):
d mu_edge / dX against central differences of the forward mobility map, F_mu against a residual
finite difference, and the R_on,ch gradient (exact and lagged) against Richardson-extrapolated
central differences of the Picard fixed point.  < 60 s."""

import numpy as np
import pytest

from widegap.adjoint_channel import ChannelRon, fmu_matrix, noguchi_dmu_dns
from widegap.adjoint_sic import richardson_fd
from widegap.adjoint_sic_oxide import OxideAdjoint
from widegap.channel_mobility import NS_MIN_CM2, NoguchiParams, noguchi_mobility
from widegap.sic import DMOSDesign, SiC4H, TracedRegistry

pytestmark = pytest.mark.filterwarnings(
    "ignore::RuntimeWarning", "ignore::scipy.sparse.linalg.MatrixRankWarning"
)

DESIGN = DMOSDesign(
    x_jfet=0.8, ch_len=0.8, src_w=1.0, bc_w=0.8, pad_w=0.4,
    d_src=0.3, d_body=0.8, w_drift=3.0, t_drain=0.3, t_ox=0.05,
    n_drift=5e16, n_jfet=1e17, n_body=3e17, n_src=1e19, n_drain=1e19,
    mesh_kw=dict(nx=(2, 1, 2, 2, 2, 1), ny_drain=1, ny_drift=5, ny_body=4, ny_surf=12, ny_ox=1,
                 q_surf=0.6),
)  # fmt: skip
TH = np.array([0.1, -0.05, 0.03, 0.02, 0.1])


@pytest.fixture(scope="module")
def stack():
    reg = TracedRegistry.from_json()
    return SiC4H(reg), NoguchiParams.from_registry(reg)


@pytest.fixture(scope="module")
def cr(stack):
    sc, p = stack
    return ChannelRon(OxideAdjoint(sc, DESIGN, 0), p)


@pytest.fixture(scope="module")
def point(cr):
    return cr.ron(TH, detail=True)


def test_noguchi_derivative_and_clamps(stack):
    sc, p = stack
    na = np.array([1e17, 3e17, 1e17, 2e17, 1e10])
    ns = np.array([3e12, 8e12, 2e11, 5e10, 1e11])  # last: E_eff clamp; 4th: N_S clamp
    mu, dmu, r = noguchi_dmu_dns(sc, p, na, ns)
    assert np.array_equal(mu, noguchi_mobility(sc, p, na, ns)["mu"])
    assert r["clamp_ns"].tolist() == [False, False, False, True, False]
    assert r["clamp_eeff"][4] and not r["clamp_eeff"][:4].any()
    h = 1e-6 * ns
    fd = (noguchi_mobility(sc, p, na, ns + h)["mu"] - noguchi_mobility(sc, p, na, ns - h)["mu"]) / (
        2 * h
    )
    assert np.allclose(dmu[:3], fd[:3], rtol=1e-7)
    assert dmu[3] == 0.0  # clamped N_S: zero derivative
    assert ns[3] < NS_MIN_CM2
    # E_eff clamp active: only the Coulomb term contributes
    hh = 1e-6 * ns[4]
    f4 = (
        noguchi_mobility(sc, p, na[4:], ns[4:] + hh)["mu"]
        - noguchi_mobility(sc, p, na[4:], ns[4:] - hh)["mu"]
    ) / (2 * hh)
    assert dmu[4] == pytest.approx(f4[0], rel=1e-7)


def test_test_device_is_not_trivial(cr, point):
    """The checks below can only pass trivially on a device that exercises the channel model."""
    cm = cr.cmap
    assert cm.col_gated.sum() >= 2 and len(cm.edges) >= 10
    assert min(len(v) for v in cm.col_nodes.values()) >= 3  # trapezoid over >= 3 nodes
    g = point["guard"]
    assert g["gated_ns_min"] > 2.0 * NS_MIN_CM2 and g["gated_eeff_raw_min"] > 2e-2
    assert point["picard"]["iterations"] >= 2 and point["picard"]["converged"]


def test_mu_x_matches_finite_difference(cr, point):
    mod, cm = cr.ev.mod, cr.cmap
    X = point["state"]["X"]
    mu_x = point["state"]["mu_x"]
    rng = np.random.default_rng(7)
    d = rng.standard_normal(len(X))
    # restrict to the node set that mu depends on (all others have exactly zero columns)
    assert mu_x.nnz > 0

    def mu_of(Xp):
        n = mod.densities(Xp)[0]
        res, _ = cm.column_mobility(cr.P, n)
        return cm.edge_mobility(res["mu"])[cm.edges]

    for h in (1e-4, 1e-5):
        fd = (mu_of(X + h * d) - mu_of(X - h * d)) / (2 * h)
        an = mu_x @ d
        assert np.max(np.abs(fd - an)) <= 1e-6 * np.max(np.abs(an))
    # node sets: only n at the nodes within d_ch of the channel columns
    used = np.unique(mu_x.nonzero()[1] // 3)
    assert np.all(cm.in_ch[used])


def test_f_mu_matches_residual_difference(cr, point):
    mod, cm = cr.ev.mod, cr.cmap
    X = point["state"]["X"]
    mu0 = point["picard"]["mu_edge"].copy()
    mod.mu_edge = mu0.copy()  # the converged channel mobility
    f_mu = fmu_matrix(mod, X, cm.edges)
    rng = np.random.default_rng(3)
    dmu = rng.standard_normal(len(cm.edges)) * 1e-3 * np.abs(mu0[cm.edges])
    mod.mu_edge = mu0.copy()
    mu1 = mu0.copy()
    mu1[cm.edges] += dmu
    mod.mu_edge = mu1
    r1 = mod._residual_core(X)
    mu2 = mu0.copy()
    mu2[cm.edges] -= dmu
    mod.mu_edge = mu2
    r2 = mod._residual_core(X)
    mod.mu_edge = cm.mu_bulk.copy()
    fd = (r1 - r2) / 2.0
    an = f_mu @ dmu
    assert np.max(np.abs(an)) > 0
    assert np.max(np.abs(fd - an)) <= 1e-4 * np.max(np.abs(an))  # O(dmu^2) curvature


def test_ron_gradient_exact_vs_fd_and_lagged_differs(cr, point):
    g, gl = point["grad"], point["grad_lagged"]
    fd = np.zeros(5)
    for k in range(5):
        u = np.zeros(5)
        u[k] = 1.0
        fd[k] = richardson_fd(lambda t: cr.ron(t, grad=False)["R"], TH, u, 1e-3)[0]
    err = np.linalg.norm(g - fd) / np.linalg.norm(fd)
    err_l = np.linalg.norm(gl - fd) / np.linalg.norm(fd)
    assert err < 1e-6
    assert err_l > 1e-3 and err_l > 1e3 * err  # the F_mu mu_X term matters, and FD resolves it
    assert all(r["converged"] for r in cr.picard_log)


def test_guard_flags_a_sign_change(cr):
    """J4 guard: a theta that flips the surface sign of N_D - N_A is flagged (here via sign0)."""
    r = cr.ron(np.zeros(5), grad=False)
    assert r["guard"]["ok"] and r["guard"]["sign_changes"] == 0
    saved = cr.sign0.copy()
    try:
        cr.sign0 = -saved
        r = cr.ron(np.zeros(5), grad=False)
        assert not r["guard"]["ok"] and r["guard"]["sign_changes"] > 0
    finally:
        cr.sign0 = saved


def test_classification_constraint_gradient_and_sign():
    from widegap.adjoint_channel import ClassificationConstraint

    sc = SiC4H(TracedRegistry.from_json())
    ev = OxideAdjoint(sc, DESIGN, 0)
    cc = ClassificationConstraint([ev], margin=0.02)
    assert cc.n >= 2
    th = np.array([0.1, -0.05, 0.03, 0.02, 0.1])
    c0 = cc.values(np.zeros(5))
    assert np.all(c0 > 0.0) and np.all(c0 < 1.0 - 0.02 + 1e-12)
    jac = cc.jac(th)
    fd = np.zeros_like(jac)
    for k in range(5):
        u = np.zeros(5)
        u[k] = 1.0
        h = 1e-5
        fd[:, k] = (cc.values(th + h * u) - cc.values(th - h * u)) / (2 * h)
    assert np.linalg.norm(jac - fd) <= 1e-8 * np.linalg.norm(fd)
    # raising the JFET level lowers the margin of the corner column
    assert cc.values(np.array([0, 0, 0, 0, 1.0])).min() < c0.min()
