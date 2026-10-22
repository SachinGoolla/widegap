"""Fast tests of the stage-13 full-physics adjoint (widegap.adjoint_sic) on a small 4H-SiC
DMOSFET mesh (187 nodes, physical units, sourced parameters): analytic pieces against finite
differences and adjoint gradients of R_on,sp and of the fixed-cut avalanche BV against
Richardson-extrapolated central differences (eps = 1e-3, 5e-4).  < 60 s."""

import numpy as np
import pytest

from widegap.adjoint_sic import SiCAdjoint, ionization_n_grad, richardson_fd
from widegap.fv.avalanche import ionization_integrals
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


@pytest.fixture(scope="module")
def sc():
    return SiC4H(TracedRegistry.from_json())


@pytest.fixture(scope="module")
def adj(sc):
    return SiCAdjoint(sc, DESIGN, 0)


def _unit(seed=0):
    u = np.random.default_rng(seed).standard_normal(5)
    return u / np.linalg.norm(u)


def test_ionization_integral_value_and_field_gradient(sc):
    rng = np.random.default_rng(3)
    x = np.cumsum(rng.uniform(0.5, 3.0, 60))
    x = np.r_[0.0, x]
    for ch in (sc.niwa(), sc.hatakeyama()):
        e = rng.uniform(40.0, 160.0, 60) * rng.choice([-1.0, 1.0], 60)  # sign-changing fields
        val, grad = ionization_n_grad(x, e, ch)
        assert val == pytest.approx(ionization_integrals(x, e, ch)[0], rel=1e-14)
        for k in rng.choice(60, 8, replace=False):
            h = 1e-6 * abs(e[k])
            ep, em = e.copy(), e.copy()
            ep[k] += h
            em[k] -= h
            fd = (ionization_integrals(x, ep, ch)[0] - ionization_integrals(x, em, ch)[0]) / (2 * h)
            assert grad[k] == pytest.approx(fd, rel=1e-7, abs=1e-12 * abs(val))


def test_theta0_model_equals_stage12_and_mobility_derivative(sc, adj):
    b = adj.base
    adj.set_theta(np.zeros(5))
    X = b.equilibrium().X
    rows, vals = b.model._all_bc(b.bc(0, 0, 0))
    assert np.array_equal(b.model.residual_bc(X, rows, vals), adj.mod.residual_bc(X, rows, vals))
    assert np.array_equal(b.model.mu_edge, adj.mod.mu_edge)
    n = np.array([1e15, 5e16, 3e17, 2e18])
    for o in ("perp", "par"):
        mu, dmu = adj._mobility_cm2_and_deriv(n, o)
        assert np.allclose(mu, sc.mobility_cm2(n, "n", o), rtol=1e-14)
        h = 1e-6 * n
        fd = (sc.mobility_cm2(n + h, "n", o) - sc.mobility_cm2(n - h, "n", o)) / (2 * h)
        assert np.allclose(dmu, fd, rtol=1e-7)


def test_dF_dtheta_analytic_vs_residual_central_difference(adj):
    rng = np.random.default_rng(1)
    adj.set_theta(TH)
    X = adj.base.equilibrium().X
    semi = adj.mod.reg.semi_node
    X[1::3] += 0.3 * rng.standard_normal(len(X) // 3) * semi  # nonzero phi_n: mobility path
    X[0::3] += 0.05 * rng.standard_normal(len(X) // 3) * semi
    w = rng.standard_normal(len(X))
    u = _unit(2)
    g = adj.dnd_dtheta.T @ adj.dF_dnd(w, X)

    def f(t):
        adj.set_theta(t)
        return w @ adj.mod._residual_core(X)

    eps = 1e-6
    fd = (f(TH + eps * u) - f(TH - eps * u)) / (2 * eps)
    adj.set_theta(TH)
    assert abs(g @ u - fd) <= 1e-7 * abs(fd)


def test_ron_adjoint_vs_richardson_fd_and_drain_current_identity(adj):
    r = adj.ron(TH)
    # I_drain = sum of the drain-node continuity rows of the core residual (identity used for dI/dX)
    F = adj.mod._residual_core(r["X"])
    dn = adj.base.dev.nodes["drain"]
    assert F[3 * dn + 1].sum() + F[3 * dn + 2].sum() == pytest.approx(r["I_scaled"], rel=1e-9)
    for seed in (0, 1):
        u = _unit(seed)
        rich, _, _ = richardson_fd(lambda t: adj.ron(t, grad=False)["R"], TH, u, 1e-3)
        assert abs(r["grad"] @ u - rich) <= 1e-6 * abs(rich)


def test_bv_adjoint_vs_richardson_fd_and_root_property(sc, adj):
    b = adj.bv(TH)
    st = b["state"]
    # dI_n/dV_D: adjoint (lam^T dvals) equals the tangent (dI/dX . dX/dV)
    assert b["di_dv_adjoint"] == pytest.approx(b["di_dv_tangent"], rel=1e-9)
    # independent root check: re-solve at the BV, library I_n along the fixed column = 1
    rows, vals = adj._off_bc(st.v_s)
    s = adj.mod.solve_bc(st.X, _bc(rows, vals), max_step=20.0)
    assert s.report.converged
    nodes = adj.cut_nodes
    depth = -adj.mesh.node_xy[nodes, 1]
    e = -np.diff(s.X[0::3][nodes]) / np.diff(depth)
    assert ionization_integrals(depth, e, adj.ch)[0] == pytest.approx(1.0, abs=1e-7)
    assert adj.col_x == pytest.approx(adj.prm.x_bc)
    for seed in (0, 1):
        u = _unit(seed)
        rich, _, _ = richardson_fd(lambda t: adj.bv(t, grad=False)["BV"], TH, u, 1e-3)
        assert abs(b["grad"] @ u - rich) <= 1e-6 * abs(rich)
    # cold march and warm Newton agree
    adj.set_theta(TH)
    cold = adj.bv(TH, grad=False, warm=False)["BV"]
    assert cold == pytest.approx(b["BV"], rel=1e-9)


def _bc(rows, vals):
    from widegap.adjoint_sic import _AffineBC

    return _AffineBC(rows, vals)


def test_design_variable_dependence_is_the_expected_physics(adj):
    """Sanity of signs (not tuned numbers): more drift donors -> lower BV, lower R_on."""
    th_hi = TH.copy()
    th_hi[0] += 0.3
    assert adj.bv(th_hi, grad=False)["BV"] < adj.bv(TH, grad=False)["BV"]
    assert adj.ron(th_hi, grad=False)["R"] < adj.ron(TH, grad=False)["R"]
    th_j = TH.copy()
    th_j[4] += 0.3
    assert adj.ron(th_j, grad=False)["R"] < adj.ron(TH, grad=False)["R"]
