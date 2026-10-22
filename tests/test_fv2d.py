"""Fast unit tests for the stage-8 edge-based 2D finite-volume core."""

import numpy as np

from widegap.fv.drift_diffusion import DDParams
from widegap.fv.qf_drift_diffusion import QFDriftDiffusion1D
from widegap.fv2d.core import edge_conductance, stiffness
from widegap.fv2d.mesh import TensorMesh
from widegap.fv2d.poisson import solve_poisson
from widegap.fv2d.qf_dd import QFDriftDiffusion, pack


def _graded(n, r):
    h = r ** np.arange(n)
    return np.r_[0.0, np.cumsum(h)] / h.sum()


def test_volumes_and_conductances():
    m = TensorMesh(_graded(7, 1.1), _graded(5, 0.9))
    assert np.isclose(m.vol.sum(), 1.0, rtol=1e-14)
    assert np.isclose(m.cell_size.sum(), 1.0, rtol=1e-14)
    # unit coefficient: row sums of L vanish and G_e = w_e / h_e
    L = stiffness(m, 1.0)
    assert np.abs(L @ np.ones(m.n_nodes)).max() < 1e-13
    assert np.allclose(edge_conductance(m, 1.0), m.w / m.h)
    # CV average of a constant is the constant; of a linear field is exact
    assert np.allclose(m.cv_average(np.full(m.n_cells, 3.0)), 3.0)
    m1 = TensorMesh(np.array([0.0, 0.3, 1.0, 1.2]))
    assert np.isclose(m1.vol.sum(), 1.2)


def test_uniform_poisson_is_five_point_stencil():
    n = 8
    m = TensorMesh(np.linspace(0, 1, n + 1), np.linspace(0, 1, n + 1))
    rng = np.random.default_rng(1)
    u = rng.standard_normal(m.n_nodes)
    Lu = (stiffness(m, 1.0) @ u).reshape(n + 1, n + 1)
    U = u.reshape(n + 1, n + 1)
    st = 4 * U[1:-1, 1:-1] - U[2:, 1:-1] - U[:-2, 1:-1] - U[1:-1, 2:] - U[1:-1, :-2]
    assert np.allclose(Lu[1:-1, 1:-1], st)  # dual area h^2 cancels the 1/h^2
    # solve: exact for quadratic u = x^2 + y^2  (f = -4)
    ue = m.node_xy[:, 0] ** 2 + m.node_xy[:, 1] ** 2
    bn = m.boundary()
    sol = solve_poisson(m, 1.0, np.full(m.n_nodes, -4.0), bn, ue[bn])
    assert np.abs(sol - ue).max() < 1e-12


def _random_state(m, rng, ld):
    xy = m.node_xy
    psi = 5 * np.sin(2 * xy[:, 0] + 0.3) + rng.standard_normal(m.n_nodes) * 0.5
    pn = rng.standard_normal(m.n_nodes)
    pp = rng.standard_normal(m.n_nodes)
    return pack(psi, pn, pp)


def test_jacobian_vs_fd_2d():
    rng = np.random.default_rng(3)
    m = TensorMesh(_graded(5, 1.2), _graded(4, 0.9))
    for ld, rec in ((np.log(1e-3), True), (np.log(1e-24), True), (np.log(1e-24), False)):
        prm = DDParams(delta=float(np.exp(ld)), tau_n=3.0, tau_p=2.0, recombination=rec)
        c = rng.standard_normal(m.n_nodes)
        dd = QFDriftDiffusion(m, c, prm, ld)
        dn = m.left
        dv = rng.standard_normal((len(dn), 3))
        X = _random_state(m, rng, ld)
        J = dd.jacobian(X, dn)
        v = rng.standard_normal(len(X))
        for k in range(3):
            v[3 * dn + k] = 0.0
        e = 1e-7
        fd = (dd.residual(X + e * v, dn, dv) - dd.residual(X - e * v, dn, dv)) / (2 * e)
        den = abs(J) @ np.abs(v)
        mask = den > 0
        assert np.max(np.abs(J @ v - fd)[mask] / den[mask]) < 1e-6


def test_1d_through_core_equals_stage71_residual():
    rng = np.random.default_rng(5)
    x = np.cumsum(np.r_[0.0, rng.uniform(0.1, 0.5, 30)])
    ld = np.log(1e-24)
    prm = DDParams(delta=float(np.exp(ld)), tau_n=1e3, tau_p=1e3, recombination=True)
    c = rng.standard_normal(len(x))
    old = QFDriftDiffusion1D(x, c, prm, ld)
    new = QFDriftDiffusion(TensorMesh(x), c, prm, ld)
    X = pack(
        np.cumsum(rng.standard_normal(len(x))),
        rng.standard_normal(len(x)),
        rng.standard_normal(len(x)),
    )
    bc = [tuple(X[:3] + 0.1), tuple(X[-3:] - 0.2)]
    dn = np.array([0, len(x) - 1])
    dv = np.array([bc[0], bc[1]])
    Fo, Fn = old.residual(X, bc), new.residual(X, dn, dv)
    assert np.allclose(Fo, Fn, rtol=1e-13, atol=1e-300)
    Jo, Jn = old.jacobian(X), new.jacobian(X, dn)
    assert abs(Jo - Jn).max() <= 1e-13 * abs(Jo).max()
