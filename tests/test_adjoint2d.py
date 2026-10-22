"""Fast tests of the stage-10 2D adjoint (tiny DMOSFET mesh): adjoint vs central finite
differences for I_D and for B_s at fixed bias, the implicit-function dBV_s/dC pieces
(dB_s/dV_D), and the analytic dC/dtheta of the device builder."""

import numpy as np
import pytest

from widegap.adjoint2d import (
    FieldIndex,
    adjoint_solve,
    breakdown_voltage,
    bv_gradient_c,
    contact_current_grad_x,
    dj_dc,
    drain_current_grad_c,
)
from widegap.devices.dmosfet import DMOSFET, DMOSParams
from widegap.fv2d.mos import continue_bias

pytestmark = pytest.mark.filterwarnings(
    "ignore::RuntimeWarning", "ignore::scipy.sparse.linalg.MatrixRankWarning"
)

TINY = DMOSParams(
    n0_drift=16.0,
    n_jfet=18.0,
    nx=(2, 1, 2, 2, 2, 1),
    ny_drain=1,
    ny_drift=4,
    ny_body=4,
    ny_surf=4,
    ny_ox=1,
)
E_CRIT, KAPPA = 75.0, 200.0
TH0 = np.array([0.1, -0.05, 0.03, 0.0, 0.1])


@pytest.fixture(scope="module")
def dev():
    return DMOSFET(TINY)


def _eq(dev, th, X0=None):
    mod, c = dev.model(th)
    bc = dev.bc(c, 0.0, 0.0, 0.0)
    s = mod.solve_bc(mod.initial_guess(bc) if X0 is None else X0, mod.equilibrium_bc(bc))
    assert s.report.converged
    return mod, c, s


def _state(dev, th, vg, vd, X0=None):
    mod, c, s0 = _eq(dev, th, X0)
    out = continue_bias(
        mod, s0.X, lambda p: dev.bc(c, p[0], 0.0, p[1]), [0, 0], [[0, vd], [vg, vd]]
    )
    s = out[-1][1]
    assert s.report.converged
    return mod, c, s


def test_doping_jacobian_matches_fd_and_is_semiconductor_only(dev):
    c0, dc = dev.doping(TH0)
    eps = 1e-6
    for k in range(5):
        e = np.zeros(5)
        e[k] = eps
        fd = (dev.doping(TH0 + e)[0] - dev.doping(TH0 - e)[0]) / (2 * eps)
        assert np.max(np.abs(fd - dc[:, k])) <= 1e-7 * max(1.0, np.max(np.abs(dc[:, k])))
    assert np.all(c0[dev.reg.oxide_node] == 0.0) and np.all(dc[dev.reg.oxide_node] == 0.0)
    # uniform drift profile: theta_0 only shifts log N_D
    assert np.isclose(dev.jfet_level(TH0), 18.0 * np.exp(0.1))
    y = np.array([-12.0, -16.0])
    assert np.allclose(dev.drift_profile(np.zeros(5), y), 16.0)
    nd = dev.drift_profile(TH0, y)
    jac = dev.drift_profile_jac(TH0, y)
    fd = (dev.drift_profile(TH0 + [1e-6, 0, 0, 0, 0], y) - nd) / 1e-6
    assert np.allclose(fd, jac[:, 0], rtol=1e-5)


def test_current_gradient_x_matches_fd(dev):
    mod, c, s = _state(dev, TH0, 200.0, 2.0)
    cur, g = contact_current_grad_x(mod, s.X, dev.nodes["drain"])
    assert np.isclose(cur, mod.contact_current(s, dev.nodes["drain"]), rtol=1e-13)
    rng = np.random.default_rng(3)
    u = rng.standard_normal(len(s.X))
    h = 1e-6
    semi = np.repeat(mod.reg.semi_node, 3)
    u[~semi] = 0.0

    def cur(X):
        sol = type(s)(X, *mod.split(X), s.n, s.p, *mod.fluxes(X), s.report)
        return mod.contact_current(sol, dev.nodes["drain"])

    fd = (cur(s.X + h * u) - cur(s.X - h * u)) / (2 * h)
    assert abs(g @ u - fd) <= 1e-6 * max(1.0, abs(fd))


def test_adjoint_drain_current_vs_fd(dev):
    mod, c, s = _state(dev, TH0, 200.0, 2.0)
    bc = dev.bc(c, 200.0, 0.0, 2.0)
    _, gc = drain_current_grad_c(mod, s.X, bc, dev.nodes["drain"])
    _, dcdth = dev.doping(TH0)
    g = dcdth.T @ gc

    def f(th):
        m2, c2, s2 = _state(dev, th, 200.0, 2.0, X0=s.X)
        return m2.contact_current(s2, dev.nodes["drain"])

    eps = 1e-3
    for seed in range(3):
        u = np.random.default_rng(seed).standard_normal(5)
        u /= np.linalg.norm(u)
        fd = (f(TH0 + eps * u) - f(TH0 - eps * u)) / (2 * eps)
        assert abs(g @ u - fd) <= 1e-5 * abs(fd)


def test_adjoint_smoothed_index_and_dbs_dvd_vs_fd(dev):
    vd = 120.0
    mod, c, s = _state(dev, TH0, 0.0, vd)
    idx = FieldIndex(E_CRIT, dev.blocking_edge_mask(), KAPPA)
    bc = dev.bc(c, 0.0, 0.0, vd)
    b0, gx = idx(mod, s.X, grad=True)
    lam = adjoint_solve(mod, s.X, bc, gx)
    _, dcdth = dev.doping(TH0)
    g = dcdth.T @ dj_dc(mod, lam, bc)

    def bs(th, v=vd):
        m2, c2, s2 = _state(dev, th, 0.0, v, X0=s.X)
        return idx(m2, s2.X)

    assert np.isclose(bs(TH0), b0, rtol=1e-12)
    eps = 1e-3
    for seed in range(3):
        u = np.random.default_rng(seed).standard_normal(5)
        u /= np.linalg.norm(u)
        fd = (bs(TH0 + eps * u) - bs(TH0 - eps * u)) / (2 * eps)
        assert abs(g @ u - fd) <= 1e-4 * abs(fd)
    # dB_s/dV_D from the drain Dirichlet rows
    _, dbdv, _ = bv_gradient_c(mod, s.X, bc, dev.nodes["drain"], idx)
    h = 1e-3
    fd = (bs(TH0, vd + h) - bs(TH0, vd - h)) / (2 * h)
    assert abs(dbdv - fd) <= 1e-5 * abs(fd)


def test_bv_implicit_function_gradient_vs_fd(dev):
    idx = FieldIndex(E_CRIT, dev.blocking_edge_mask(), KAPPA)

    def bv(th):
        mod, c, s = _eq(dev, th)
        return breakdown_voltage(mod, lambda v: dev.bc(c, 0.0, 0.0, v), idx, s.X, v_step=50.0), mod

    r, mod = bv(TH0)
    assert abs(idx(mod, r.X) - 1.0) <= 1e-9
    g_c, _, _ = bv_gradient_c(mod, r.X, r.bc, dev.nodes["drain"], idx)
    g = dev.doping(TH0)[1].T @ g_c
    eps = 1e-3
    u = np.random.default_rng(1).standard_normal(5)
    u /= np.linalg.norm(u)
    fd = (bv(TH0 + eps * u)[0].bv - bv(TH0 - eps * u)[0].bv) / (2 * eps)
    assert abs(g @ u - fd) <= 1e-4 * abs(fd)
