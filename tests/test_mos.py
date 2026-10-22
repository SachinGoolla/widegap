"""Fast tests for the stage-8.1 MOS machinery (regions, region-aware QF-DD, 1D reference)."""

import math

import numpy as np
from scipy.integrate import quad

from widegap import mos_reference as ref
from widegap.fv.drift_diffusion import DDParams
from widegap.fv2d.mesh import TensorMesh
from widegap.fv2d.mos import Dirichlet, MOSDriftDiffusion, continue_bias
from widegap.fv2d.regions import Regions

DP = DDParams(delta=1e-24, recombination=False)


def test_qs_first_integral_vs_quadrature():
    # (1/2) Q_s^2 = -int_0^u rho du',  rho = p0 (e^-u - 1) - n0 (e^u - 1)
    for a, delta, us in (
        (10.0, 1e-24, (-5, -0.3, 0.02, 1, 7, 40)),
        (2.0, 0.3, (-4, -1, 0.5, 3, 6)),
    ):
        _, p0, n0 = ref.bulk(a, delta)
        assert math.isclose(p0 - n0, a, rel_tol=1e-14) and math.isclose(
            p0 * n0, delta**2, rel_tol=1e-12
        )
        for u in us:
            integral, _ = quad(
                lambda t, p0=p0, n0=n0: -(p0 * math.expm1(-t) - n0 * math.expm1(t)), 0.0, u
            )
            assert math.isclose(0.5 * ref.Qs(u, a, delta) ** 2, integral, rel_tol=1e-8)
            assert ref.Qs(u, a, delta) * u < 0  # accumulation +, depletion/inversion -


def test_reference_root_find_roundtrip_and_overflow_safe():
    a, delta, tox, er = 10.0, 1e-24, 0.5, 0.38
    for u in (-5.0, -0.01, 0.0, 0.3, 30.0, 2 * 57.56 + 10):
        v = ref.gate_voltage(u, a, delta, tox, er)
        assert math.isfinite(v)
        assert abs(ref.psi_s_of_vg(v, a, delta, tox, er) - u) < 1e-9


def _small_mos_mesh(nx=3, depth=3.0, tox=0.5):
    ys = np.r_[np.linspace(-depth, -0.3, 4), np.linspace(-0.3, 0.0, 5)[1:], [0.25, tox]]
    m = TensorMesh(np.linspace(0.0, 1.0, nx + 1), ys)
    reg = Regions.from_function(m, lambda xc, yc: yc < 0)
    return m, reg


def test_region_masks():
    m, reg = _small_mos_mesh()
    y = m.node_xy[:, 1]
    assert np.array_equal(reg.semi_node, y <= 0.0)  # interface nodes touch a semiconductor cell
    assert np.array_equal(reg.oxide_node, y > 0.0)
    # semiconductor area is conserved; interface edges (both nodes semi, cells oxide) carry no flux
    assert np.isclose(reg.vol_s.sum(), 1.0 * 3.0)
    assert np.all(reg.vol_s[reg.oxide_node] == 0.0)
    horiz_if = (m.edge_dir == 0) & (m.node_xy[m.ea, 1] == 0.0)
    assert horiz_if.any()
    # the interface x-edge sees only the semiconductor cell below: half of the full face
    assert np.allclose(reg.w_s[horiz_if], m.w[horiz_if] - m.wc[horiz_if, 1])
    ox_edge = (y[m.ea] > 0.0) | (y[m.eb] > 0.0)
    assert np.all(reg.w_s[ox_edge & (m.edge_dir == 1)] == 0.0)
    # oxide-only nodes: identity rows for phi_n, phi_p in the Jacobian
    c = np.where(reg.semi_node, -10.0, 0.0)
    mod = MOSDriftDiffusion(reg, c, DP, np.log(1e-24))
    gate = np.flatnonzero(y == y.max())
    bc = Dirichlet().ohmic(c, 1e-24, np.flatnonzero(y == y.min())).gate(gate, 3.0, 0.0)
    X = mod.initial_guess(bc)
    rows, vals = mod._all_bc(bc)
    J = mod.jacobian_bc(X, rows).tocsr()
    for k in np.flatnonzero(reg.oxide_node):
        for v in (1, 2):
            r = J.getrow(3 * k + v)
            assert r.nnz == 1 and r[0, 3 * k + v] == 1.0
    # nothing flows or is stored in the oxide: n = p = 0 there
    sol = mod.solve_bc(X, mod.equilibrium_bc(bc))
    assert sol.report.converged
    assert np.all(sol.n[reg.oxide_node] == 0) and np.all(sol.p[reg.oxide_node] == 0)


def test_jacobian_vs_fd_with_oxide():
    m, reg = _small_mos_mesh(nx=2)
    delta = 1e-6
    cell = m.cell_centres
    c = reg.semi_average(np.where(cell[:, 0] < 0.5, 8.0, -3.0))
    mod = MOSDriftDiffusion(reg, c, DDParams(delta=delta, recombination=False), np.log(delta))
    xy = m.node_xy
    bc = (
        Dirichlet()
        .ohmic(c, delta, np.flatnonzero((xy[:, 1] == xy[:, 1].min()) & (xy[:, 0] == 0.0)), 0.2)
        .ohmic(c, delta, np.flatnonzero((xy[:, 1] == xy[:, 1].min()) & (xy[:, 0] == 1.0)), 0.0)
        .gate(np.flatnonzero(xy[:, 1] == xy[:, 1].max()), 2.0, 0.0)
    )
    rng = np.random.default_rng(3)
    X = mod.initial_guess(bc) + 0.3 * rng.standard_normal(3 * m.n_nodes)
    rows, vals = mod._all_bc(bc)
    X[rows] = vals
    Jan = mod.jacobian_bc(X, rows).toarray()
    Jfd = mod.fd_jacobian(X, bc, h=1e-6).toarray()
    err = np.abs(Jan - Jfd).max(axis=1) / np.abs(Jan).max(axis=1)
    assert err.max() < 1e-6


def _cap_solutions(mesh, bottom, top, c, reg, vgs):
    mod = MOSDriftDiffusion(reg, c, DP, np.log(1e-24))

    def bc_of(p):
        return Dirichlet().ohmic(c, 1e-24, bottom).gate(top, p[0], 0.0)

    X0 = mod.initial_guess(bc_of([vgs[0]]))
    s0 = mod.solve_bc(X0, mod.equilibrium_bc(bc_of([vgs[0]])))
    out = continue_bias(mod, s0.X, bc_of, [vgs[0]], [[v] for v in vgs[1:]], equilibrium=True)
    return [s0] + [o[1] for o in out]


def test_2d_capacitor_column_equals_1d():
    y = np.r_[np.linspace(-6.0, -1.0, 6), -0.4, -0.12, -0.03, 0.0, np.linspace(0.0, 0.5, 3)[1:]]
    m1 = TensorMesh(y)
    reg1 = Regions.from_function(m1, lambda yc: yc < 0)
    c1 = np.where(reg1.semi_node, -10.0, 0.0)
    m2 = TensorMesh(np.array([0.0, 1.0, 2.5]), y)
    reg2 = Regions.from_function(m2, lambda xc, yc: yc < 0)
    c2 = np.where(reg2.semi_node, -10.0, 0.0)
    vgs = [ref.gate_voltage(u, 10.0, 1e-24, 0.5, 0.38) for u in (0.0, 5.0, 60.0, 119.0)]
    s1 = _cap_solutions(m1, [0], [m1.n_nodes - 1], c1, reg1, vgs)
    ny1 = m2.ny + 1
    bottom = [m2.node(i, 0) for i in range(m2.nx + 1)]
    top = [m2.node(i, m2.ny) for i in range(m2.nx + 1)]
    s2 = _cap_solutions(m2, bottom, top, c2, reg2, vgs)
    for a, b in zip(s1, s2, strict=True):
        assert a.report.converged and b.report.converged
        for i in range(m2.nx + 1):
            assert np.abs(b.psi[i * ny1 : (i + 1) * ny1] - a.psi).max() < 1e-10
