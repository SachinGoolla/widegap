"""Fast tests for the stage-9a physics models (closed forms, derivatives, small solves)."""

import math

import numpy as np
from scipy.integrate import quad

from widegap import mos_reference as ref
from widegap.fv.drift_diffusion import DDParams
from widegap.fv.ionized_equilibrium import IonizedJunction
from widegap.fv2d.mesh import TensorMesh
from widegap.fv2d.mos import Dirichlet, MOSDriftDiffusion, continue_bias
from widegap.fv2d.physics import (
    InterfaceCharge,
    Ionization,
    Mobility,
    PhysicsMOSDD,
    ct_factor,
    gate_voltage_it,
    interface_length,
    neutral_psi,
    psi_s_of_vg_it,
    trap_charge,
)
from widegap.fv2d.regions import Regions
from widegap.materials.ionization import Level


def _all_semi(mesh):
    return Regions.from_function(mesh, lambda *c: np.ones_like(c[0], dtype=bool))


def _fd_dq(s, dit, e1, e2, h=1e-6):
    return (trap_charge(s + h, dit, e1, e2)[0] - trap_charge(s - h, dit, e1, e2)[0]) / (2 * h)


def test_trap_charge_closed_form_vs_quadrature():
    for e1, e2 in ((-20.0, 20.0), (-25.0, 15.0), (-5.0, -1.0), (2.0, 9.0), (-3.0, 0.0)):
        for s in (-60.0, -22.0, -3.0, 0.0, 0.4, 7.0, 18.0, 55.0):
            f = lambda e, s=s: 1.0 / (1.0 + math.exp(min(e - s, 700.0)))  # noqa: E731
            acc = (
                quad(f, max(e1, 0.0), e2, points=[s] if max(e1, 0) < s < e2 else None)[0]
                if e2 > 0
                else 0.0
            )
            don = quad(lambda e, f=f: 1.0 - f(e), e1, min(e2, 0.0))[0] if e1 < 0 else 0.0
            q_ref = 0.7 * (don - acc)
            q, dq = trap_charge(s, 0.7, e1, e2)
            assert math.isclose(q, q_ref, rel_tol=1e-9, abs_tol=1e-12)
            assert math.isclose(float(dq), _fd_dq(s, 0.7, e1, e2), rel_tol=1e-6, abs_tol=1e-9)
            assert dq <= 1e-12 * 0.7


def test_trap_charge_limits():
    # all empty far below, all full far above: donors +D*|E1|, acceptors -D*E2
    assert math.isclose(float(trap_charge(-200.0, 0.3, -20.0, 10.0)[0]), 0.3 * 20.0, rel_tol=1e-12)
    assert math.isclose(float(trap_charge(200.0, 0.3, -20.0, 10.0)[0]), -0.3 * 10.0, rel_tol=1e-12)


def test_mobility_factor_derivative_and_limits():
    for beta in (1.0, 1.5, 2.0):
        for d in (-3.0, -0.2, 0.7, 4.0):
            m, k = ct_factor(d, 0.5, 1.3, 2.0, beta)
            h = 1e-6
            fd = (
                ct_factor(d + h, 0.5, 1.3, 2.0, beta)[0] - ct_factor(d - h, 0.5, 1.3, 2.0, beta)[0]
            ) / (2 * h)
            assert 0 < m <= 1 and math.isclose(float(k), float(fd), rel_tol=1e-6)
        m0, k0 = ct_factor(0.0, 0.5, 1.3, 2.0, beta)
        assert m0 == 1.0 and abs(k0) <= 1e-300 or beta == 1.0 and k0 == 0.0
    m, k = ct_factor(np.array([5.0]), 1.0, 1.0, 1e30, 2.0)
    assert m[0] == 1.0 and abs(k[0]) < 1e-50
    # saturation: mu F -> vsat
    m, _ = ct_factor(1e6, 1.0, 1.0, 0.5, 2.0)
    assert math.isclose(m * 1e6, 0.5, rel_tol=1e-9)


def test_interface_length_1d_and_2d():
    ys = np.r_[np.linspace(-1, 0, 3), [0.5]]
    m1 = TensorMesh(ys)
    reg1 = Regions.from_function(m1, lambda yc: yc < 0)
    L = interface_length(reg1)
    assert L.sum() == 1.0 and L[2] == 1.0
    m2 = TensorMesh(np.array([0.0, 0.5, 2.0]), ys)
    reg2 = Regions.from_function(m2, lambda xc, yc: yc < 0)
    L2 = interface_length(reg2)
    assert math.isclose(L2.sum(), 2.0)
    assert np.all(L2[m2.node_xy[:, 1] != 0.0] == 0.0)
    assert math.isclose(L2[m2.node(1, 2)], 1.0)


def test_a1_equilibrium_matches_ionized_junction():
    lev = (Level(1.0, -58.5),)
    t = np.linspace(-1, 1, 81)
    x = 12.0 * np.sinh(3 * t) / np.sinh(3.0)
    nd = np.where(x > 0, 10.0, 0.0)
    na = np.where(x < 0, 10.0, 0.0)
    nd[x == 0] = na[x == 0] = 5.0
    mesh = TensorMesh(x)
    dd = PhysicsMOSDD(
        _all_semi(mesh), nd - na, DDParams(1e-24, recombination=False), math.log(1e-24),
        ionization=Ionization(nd, na, lev, lev),
    )  # fmt: skip
    pl, pr = neutral_psi(0.0, 10.0, 1e-24, lev, lev), neutral_psi(10.0, 0.0, 1e-24, lev, lev)
    bc = dd.equilibrium_bc(Dirichlet().psi([0, len(x) - 1], [pl, pr]))
    sol = dd.solve_bc(dd.initial_guess(bc), bc)
    psi, rep = IonizedJunction(x, nd, na, 1e-24, lev, lev).solve(np.where(x > 0, pr, pl), pl, pr)
    assert sol.report.converged and rep.converged
    assert np.abs(sol.psi - psi).max() <= 1e-10


def _small_mos(delta, rng):
    ys = np.r_[np.linspace(-3, -0.3, 4), np.linspace(-0.3, 0.0, 5)[1:], [0.25, 0.5]]
    m = TensorMesh(np.linspace(0, 1, 4), ys)
    reg = Regions.from_function(m, lambda xc, yc: yc < 0)
    nd = np.where(reg.semi_node, rng.uniform(1, 5, m.n_nodes), 0.0)
    na = np.where(reg.semi_node, rng.uniform(1, 5, m.n_nodes), 0.0)
    lg = math.log(delta)
    ion = Ionization(nd, na, (Level(0.6, lg + 6), Level(0.4, lg + 14)), (Level(1.0, lg + 4),))
    dd = PhysicsMOSDD(
        reg, nd - na, DDParams(delta, recombination=False), lg, ionization=ion,
        mobility=Mobility(0.7, 0.4, 2.0, 1.0), interface=InterfaceCharge(1.0, 0.1, -10.0, 15.0),
    )  # fmt: skip
    return m, reg, dd


def test_jacobian_all_models_vs_fd():
    rng = np.random.default_rng(3)
    for delta, scale in ((1e-24, 55.0), (1e-3, 3.0)):
        m, reg, dd = _small_mos(delta, rng)
        for _ in range(2):
            X = np.zeros(3 * m.n_nodes)
            sm = reg.semi_node
            X[0::3] = np.where(
                sm, rng.uniform(-scale, scale, m.n_nodes), rng.uniform(-3, 3, m.n_nodes)
            )
            X[1::3] = np.where(sm, rng.uniform(-1.5, 1.5, m.n_nodes), 0.0)
            X[2::3] = np.where(sm, rng.uniform(-1.5, 1.5, m.n_nodes), 0.0)
            bc = Dirichlet()
            rows, _ = dd._all_bc(bc)
            Ja = dd.jacobian_bc(X, rows).toarray()
            Jf = dd.fd_jacobian(X, bc, h=1e-6).toarray()
            err = np.abs(Jf - Ja).max(axis=1) / np.abs(Ja).max(axis=1)
            assert err.max() <= 1e-6


def test_vsat_infinite_reproduces_constant_mobility():
    rng = np.random.default_rng(4)
    m, reg, _ = _small_mos(1e-3, rng)
    c = np.where(reg.semi_node, -2.0, 0.0)
    dp = DDParams(1e-3, recombination=False)
    base = MOSDriftDiffusion(reg, c, dp, math.log(1e-3))
    new = PhysicsMOSDD(reg, c, dp, math.log(1e-3), mobility=Mobility(1e30, 1e30, 2.0, 1.0))
    y = m.node_xy[:, 1]
    bot, gate = np.flatnonzero(y == y.min()), np.flatnonzero(y == y.max())
    left, right = (
        np.flatnonzero((m.node_xy[:, 0] == 0) & (y == 0.0)),
        np.flatnonzero((m.node_xy[:, 0] == 1) & (y == 0.0)),
    )
    cs = []
    for mod in (base, new):

        def bc_of(p, mod=mod):
            return (Dirichlet().ohmic(c, 1e-3, bot, 0.0).ohmic(c, 1e-3, left, 0.0)
                    .ohmic(c, 1e-3, right, p[1]).gate(gate, p[0], 0.0))  # fmt: skip

        eq = mod.solve_bc(mod.initial_guess(bc_of([0, 0])), mod.equilibrium_bc(bc_of([0, 0])))
        out = continue_bias(mod, eq.X, bc_of, [0.0, 0.0], [[3.0, 4.0]])
        assert out[-1][1].report.converged
        cs.append(out[-1][1].j)
    gated = np.abs(cs[0]) >= 1e-6
    assert gated.any()
    assert np.max(np.abs(cs[1] - cs[0])[gated] / np.abs(cs[0])[gated]) <= 1e-12


def test_uniform_bar_saturating_current():
    nd0, dl, L, vs = 10.0, 1e-3, 40.0, 0.5
    mesh = TensorMesh(np.linspace(0.0, L, 41))
    c = np.full(mesh.n_nodes, nd0)
    dp = DDParams(dl, recombination=False)
    dd = PhysicsMOSDD(_all_semi(mesh), c, dp, math.log(dl), mobility=Mobility(vs, vs, 2.0, 1.0))
    xc = 0.5 * (mesh.x[:-1] + mesh.x[1:])
    inner = (xc >= 5) & (xc <= L - 5)
    for fct in (1.0, 20.0):
        v = fct * vs * L
        bc = Dirichlet().ohmic(c, dl, [0], 0.0).ohmic(c, dl, [mesh.n_nodes - 1], v)
        X0 = np.zeros(3 * mesh.n_nodes)
        X0[0::3] = np.arcsinh(nd0 / (2 * dl)) + v * mesh.x / L
        X0[1::3] = X0[2::3] = v * mesh.x / L
        sol = dd.solve_bc(X0, bc)
        assert sol.report.converged
        F = v / L
        j = nd0 * F * (1.0 + (F / vs) ** 2) ** -0.5
        assert np.max(np.abs(-sol.j[inner] - j)) / j <= 1e-6


def test_mos_cap_with_interface_charge_matches_reference():
    itf = InterfaceCharge(3.0, 0.2, -25.0, 15.0)
    a, delta, tox, er = 10.0, 1e-24, 0.5, 0.38
    # reference root-find round trip, and a Q_it-dependent shift relative to Q_f only
    for u in (-3.0, 0.0, 30.0, 60.0, 100.0):
        v = gate_voltage_it(u, a, delta, tox, er, 0.0, itf)
        assert abs(psi_s_of_vg_it(v, a, delta, tox, er, 0.0, itf) - u) < 1e-9
    h = 1.2 ** np.arange(30)[::-1]
    h = h * 30.0 / h.sum()
    ys = np.r_[-30.0, -30.0 + np.cumsum(h)]
    ys[-1] = 0.0
    y = np.r_[ys, np.linspace(0.0, tox, 5)[1:]]
    m = TensorMesh(y)
    reg = Regions.from_function(m, lambda *c: c[-1] < 0)
    c = np.where(reg.semi_node, -a, 0.0)
    mod = PhysicsMOSDD(reg, c, DDParams(delta, recombination=False), math.log(delta), interface=itf)
    psi_b = ref.bulk(a, delta)[0]
    isurf = int(np.argmin(np.abs(y)))
    bc_of = lambda p: Dirichlet().ohmic(c, delta, [0], 0.0).gate([m.n_nodes - 1], p[0], 0.0)  # noqa: E731
    targets = [0.0, 20.0, 60.0, 115.0]
    vgs = [gate_voltage_it(u, a, delta, tox, er, 0.0, itf) for u in targets]
    X0 = mod.initial_guess(bc_of([vgs[0]]))
    s0 = mod.solve_bc(X0, mod.equilibrium_bc(bc_of([vgs[0]])))
    out = continue_bias(mod, s0.X, bc_of, [vgs[0]], [[v] for v in vgs[1:]], equilibrium=True)
    for u, (_, s, _, _) in zip(targets[1:], out, strict=True):
        assert s.report.converged
        assert (
            abs(s.psi[isurf] - psi_b - u) <= 5e-3 * u
        )  # single mesh: loose; benchmark A5 uses Richardson
