"""Stage-12 infrastructure tests (widegap.sic).  Fast (< 30 s); the SiC values come from the
registry, the checks are identities / independent re-computations, not tuned numbers."""

import copy
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

from widegap.fv.avalanche import (
    Chynoweth,
    avalanche_bv,
    ionization_integrals,
    uniform_field_integral,
)
from widegap.fv.drift_diffusion import DDParams
from widegap.fv.equilibrium import cv_average_piecewise
from widegap.fv.ionized_equilibrium import IonizedJunction
from widegap.fv2d.mesh import TensorMesh
from widegap.fv2d.mos import Dirichlet
from widegap.fv2d.physics import Ionization
from widegap.fv2d.regions import Regions
from widegap.materials.ionization import K_B, Q, neutral_potential
from widegap.materials.registry import ProvenanceError, Registry
from widegap.sic import (
    REGISTRY_PATH,
    AnisoMobilityDD,
    DMOSDesign,
    IonizedReverseDiode,
    SiC4H,
    SiCDMOS,
    TracedRegistry,
    build_pn_diode,
    tangent_predictor,
    triangular_breakdown,
)

REPO = Path(__file__).resolve().parents[1]
DESIGN = DMOSDesign(
    x_jfet=0.8, ch_len=0.8, src_w=1.0, bc_w=0.8, pad_w=0.4,
    d_src=0.3, d_body=0.8, w_drift=3.0, t_drain=0.3, t_ox=0.05,
    n_drift=5e16, n_jfet=1e17, n_body=3e17, n_src=1e19, n_drain=1e19,
    mesh_kw=dict(ny_surf=10, q_surf=0.6),
)  # fmt: skip


@pytest.fixture(scope="module")
def sc():
    return SiC4H(TracedRegistry.from_json())


# ----------------------------------------------------------------------------- S1 / mapping
def test_registry_reads_are_traced_and_complete(sc):
    sc.niwa(), sc.hatakeyama(), sc.vfb(3e17), sc.mobility_cm2(1e17, "n", "par")
    rows = sc.traceability()
    names = {r["name"] for r in rows}
    assert set(sc.reg.used) <= names
    for r in rows:
        if r["class"] not in ("constant", "spec", "structural"):
            assert r["source"].startswith(("doi:", "owner approval")), r
            assert sc.reg.usable(r["name"])
    SiCDMOS(sc, DESIGN)  # reads the owner-approved oxide permittivity
    assert any(
        r["name"] == "eps_r_SiO2_thermal" and "owner approval" in r["source"]
        for r in sc.traceability()
    )


def test_unusable_or_missing_entry_raises_no_fallback(tmp_path):
    data = json.loads(Path(REGISTRY_PATH).read_text())
    bad = copy.deepcopy(data)
    for p in bad["parameters"]:
        if p["name"] == "eps_r_parallel_c":
            p["access"] = "abstract-only"
    reg = TracedRegistry(Registry(bad["parameters"]))
    with pytest.raises(ProvenanceError):
        SiC4H(reg)
    miss = copy.deepcopy(data)
    miss["parameters"] = [p for p in miss["parameters"] if p["name"] != "oc_a_n"]
    sc2 = SiC4H(TracedRegistry(Registry(miss["parameters"])))
    with pytest.raises(KeyError):
        sc2.niwa()


def test_derived_statistics_equal_stage7_benchmark(sc):
    spec = importlib.util.spec_from_file_location("s7", REPO / "benchmarks" / "stage7" / "run.py")
    s7 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(s7)
    reg = Registry.from_json(REGISTRY_PATH)
    ni, nc, nv = s7.derived(reg, 300.0)
    assert sc.n_i == pytest.approx(ni, rel=1e-14)
    assert sc.n_c == pytest.approx(nc, rel=1e-14)
    assert sc.n_v == pytest.approx(nv, rel=1e-14)
    assert sc.e_g == pytest.approx(s7.band_gap(reg, 300.0), rel=1e-14)
    assert 1e-25 < sc.delta < 1e-23


# ----------------------------------------------------------------------------- unit conversions
def test_unit_conversions_explicit(sc):
    vt = K_B * 300.0 / Q
    assert sc.vt == pytest.approx(vt, rel=1e-14)
    eps = 8.8541878128e-14 * sc.eps_r  # F/cm
    ld = np.sqrt(eps * vt / (Q * 1e16))
    assert sc.l_d == pytest.approx(ld, rel=1e-12)
    assert 3.5e-6 < sc.l_d < 4.0e-6
    e_phys = 2.5e6  # V/cm -> scaled V_T/L_D
    assert float(sc.field_s(e_phys)) == pytest.approx(e_phys * ld / vt, rel=1e-12)
    assert float(sc.field_phys(sc.field_s(e_phys))) == pytest.approx(e_phys, rel=1e-14)
    assert float(sc.alpha_s(1.0e4)) == pytest.approx(1.0e4 * ld, rel=1e-14)
    assert float(sc.mu_s(520.0)) == pytest.approx(520.0 / sc.mu0, rel=1e-14)
    assert float(sc.dens_s(2.5e17)) == pytest.approx(25.0, rel=1e-14)
    assert float(sc.sheet_s(1e12)) == pytest.approx(1e12 / (1e16 * ld), rel=1e-12)
    assert float(sc.volt_s(1.0)) == pytest.approx(1.0 / vt, rel=1e-14)
    # velocity: mu_s E_s / v_s = mu E / v for any physical triple
    mu, e, v = 450.0, 1.0e5, 2.2e7
    lhs = float(sc.mu_s(mu)) * float(sc.field_s(e)) / float(sc.vsat_s(v))
    assert lhs == pytest.approx(mu * e / v, rel=1e-12)
    # physical drift velocity of one scaled mobility unit in one scaled field unit
    assert sc.mu0 * sc.e_unit == pytest.approx(sc.mu0 * vt / ld, rel=1e-14)


def test_current_conversion_a_per_cm(sc):
    dev = SiCDMOS(sc, DESIGN)
    # J = q n mu E with n = N0, mu = mu0 (scaled 1), E = V_T/L_D (scaled 1) is J0 [A/cm^2];
    # one scaled current unit (per L_D of contact length) is therefore J0 * L_D [A/cm]
    j0 = Q * 1e16 * sc.mu0 * (sc.vt / sc.l_d)
    assert dev.current_a_per_cm(1.0) == pytest.approx(j0 * sc.l_d, rel=1e-12)
    assert dev.current_a_per_cm(1.0) == pytest.approx(Q * 1e16 * sc.mu0 * sc.vt, rel=1e-12)
    assert sc.current_unit_a_per_cm2() == pytest.approx(j0, rel=1e-12)


def test_ionization_coefficients_scaled_vs_physical(sc):
    ni, ha = sc.niwa(), sc.hatakeyama()
    r = sc.reg
    for e_phys in (1.2e6, 2.0e6, 3.0e6):
        es = float(sc.field_s(e_phys))
        an = r.value("oc_a_n") * np.exp(-((r.value("oc_b_n") / e_phys) ** r.value("oc_delta_n")))
        ap = r.value("oc_a_p") * np.exp(-((r.value("oc_b_p") / e_phys) ** r.value("oc_delta_p")))
        assert ni.alpha_n(es)[()] == pytest.approx(an * sc.l_d, rel=1e-12)
        assert ni.alpha_p(es)[()] == pytest.approx(ap * sc.l_d, rel=1e-12)
        hn = r.value("ii_a_n_par_c") * np.exp(-r.value("ii_b_n_par_c") / e_phys)
        assert ha.alpha_n(es)[()] == pytest.approx(hn * sc.l_d, rel=1e-12)
        # identical to the stage-9b Chynoweth with scaled constants
        ch = Chynoweth(ha.an * sc.l_d, ha.bn / sc.e_unit, ha.ap * sc.l_d, ha.bp / sc.e_unit)
        assert ch.alpha_p(es)[()] == pytest.approx(ha.alpha_p(es)[()], rel=1e-12)
        # analytic derivatives vs central differences
        for obj in (ni, ha):
            for f, df in ((obj.alpha_n, obj.dalpha_n), (obj.alpha_p, obj.dalpha_p)):
                h = 1e-5 * es
                fd = (f(es + h) - f(es - h)) / (2 * h)
                assert df(es)[()] == pytest.approx(fd[()], rel=1e-6)
    assert ni.alpha_n(0.0)[()] == 0.0 and ha.alpha_p(0.0)[()] == 0.0  # zero-field floor


def test_caughey_thomas_hall_mobility_identities(sc):
    r = sc.reg
    mx, mn, nref = (r.value(k) for k in ("mu_n_max_par_c", "mu_n_min_par_c", "N_n_ref_par_c"))
    assert sc.mobility_cm2(1e-30, "n", "par") == pytest.approx(mx, rel=1e-14)
    assert sc.mobility_cm2(nref, "n", "par") == pytest.approx(0.5 * (mx + mn), rel=1e-12)
    assert sc.mobility_cm2(1e30, "n", "par") == pytest.approx(mn, rel=1e-6)
    assert sc.mobility_cm2(1e17, "n", "perp") < r.value("mu_n_max_perp_c")


# ----------------------------------------------------------------------------- flat band
def test_flat_band_derivation_identities(sc):
    v_fb, v_code, psi_b = sc.vfb(3e17)
    dphi = sc.reg.value("PhiC_Si_SiO2") - sc.reg.value("PhiC_4HSiC_SiO2")
    assert dphi == pytest.approx(0.45, abs=1e-12)
    # (E_c - E_F)_bulk = (E_c - E_i) - (E_F - E_i), E_i from n_i and N_c, E_F - E_i = V_T psi_b
    ec_ei = sc.vt * np.log(sc.n_c / sc.n_i)
    ec_ef = ec_ei - sc.vt * psi_b
    assert v_fb == pytest.approx(dphi - ec_ef, rel=1e-12)
    # solver convention psi_gate = V_G - V_FB,code with V_FB,code = V_FB/V_T - psi_b
    assert v_code == pytest.approx(v_fb / sc.vt - psi_b, rel=1e-12)
    # independent of the body doping
    assert sc.vfb(1e17)[1] == pytest.approx(v_code, rel=1e-12)
    assert sc.vfb(1e17)[0] != pytest.approx(v_fb, abs=0.01)  # the physical V_FB does depend on N_A
    # the body is p-type, and the flat-band voltage is negative for an n+ poly gate
    assert psi_b < -30.0 and v_fb < 0.0


def test_flat_band_matches_mos_gate_convention(sc):
    # at V_G = V_FB the oxide-side potential equals the neutral body potential psi_b
    v_fb, v_code, psi_b = sc.vfb(3e17)
    assert (v_fb / sc.vt) - v_code == pytest.approx(psi_b, rel=1e-12)


# ----------------------------------------------------------------------------- 12a machinery
def _p_n_junction(sc, nd_cm=1e17, na_cm=1e19, n=801, lo=-3.0, hi=60.0):
    x = np.linspace(lo, hi, n)
    nd = cv_average_piecewise(x, [0.0], [0.0, sc.dens_s(nd_cm)])
    na = cv_average_piecewise(x, [0.0], [sc.dens_s(na_cm), 0.0])
    return x, nd, na


def test_reverse_diode_at_zero_bias_equals_ionized_junction(sc):
    x, nd, na = _p_n_junction(sc)
    lev_d, lev_a = sc.donor_levels(), sc.acceptor_levels(1e19)
    d = IonizedReverseDiode(x, nd, na, sc.delta, lev_d, lev_a)
    psi, rep = d.solve(0.0)
    assert rep.converged
    j = IonizedJunction(x, nd, na, sc.delta, lev_d, lev_a)
    ref, rep2 = j.solve(psi, d.psi_anode_eq, d.psi_cathode_eq)
    assert rep2.converged
    assert np.max(np.abs(psi - ref)) < 1e-9
    # contact potentials are the incomplete-ionization neutral roots
    assert d.psi_cathode_eq == pytest.approx(
        neutral_potential(nd[-1], na[-1], sc.delta, lev_d, lev_a), abs=1e-12
    )


def test_reverse_diode_field_vs_depletion_and_jacobian(sc):
    x, nd, na = _p_n_junction(sc, hi=120.0, n=1201)
    d = IonizedReverseDiode(x, nd, na, sc.delta, sc.donor_levels(), sc.acceptor_levels(1e19))
    v = float(sc.volt_s(100.0))
    psi, rep = d.solve(v)
    assert rep.converged and rep.iterations < 30
    # one-sided abrupt junction: E_max ~ sqrt(2 N_D (V + V_bi - 2))  (all scaled, N_D in N0)
    vbi = d.psi_cathode_eq - d.psi_anode_eq
    e_dep = np.sqrt(2.0 * float(sc.dens_s(1e17)) * (v + vbi - 2.0))
    assert d.max_field(psi) == pytest.approx(e_dep, rel=0.03)
    # analytic Jacobian vs central FD on the residual
    jb = d.jacobian_banded(psi, v)
    n = len(psi)
    col = n // 2
    h = 1e-6
    e = np.zeros(n)
    e[col] = h
    fd = (d.residual(psi + e, v) - d.residual(psi - e, v)) / (2 * h)
    an = np.zeros(n)
    an[col] = jb[1, col]
    an[col - 1] = jb[0, col]
    an[col + 1] = jb[2, col]
    sel = slice(col - 2, col + 3)
    assert np.max(np.abs(fd[sel] - an[sel])) < 1e-6 * np.max(np.abs(an[sel]))


def test_large_bias_jump_uses_depletion_guess(sc):
    x, nd, na = _p_n_junction(sc, hi=400.0, n=801)
    d = IonizedReverseDiode(x, nd, na, sc.delta, sc.donor_levels(), sc.acceptor_levels(1e19))
    p0, r0 = d.solve(float(sc.volt_s(5.0)))
    p1, r1 = d.solve(float(sc.volt_s(300.0)), psi0=p0)  # huge jump: guess is rebuilt
    assert r0.converged and r1.converged


def test_ionization_integral_identity_and_bv_forms(sc):
    ha = sc.hatakeyama()
    nd_cm = 1e17
    nd_s = float(sc.dens_s(nd_cm))
    e_m, w_d, v_bd = triangular_breakdown(nd_s, ha)
    d, info = build_pn_diode(sc, nd_cm, 1.25 * w_d, level=0, rule={"w_dep": w_d})
    bn = avalanche_bv(d, ha, "n", v_step=v_bd / 3, v_max=6 * v_bd)
    bp = avalanche_bv(d, ha, "p", v_step=v_bd / 3, v_max=6 * v_bd)
    assert abs(bn.bv - bp.bv) / bn.bv < 1e-8  # stage-9b identity I_n = 1 <=> I_p = 1
    i_n, i_p, _ = ionization_integrals(d.x, d.element_field(bn.psi), ha)
    assert abs(i_n - i_p) / i_n < 1e-8
    # design estimate within 5% of the solved BV (triangular field vs full electrostatics)
    assert sc.volt_phys(bn.bv) == pytest.approx(sc.volt_phys(v_bd), rel=0.05)
    assert info["n_nodes"] == len(d.x)


def test_uniform_field_closed_form_on_a_column(sc):
    dev = SiCDMOS(sc, DESIGN)
    ch = sc.niwa()
    x_col, nodes = dev.column_nodes(dev.prm.x_max)
    depth = -dev.dev.mesh.node_xy[nodes, 1]
    e0 = float(sc.field_s(2.5e6))
    psi = np.zeros(dev.dev.mesh.n_nodes)
    psi[nodes] = e0 * depth  # uniform field -dpsi/dx = -e0 (magnitude e0)
    i_n, i_p, a_w = dev.column_integrals(psi, nodes, ch)
    an, ap = ch.alpha_n(e0)[()], ch.alpha_p(e0)[()]
    w = depth[-1] - depth[0]
    assert i_n == pytest.approx(uniform_field_integral(an, ap, w), rel=1e-12)
    assert a_w == pytest.approx((an - ap) * w, rel=1e-12)
    assert x_col == dev.dev.mesh.x[-1] and depth[0] == 0.0  # anode at index 0


# ----------------------------------------------------------------------------- 12b machinery
@pytest.fixture(scope="module")
def device(sc):
    return SiCDMOS(sc, DESIGN)


def test_device_construction_uses_registry_values(sc, device):
    p = device.prm
    assert p.eps_ratio == pytest.approx(3.9 / 10.2, rel=1e-14)
    assert p.delta == sc.delta
    assert p.a_body == pytest.approx(30.0) and p.n0_drift == pytest.approx(5.0)
    assert device.dev.vfb != device.vfb_code  # the MOSFET uses the sourced code V_FB
    assert device.itf.dit > 0 and device.itf.q_f == 0.0
    assert device.itf.e2 > device.itf.e1 > 0.0
    # interface window lies at E_c - [0.2, 0] eV (scaled V_T from midgap)
    assert (device.itf.e2 - device.itf.e1) * sc.vt == pytest.approx(0.2, rel=1e-12)
    # per-edge mobility: perp c on x-edges, par c on y-edges, between mu_min and mu_max
    mesh = device.dev.mesh
    assert 0.0 < device.mu_edge.min() and device.mu_edge.max() <= 1260.0 / sc.mu0 + 1e-12
    pure = (np.abs(device.nd - 5.0) < 1e-9) & (device.na == 0.0)
    drift_y = (mesh.edge_dir == 1) & pure[mesh.ea] & pure[mesh.eb]
    assert drift_y.any()
    assert device.mu_edge[drift_y].mean() == pytest.approx(
        float(sc.mu_s(sc.mobility_cm2(5e16, "n", "par"))), rel=1e-6
    )


def test_equilibrium_converges_and_offstate_tangent_predictor(sc, device):
    eq = device.equilibrium()
    assert eq.report.converged
    mod = device.model
    s0 = mod.solve_bc(eq.X, device.offstate_bc(0.0), max_step=20.0)
    assert s0.report.converged
    # zero-current electrostatics at 0 V equals the full equilibrium solution
    assert np.max(np.abs(s0.psi - eq.psi)) < 1e-9
    v1 = float(sc.volt_s(5.0))
    pred = tangent_predictor(mod, s0.X, device.offstate_bc, 0.0, v1)
    s1 = mod.solve_bc(pred, device.offstate_bc(v1), max_step=20.0)
    assert s1.report.converged
    # predictor is a good guess: closer than the unshifted state
    assert np.max(np.abs(pred - s1.X)) < 0.75 * np.max(np.abs(s0.X - s1.X))
    # zero-current: contact currents (scaled) are tiny
    cur = device.currents(s1)
    assert abs(cur["I_drain"]) < 1e-12


def test_on_state_current_balance_and_gate_physics(sc, device):
    from widegap.fv2d.mos import continue_bias

    eq = device.equilibrium()
    mod = device.model
    out = continue_bias(
        mod, eq.X, lambda p: device.bc(p[0], 0.0, p[1]), [0.0, 0.0],
        [[0.0, float(sc.volt_s(0.1))], [float(sc.volt_s(12.0)), float(sc.volt_s(0.1))]],
        max_step=20.0,
    )  # fmt: skip
    sol = out[-1][1]
    assert sol.report.converged
    c = device.currents(sol)
    assert abs(c["I_drain"]) > 1e-4
    assert abs(c["I_source"] + c["I_drain"] + c["I_body"]) <= 1e-8 * abs(c["I_drain"])
    assert c["I_drain"] * c["I_source"] < 0.0  # current leaves through the other terminal


def test_jacobian_with_all_models_and_sourced_parameters(sc, device):
    rng = np.random.default_rng(20261021)
    ys = np.r_[np.linspace(-1.0, 0.0, 4), 0.5]
    m = TensorMesh(np.linspace(0, 1, 3), ys)
    reg = Regions.from_function(m, lambda xc, yc: yc < 0)
    semi = reg.semi_node
    lg = float(np.log(sc.delta))
    worst = 0.0
    for _ in range(2):
        nd = np.where(semi, rng.uniform(1.0, 100.0, m.n_nodes), 0.0)
        na = np.where(semi, rng.uniform(1.0, 100.0, m.n_nodes), 0.0)
        ntot = (nd + na) * sc.n0
        en = 0.5 * (ntot[m.ea] + ntot[m.eb])
        mu_e = sc.mu_s(
            np.where(
                m.edge_dir == 0,
                sc.mobility_cm2(en, "n", "perp"),
                sc.mobility_cm2(en, "n", "par"),
            )
        )
        dd = AnisoMobilityDD(
            reg, nd - na, DDParams(delta=sc.delta, recombination=False), lg,
            ionization=Ionization(nd, na, device.donors, device.acceptors),
            mobility=device.mob, interface=device.itf, mu_edge=mu_e,
        )  # fmt: skip
        X = np.zeros(3 * m.n_nodes)
        X[0::3] = np.where(semi, rng.uniform(-abs(lg), abs(lg), m.n_nodes), 0.0)
        X[1::3] = np.where(semi, rng.uniform(-1.5, 1.5, m.n_nodes), 0.0)
        X[2::3] = np.where(semi, rng.uniform(-1.5, 1.5, m.n_nodes), 0.0)
        bc = Dirichlet()
        rows, _ = dd._all_bc(bc)
        ja = dd.jacobian_bc(X, rows).toarray()
        jf = dd.fd_jacobian(X, bc, h=1e-6).toarray()
        worst = max(worst, float((np.abs(jf - ja).max(axis=1) / np.abs(ja).max(axis=1)).max()))
    assert worst <= 1e-6
