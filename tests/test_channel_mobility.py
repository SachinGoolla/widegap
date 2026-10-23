"""Stage 16: Noguchi inversion-layer mobility (formulas, clamps, edge selection, Picard driver)."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from widegap.channel_mobility import (
    EEFF_MIN_MV_CM,
    NOGUCHI_NAMES,
    NS_MIN_CM2,
    ChannelMap,
    NoguchiParams,
    depletion_sheet,
    noguchi_mobility,
    picard_solve,
)
from widegap.materials.registry import ProvenanceError, Registry, check_entry
from widegap.sic import REGISTRY_PATH, DMOSDesign, SiC4H, SiCDMOS, TracedRegistry

PENDING = Path(REGISTRY_PATH).with_name("sic4h_pending_16.json")
DESIGN = DMOSDesign(
    x_jfet=0.8, ch_len=0.8, src_w=1.0, bc_w=0.8, pad_w=0.4, d_src=0.3, d_body=0.8, w_drift=3.0,
    t_drain=0.3, t_ox=0.05, n_drift=5e16, n_jfet=1e17, n_body=3e17, n_src=1e19, n_drain=1e19,
    mesh_kw=dict(ny_surf=10, q_surf=0.6),
)  # fmt: skip


def _entries():
    """Registry entries; the Noguchi ones come from sic4h.json once verified, else the pending
    file with a TEST-ONLY stub verification (check_entry refuses the pending entries)."""
    ents = json.loads(Path(REGISTRY_PATH).read_text())["parameters"]
    have = {e["name"] for e in ents}
    pend = json.loads(PENDING.read_text())["parameters"] if PENDING.exists() else []
    ents += [
        {**e, "verified_by": "test-stub", "verified_value_match": True}
        for e in pend
        if e["name"] not in have
    ]
    return ents


@pytest.fixture(scope="module")
def stack():
    reg = TracedRegistry(Registry(_entries()))
    sc = SiC4H(reg)
    return sc, NoguchiParams.from_registry(reg)


def test_pending_entries_are_refused_and_complete():
    pend = {e["name"]: e for e in json.loads(PENDING.read_text())["parameters"]}
    assert set(NOGUCHI_NAMES) == set(pend)
    main = {e["name"] for e in json.loads(Path(REGISTRY_PATH).read_text())["parameters"]}
    for n, e in pend.items():
        assert e["extracted_by"] == "literature-agent-16 (sonnet)"
        if n in main:  # already merged after independent verification
            continue
        assert e["verified_by"] is None
        with pytest.raises(ProvenanceError):
            check_entry(e)
        # the two inferences are recorded in every entry
        assert "UNITS INFERENCE" in e["validity"] and "LOG10 INFERENCE" in e["validity"]


def test_registered_constants(stack):
    _, p = stack
    assert (p.ph_coeff, p.ph_exp, p.sr_coeff, p.sr_exp) == (66.5, -0.39, 76.7, -2.0)
    assert (p.A, p.B, p.C, p.D) == (59.0, -0.939, 0.428, 0.334)
    assert (p.ns0, p.ndep0) == (1.5e12, 1.0e12)
    assert p.eta == pytest.approx(1.0 / 3.0)


def test_gamma_table_log10_reproduces_spec(stack):
    """Spec section 1 pre-registration check: log10 gives 0.229/0.333/0.447/0.578."""
    sc, p = stack
    na = np.array([1e16, 3e16, 1e17, 4e17])
    r = noguchi_mobility(sc, p, na, np.full(4, 1.5e12))
    assert np.allclose(r["gamma"], [0.229, 0.333, 0.447, 0.578], atol=6e-4)
    ln = p.C * np.log(depletion_sheet(sc, na) / p.ndep0) + p.D  # the rejected reading
    assert np.allclose(ln, [0.092, 0.331, 0.593, 0.895], atol=6e-4)
    # N_dep itself (N20 Eq. 3) at the registry eps_r, n_i
    assert np.allclose(
        depletion_sheet(sc, na), [5.675e11, 9.927e11, 1.8317e12, 3.7073e12], rtol=2e-3
    )


def test_phonon_anchor_fig6_units(stack):
    """N17 Fig. 6: about 163 cm^2/Vs at E_eff = 0.1 MV/cm (V/cm would give ~0.75)."""
    _, p = stack
    assert p.ph_coeff * 0.1**p.ph_exp == pytest.approx(163.2, abs=0.1)
    assert p.ph_coeff * (0.1e6) ** p.ph_exp == pytest.approx(0.75, abs=0.02)


def test_matthiessen_and_limit(stack):
    sc, p = stack
    na, ns = np.array([1e17, 3e17]), np.array([3e12, 1e12])
    r = noguchi_mobility(sc, p, na, ns)
    inv = 1 / r["mu_ph"] + 1 / r["mu_sr"] + 1 / r["mu_c"]
    assert np.allclose(r["mu"], 1 / inv, rtol=1e-14)
    assert np.all(r["mu"] < np.minimum(np.minimum(r["mu_ph"], r["mu_sr"]), r["mu_c"]))
    # limit: Coulomb scattering switched off -> phonon + roughness only
    big = NoguchiParams(**{**p.__dict__, "A": 1e30})
    rb = noguchi_mobility(sc, big, na, ns)
    assert np.allclose(rb["mu"], 1 / (1 / rb["mu_ph"] + 1 / rb["mu_sr"]), rtol=1e-12)
    # E_eff = (q/eps)(N_dep + N_S/3) in MV/cm
    eps = 8.8541878128e-12 * sc.eps_r * 1e-2
    want = 1.602176634e-19 / eps * (depletion_sheet(sc, na) + ns / 3.0) * 1e-6
    assert np.allclose(r["e_eff"], want, rtol=1e-9)


def test_clamps(stack):
    sc, p = stack
    lo = noguchi_mobility(sc, p, np.array([1e17, 1e17]), np.array([1e9, NS_MIN_CM2]))
    assert lo["mu"][0] == lo["mu"][1]
    assert lo["clamp_ns"].tolist() == [True, False]
    # tiny N_A: E_eff clamped at 1e-2 MV/cm; N_S above 1e13 and gamma < 0 are flagged, not clamped
    r = noguchi_mobility(sc, p, np.array([1e12, 1e17, 5e14]), np.array([1e11, 5e13, 1e13]))
    assert r["e_eff"][0] == EEFF_MIN_MV_CM and bool(r["clamp_eeff"][0])
    assert r["ns_extrap"].tolist() == [False, True, False]
    assert r["gamma_neg"].tolist() == [True, False, True]
    assert r["na_out_of_range"].tolist() == [True, False, True]
    off = noguchi_mobility(sc, p, np.array([1e17]), np.array([5e10]), clamp=False)
    assert off["ns_used"][0] == 5e10


@pytest.fixture(scope="module")
def device(stack):
    sc, p = stack
    dev = SiCDMOS(sc, DESIGN, 0)
    mu0 = dev.model.mu_edge.copy()
    cmap = ChannelMap(sc, dev.model, dev.nd, dev.na, 20.0e-3, gate_x_max=dev.prm.gate_x1)
    return sc, p, dev, cmap, mu0


def test_edge_selection_and_default_path_untouched(device):
    sc, _, dev, cmap, mu0 = device
    mesh = dev.dev.mesh
    assert np.array_equal(dev.model.mu_edge, mu0)  # constructing the map changes nothing
    e = cmap.edges
    assert len(e) > 0 and np.all(mesh.edge_dir[e] == 0)
    xy = mesh.node_xy
    a, b = dev.model.a[e], dev.model.b[e]
    d_ch = cmap.d_ch_s
    assert np.all(xy[a, 1] >= -d_ch * (1 + 1e-9)) and np.all(xy[a, 1] <= 0.0)
    assert np.all(cmap.col_p[cmap.node_col[a]]) and np.all(cmap.col_p[cmap.node_col[b]])
    # n-type surface columns (source / JFET) are excluded
    assert not cmap.col_p[mesh.node_xy[:, 0].searchsorted(0.0)]
    # replacing mu_ch by the bulk value leaves the edge array untouched
    bulk_cols = np.full(len(cmap.cols), 1.0)
    new = cmap.edge_mobility(bulk_cols)
    assert np.array_equal(np.delete(new, e), np.delete(mu0, e))


def test_bulk_mode_code_path_is_identical_and_channel_lowers_current(device):
    sc, p, dev, cmap, mu0 = device
    vg, vd = float(sc.volt_s(10.0)), float(sc.volt_s(0.1))
    eq = dev.equilibrium()
    from widegap.fv2d.mos import continue_bias

    out = continue_bias(
        dev.model, eq.X, lambda q: dev.bc(q[0], 0.0, q[1]), [0.0, 0.0], [[0.0, vd], [vg, vd]],
        max_step=20.0,
    )  # fmt: skip
    sol = out[-1][1]
    assert sol.report.converged
    i_bulk = abs(dev.model.contact_current(sol, dev.dev.nodes["drain"]))
    rb = picard_solve(cmap, p, sol.X, dev.bc(vg, 0.0, vd), mode="bulk")
    assert rb["converged"] and np.array_equal(dev.model.mu_edge, mu0)
    i_path = abs(dev.model.contact_current(rb["sol"], dev.dev.nodes["drain"]))
    assert i_path == pytest.approx(i_bulk, rel=1e-12)
    rc = picard_solve(cmap, p, sol.X, dev.bc(vg, 0.0, vd))
    assert rc["converged"] and rc["final_rel"] < 1e-10 and rc["iterations"] <= 60
    i_ch = abs(dev.model.contact_current(rc["sol"], dev.dev.nodes["drain"]))
    assert i_ch < i_bulk  # inversion-layer mobility is far below the bulk Hall value
    dev.model.mu_edge = mu0.copy()  # leave the shared fixture clean
