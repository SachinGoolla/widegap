"""Stage 15f: the Ishikawa v_sat option of SiC4H / SiCDMOS (default behaviour unchanged)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from widegap.fv2d.physics import ct_factor
from widegap.materials.registry import ProvenanceError, Registry, check_entry
from widegap.sic import REGISTRY_PATH, DMOSDesign, SiC4H, SiCDMOS, TracedRegistry

PENDING = Path(REGISTRY_PATH).with_name("sic4h_pending_15f.json")
DESIGN = DMOSDesign(
    x_jfet=0.8, ch_len=0.8, src_w=1.0, bc_w=0.8, pad_w=0.4, d_src=0.3, d_body=0.8, w_drift=3.0,
    t_drain=0.3, t_ox=0.05, n_drift=5e16, n_jfet=1e17, n_body=3e17, n_src=1e19, n_drain=1e19,
    mesh_kw=dict(ny_surf=10, q_surf=0.6),
)  # fmt: skip
ALPHA = {"par": {"2e15": 0.62, "5e15": 0.78, "1e16": 0.87, "3e16": 1.24},
         "perp": {"2e15": 0.58, "5e15": 1.02, "1e16": 0.93, "3e16": 1.06}}  # fmt: skip


def _pending():
    return json.loads(PENDING.read_text())["parameters"]


def _stub_registry():
    """DEV-only merged registry: pending entries given a stub verification (test mechanics only)."""
    ents = json.loads(Path(REGISTRY_PATH).read_text())["parameters"]
    have = {e["name"] for e in ents}
    ents += [{**e, "verified_by": "test-stub", "verified_value_match": True}
             for e in _pending() if e["name"] not in have]  # fmt: skip
    return TracedRegistry(Registry(ents))


def test_default_set_is_khan_cooper_and_reads_only_its_entries():
    sc = SiC4H(TracedRegistry.from_json())
    assert sc.vsat_set == "khan_cooper"
    reg = sc.reg
    reg.used.clear()
    assert sc.saturation() == (2.2e7, 1.2)
    assert set(reg.used) == {"vsat_n_perp_c", "ct_beta_n_perp_c"}


def test_default_device_mobility_unchanged():
    sc = SiC4H(TracedRegistry.from_json())
    dev = SiCDMOS(sc, DESIGN, 0)
    assert dev.mob.beta_n == dev.mob.beta_p == 1.2
    assert dev.mob.vsat_n == dev.mob.vsat_p == float(sc.vsat_s(2.2e7))


def test_unknown_set_raises():
    with pytest.raises(ValueError, match="unknown vsat_set"):
        SiC4H(TracedRegistry.from_json(), vsat_set="bogus").saturation()


def test_ishikawa_without_verified_entries_raises_no_fallback():
    reg = TracedRegistry.from_json()
    if reg.usable("vsat_n_perp_c_ishikawa"):  # after the independent verification is merged
        pytest.skip("Ishikawa entries are registered and verified")
    with pytest.raises((KeyError, ProvenanceError)):
        SiC4H(reg, vsat_set="ishikawa", vsat_beta=1.0).saturation()


def test_pending_entries_follow_schema_and_are_refused_by_check_entry():
    ents = _pending()
    assert [e["name"] for e in ents] == [
        "vsat_n_par_c", "vsat_n_perp_c_ishikawa",
        *[f"ct_beta_n_{o}_c_ishikawa_{n}" for o in ("par", "perp") for n in ALPHA["par"]],
    ]  # fmt: skip
    ref = next(
        e
        for e in json.loads(Path(REGISTRY_PATH).read_text())["parameters"]
        if e["name"] == "mu_n_max_par_c"
    )
    for e in ents:
        for f in ("doi", "doi_agency", "citation", "access", "doi_metadata_match"):
            assert e[f] == ref[f]
        assert e["extracted_by"] == "literature-agent-15f (sonnet)"
        verified = e["verified_by"] is not None
        assert verified or e["verified_value_match"] is False
        assert e["locator"] and e["validity"] and e["param_class"] == "fit-to-data"
        if e["name"].startswith("vsat"):
            assert e["value"] == 1.7e7 and e["unit"] == "cm/s"
        else:
            o, n = e["name"].split("_c_ishikawa_")[0].split("_")[-1], e["name"].rsplit("_", 1)[1]
            assert e["value"] == ALPHA[o][n] and e["unit"] == "1"
        if not verified:
            with pytest.raises(ProvenanceError):  # the provenance check is not weakened
                check_entry(e)
        elif e["verified_value_match"] is True:
            assert e["verified_by"] != e["extracted_by"]
            check_entry(e)


def test_khan_cooper_entries_untouched():
    main = {e["name"]: e for e in json.loads(Path(REGISTRY_PATH).read_text())["parameters"]}
    assert main["vsat_n_perp_c"]["value"] == 2.2e7 and main["ct_beta_n_perp_c"]["value"] == 1.2
    assert all(e["name"] not in main for e in _pending()) or all(
        main[e["name"]]["verified_by"] for e in _pending()
    )  # either still pending or merged with a verifier


def test_ishikawa_requires_explicit_beta():
    with pytest.raises(ValueError, match="explicit vsat_beta"):
        SiC4H(_stub_registry(), vsat_set="ishikawa").saturation()


@pytest.mark.parametrize("beta", [1.0, 1.2])
def test_ishikawa_set_values_with_stub_registry(beta):
    reg = _stub_registry()
    sc = SiC4H(reg, vsat_set="ishikawa", vsat_beta=beta)
    reg.used.clear()
    assert sc.saturation() == (1.7e7, beta)
    assert not any(n.startswith("ct_beta") for n in reg.used)  # beta is not a registry value
    dev = SiCDMOS(sc, DESIGN, 0)
    assert dev.mob.beta_n == dev.mob.beta_p == beta
    assert dev.mob.vsat_n == dev.mob.vsat_p == float(sc.vsat_s(1.7e7))
    # the default set on the same registry is still Khan & Cooper
    assert SiC4H(reg).saturation() == (2.2e7, 1.2)


def test_ishikawa_requires_isotropic_vsat():
    ents = json.loads(Path(REGISTRY_PATH).read_text())["parameters"]
    ents += [{**e, "verified_by": "test-stub", "verified_value_match": True} for e in _pending()]
    for e in ents:
        if e["name"] == "vsat_n_par_c":
            e["value"] = 1.9e7
    with pytest.raises(ValueError, match="not isotropic"):
        SiC4H(TracedRegistry(Registry(ents)), vsat_set="ishikawa", vsat_beta=1.0).saturation()


def test_stage9a_caughey_thomas_rejects_beta_below_one():
    """Tripwire for the 15f finding: the registered Ishikawa alpha at N_D = 1e16 (0.93, perp c)
    is < 1 and the stage-9a Caughey-Thomas implementation rejects beta < 1, so the pre-registered
    V2 solves cannot run until that is amended (reported, not worked around)."""
    assert ALPHA["perp"]["1e16"] < 1.0
    with pytest.raises(ValueError, match="beta"):
        ct_factor(1.0, 1.0, 1.0, 1.0, 0.93)
