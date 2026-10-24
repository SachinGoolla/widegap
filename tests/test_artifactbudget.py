"""Solve-free tests of the stage-23 artifact-budget bookkeeping (widegap.artifactbudget)."""

import math

import pytest

from widegap.artifactbudget import (
    Cites,
    artifact_fraction,
    blocking_voltage,
    entry,
    rel_change,
)


def test_artifact_fraction():
    assert artifact_fraction(100.0, 40.0) == pytest.approx(0.6)
    assert artifact_fraction(100.0, 100.0) == 0.0
    assert artifact_fraction(10.0, 90.0) == pytest.approx(-8.0)  # real larger than claimed
    assert math.isnan(artifact_fraction(0.0, 1.0))


def test_blocking_voltage_cases():
    a = blocking_voltage(600.0, 95.0, 225.0)  # punch-through before avalanche
    assert a == {"blocking_V": 95.0, "lower_bound_only": False, "limited_by": "punch-through"}
    b = blocking_voltage(90.0, 400.0, 400.0)  # avalanche first
    assert b["blocking_V"] == 90.0 and b["limited_by"] == "avalanche"
    c = blocking_voltage(600.0, None, 1000.0)  # DD reached beyond BV without onset
    assert c["blocking_V"] == 600.0 and not c["lower_bound_only"]
    d = blocking_voltage(600.0, None, 230.0)  # DD failed before the onset
    assert d["blocking_V"] == 230.0 and d["lower_bound_only"]


def test_rel_change():
    assert rel_change(100.0, 103.0) == pytest.approx(0.03)
    assert rel_change(None, None) == 0.0
    assert rel_change(None, 5.0) is None


def test_cites_resolve_and_raise(tmp_path):
    fake = {"a.json": {"x": {"y": [1.5, 2.5]}}}
    c = Cites(tmp_path, loader=lambda p: fake[p.name])
    assert c.get("a.json", "x", "y", 1) == 2.5
    r = c.ref("a.json", "x", "y", 0)
    assert r == {"value": 1.5, "file": "a.json", "key": "x/y/0"}
    assert [u["key"] for u in c.used] == ["x/y/1", "x/y/0"]
    with pytest.raises(KeyError, match="missing key path"):
        c.get("a.json", "x", "z")
    with pytest.raises(KeyError):
        c.get("a.json", "x", "y", 5)


def test_entry_row():
    r = entry("AB9", "lab", "gain", 90.0, 54.0, "V", {"claimed": "c"}, note="n", extra_key=1)
    assert r["artifact_fraction"] == pytest.approx(0.4)
    assert r["extra_key"] == 1 and r["sources"] == {"claimed": "c"}
    z = entry("AB0", "lab", "gain", 0.0, 1.0, "V", {})
    assert z["artifact_fraction"] is None


def test_recorded_sources_exist_if_results_present():
    """The runner's recorded-source keys resolve (skipped when the result files are absent)."""
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    if not (root / "results/stage22/results.json").exists():
        pytest.skip("result files not present")
    c = Cites(root)
    assert c.get("results/stage22/results.json", "p22b", "runs", "D", "nit") == 37
    assert c.get("results/stage17b/results.json", "R4", "bv_true_min_h", "17b") > 500.0
