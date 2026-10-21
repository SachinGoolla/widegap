"""Stage-7 infrastructure tests on SYNTHETIC parameters (no SiC values here)."""

import numpy as np
import pytest
from scipy.integrate import quad

from widegap.materials.ionization import (
    Level,
    charge,
    charge_antiderivative,
    dcharge_dpsi,
    neutral_potential,
    single_donor_electrons,
)
from widegap.materials.registry import ProvenanceError, Registry

GOOD = {
    "name": "x",
    "value": 1.0,
    "unit": "eV",
    "param_class": "measured",
    "doi": "10.0/x",
    "crossref_title_match": True,
    "access": "full-text-read",
    "locator": "Table 1",
    "extracted_by": "agent-A",
    "verified_by": "agent-V",
    "verified_value_match": True,
}


@pytest.mark.parametrize(
    "patch",
    [
        {"access": "abstract-only"},
        {"access": "paywalled"},
        {"crossref_title_match": False},
        {"verified_by": "agent-A"},
        {"verified_value_match": False},
        {"locator": ""},
        {"param_class": "guess"},
    ],
)
def test_registry_refuses_unprovenanced(patch):
    reg = Registry([{**GOOD, **patch}])
    with pytest.raises(ProvenanceError):
        reg.get("x")
    assert not reg.usable("x")


def test_registry_accepts_complete_entry_and_reports_missing():
    reg = Registry([GOOD])
    assert reg.value("x") == 1.0
    assert reg.missing(["x", "y"]) == ["y"]


LV_D = (Level(0.5, -3.0), Level(0.5, -1.0))
LV_A = (Level(1.0, 2.0),)


def test_charge_derivative_matches_fd():
    psi = np.linspace(-30, 30, 61)
    e = 1e-6
    fd = (
        charge(psi + e, 1.0, 0.3, 1e-10, LV_D, LV_A) - charge(psi - e, 1.0, 0.3, 1e-10, LV_D, LV_A)
    ) / (2 * e)
    assert np.allclose(dcharge_dpsi(psi, 1.0, 0.3, 1e-10, LV_D, LV_A), fd, rtol=1e-6, atol=1e-9)


def test_antiderivative_matches_quadrature():
    a, b = -5.0, 7.0
    num, _ = quad(lambda p: charge(p, 1.0, 0.3, 1e-4, LV_D, LV_A), a, b, epsabs=1e-13, epsrel=1e-13)
    ex = charge_antiderivative(b, 1.0, 0.3, 1e-4, LV_D, LV_A) - charge_antiderivative(
        a, 1.0, 0.3, 1e-4, LV_D, LV_A
    )
    assert abs(num - ex) <= 1e-10 * max(1.0, abs(ex))


def test_neutrality_matches_single_donor_closed_form():
    # one donor level, n-type, delta tiny so holes are negligible
    delta = 1e-20
    for nd in (1e-3, 1.0, 1e3):
        for lg in (-55.0, -50.0, -45.0):
            psi = neutral_potential(nd, 0.0, delta, (Level(1.0, lg),), ())
            n = delta * np.exp(psi)
            # n1 such that N_D^+ = N_D / (1 + n / n1):  n1 = delta / e^{lg}
            n1 = delta * np.exp(-lg)
            assert abs(n / single_donor_electrons(nd, n1) - 1.0) <= 1e-10
