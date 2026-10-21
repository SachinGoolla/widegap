"""Unit tests for stage-4 reverse-bias electrostatics (properties, not the gate)."""

import numpy as np

from widegap.fv.equilibrium import cv_average_piecewise
from widegap.fv.equilibrium import residual as residual_stage2
from widegap.fv.reverse_bias import Diode1D, breakdown_voltage, residual, sweep


def _b1(h=0.05):
    x = np.linspace(-5.0, 60.0, int(round(65.0 / h)) + 1)
    return Diode1D(x, cv_average_piecewise(x, [0.0], [-100.0, 1.0]), 1e-3)


def test_zero_bias_operator_is_bitwise_stage2():
    d = _b1()
    psi = np.random.default_rng(0).standard_normal(len(d.x))
    assert np.array_equal(
        residual(psi, d.x, d.c_cv, d.delta),
        residual_stage2(psi, d.x, d.c_cv, d.delta, np.zeros_like(psi)),
    )


def test_reverse_bias_raises_barrier_and_converges():
    d = _b1()
    sw = sweep(d, range(0, 31, 1))
    assert all(r.converged for r in sw.reports)
    drops = [p[-1] - p[0] for p in sw.psi]
    assert np.allclose(np.diff(drops), 1.0, atol=1e-9)  # contact BCs: barrier grows by V_R
    assert np.all(np.diff(sw.e_max) > 0)  # peak field grows monotonically


def test_breakdown_extraction_hits_b_equals_one():
    d = _b1()
    sw = sweep(d, range(0, 121))
    bv, psi, _ = breakdown_voltage(d, sw, 15.0)
    assert abs(d.max_field(psi) / 15.0 - 1.0) < 1e-9
    assert sw.v[0] < bv < sw.v[-1]
