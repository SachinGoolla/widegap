"""Unit tests for the stage-2 equilibrium solver (properties, not the gate itself)."""

import numpy as np
import pytest

from widegap.fv.equilibrium import cv_average_piecewise, neutral_potential, solve


def test_neutral_potential_solves_neutrality():
    delta = 1e-8
    c = np.array([-1.0, -1e-3, 0.0, 2.0])
    psi = neutral_potential(c, delta)
    assert np.allclose(delta * (np.exp(-psi) - np.exp(psi)) + c, 0.0, atol=1e-12)


def test_cv_average_is_exact_for_step():
    x = np.linspace(-1.0, 1.0, 5)  # node at 0 -> its CV is split half/half
    c = cv_average_piecewise(x, [0.0], [-1.0, 1.0])
    assert np.allclose(c, [-1.0, -1.0, 0.0, 1.0, 1.0])


def test_uniform_doping_is_flat_and_converges_immediately():
    x = np.linspace(0.0, 10.0, 51)
    delta = 1e-6
    psi_b = neutral_potential(1.0, delta)
    sol = solve(x, np.ones_like(x), delta, float(psi_b), float(psi_b))
    assert sol.report.converged and sol.report.iterations <= 1
    assert np.allclose(sol.psi, psi_b, atol=1e-12)


def test_mms_second_order():
    delta = 1e-3

    def exact(x):
        return 3.0 * np.tanh(4.0 * (x - 0.5))

    def source(x):
        t = np.tanh(4.0 * (x - 0.5))
        return 96.0 * (1.0 - t**2) * t - delta * (np.exp(-exact(x)) - np.exp(exact(x)))

    errs = []
    for n in (64, 128, 256):
        x = np.linspace(0.0, 1.0, n + 1)
        sol = solve(x, np.zeros_like(x), delta, exact(0.0), exact(1.0), s_cv=source(x))
        assert sol.report.converged
        errs.append(np.max(np.abs(sol.psi - exact(x))))
    orders = np.log2(np.array(errs[:-1]) / np.array(errs[1:]))
    assert np.all((orders > 1.9) & (orders < 2.1)), orders


def test_symmetric_junction_is_antisymmetric():
    delta = 1e-8
    x = np.linspace(-40.0, 40.0, 513)
    c = cv_average_piecewise(x, [0.0], [-1.0, 1.0])
    sol = solve(x, c, delta, neutral_potential(-1.0, delta), neutral_potential(1.0, delta))
    assert sol.report.converged
    assert np.allclose(sol.psi, -sol.psi[::-1], atol=1e-9)


def test_rejects_bad_mesh():
    with pytest.raises(ValueError):
        solve(np.array([0.0, 1.0, 0.5]), np.zeros(3), 1e-3, 0.0, 0.0)
