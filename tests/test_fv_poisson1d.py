"""Acceptance checks R1-R5b for the reference solver (docs/benchmarks/poisson_1d.md)."""

import numpy as np
import pytest

from wbgtcad.fv.poisson1d import solve, uniform_mesh
from wbgtcad.problems.poisson1d import (
    evaluate,
    linear_charge_free,
    p1_smooth,
    p2_variable_eps,
    p3_interface,
)
from wbgtcad.verification import max_error, observed_orders

LEVELS = [32, 64, 128, 256, 512]


def _errors(problem):
    errs = []
    for n in LEVELS:
        sol = solve(problem, uniform_mesh(n))
        errs.append(max_error(sol.u, evaluate(problem.u_exact, sol.x)))
    return errs


@pytest.mark.parametrize("make", [p1_smooth, p2_variable_eps])
def test_R1_second_order_convergence(make):
    orders = observed_orders(_errors(make()))
    assert all(1.9 <= p <= 2.1 for p in orders[-3:]), orders


def test_R2_p1_accuracy():
    assert _errors(p1_smooth())[-1] <= 1e-5


def test_R3_p2_accuracy():
    assert _errors(p2_variable_eps())[-1] <= 5e-5


@pytest.mark.parametrize("n", [64, 63])
def test_R4_dielectric_interface_exact(n):
    prob = p3_interface()
    sol = solve(prob, uniform_mesh(n))
    assert max_error(sol.u, evaluate(prob.u_exact, sol.x)) <= 1e-12


def test_R4_arithmetic_mean_would_fail():
    """Control: the interface case really discriminates the averaging rule."""
    prob = p3_interface()
    x = uniform_mesh(63)
    sol = solve(prob, x)
    # Replace the straddling element's eps by the arithmetic mean and re-solve.
    mid = np.argmax((x[:-1] < 0.5) & (x[1:] > 0.5))
    eps_a = sol.eps_elem.copy()
    eps_a[mid] = 0.5 * (1.0 + 2.5)
    assert not np.isclose(eps_a[mid], sol.eps_elem[mid])
    c = eps_a / np.diff(x)
    # Series conductance across [0, 1] determines the left slope exactly.
    left_slope_flux = 1.0 / np.sum(1.0 / c)
    exact_flux = 1.0 / (0.5 + 0.5 / 2.5)
    assert abs(left_slope_flux - exact_flux) > 1e-3


def test_R5_linear_solution_exact():
    prob = linear_charge_free()
    for n in (7, 64):
        sol = solve(prob, uniform_mesh(n))
        assert max_error(sol.u, evaluate(prob.u_exact, sol.x)) <= 1e-12


def test_R5_linear_solution_exact_nonuniform():
    prob = linear_charge_free()
    x = np.sort(np.r_[0.0, np.random.default_rng(0).uniform(0, 1, 40), 1.0])
    sol = solve(prob, x)
    assert max_error(sol.u, evaluate(prob.u_exact, sol.x)) <= 1e-12


def test_R5b_discrete_conservation():
    prob = p1_smooth()
    x = uniform_mesh(512)
    sol = solve(prob, x)
    h = np.diff(x)
    vol = np.zeros_like(x)
    vol[:-1] += 0.5 * h
    vol[1:] += 0.5 * h
    fv = evaluate(prob.f, x)[1:-1] * vol[1:-1]
    balance = sol.flux[-1] - sol.flux[0]
    assert abs(balance - fv.sum()) <= 1e-10 * np.abs(fv).sum()


def test_symmetry():
    """P1 minus its linear part is symmetric about x = 1/2."""
    sol = solve(p1_smooth(), uniform_mesh(128))
    v = sol.u - sol.x
    assert np.allclose(v, v[::-1], atol=1e-13)


def test_rejects_bad_mesh():
    with pytest.raises(ValueError):
        solve(p1_smooth(), np.array([0.0, 0.5, 0.4, 1.0]))
