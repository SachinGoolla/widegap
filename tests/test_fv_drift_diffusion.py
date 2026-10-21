"""Unit tests for the stage-3 drift-diffusion solver (properties, not the gate)."""

import numpy as np
import pytest

from widegap.fv.drift_diffusion import DDParams, DriftDiffusion1D, bernoulli, bernoulli_prime
from widegap.fv.equilibrium import cv_average_piecewise, neutral_potential
from widegap.problems.drift_diffusion import D1_PARAMS, d1_fields


def test_bernoulli_identity_and_derivative():
    t = np.r_[-30.0, -1.0, -1e-6, 0.0, 1e-6, 1e-3, 2.0, 30.0]
    # B(-t) = B(t) + t  (exact identity)
    assert np.allclose(bernoulli(-t), bernoulli(t) + t, rtol=1e-12, atol=1e-14)
    e = 1e-6
    fd = (bernoulli(t + e) - bernoulli(t - e)) / (2 * e)
    assert np.allclose(bernoulli_prime(t), fd, rtol=1e-6, atol=1e-9)


def _d1(n_el):
    x = np.linspace(0.0, 1.0, n_el + 1)
    psi, n, p, s_psi, s_n, s_p = d1_fields(x)
    dd = DriftDiffusion1D(
        x, np.zeros_like(x), DDParams(**D1_PARAMS), {"psi": s_psi, "n": s_n, "p": s_p}
    )
    bc = [(psi[0], n[0], p[0]), (psi[-1], n[-1], p[-1])]
    X0 = np.empty(3 * len(x))
    X0[0::3], X0[1::3], X0[2::3] = psi[0] + (psi[-1] - psi[0]) * x, 1.0, 0.7
    return dd, bc, X0, (psi, n, p)


def test_jacobian_matches_finite_differences():
    dd, bc, X0, _ = _d1(40)
    rng = np.random.default_rng(1)
    J = dd.jacobian(X0)
    for _ in range(3):
        v = rng.standard_normal(len(X0))
        v[[0, 1, 2, -3, -2, -1]] = 0.0
        e = 1e-7
        fd = (dd.residual(X0 + e * v, bc) - dd.residual(X0 - e * v, bc)) / (2 * e)
        assert np.max(np.abs(J @ v - fd)) <= 1e-6 * np.max(np.abs(J @ v))


@pytest.mark.parametrize("tol_mode", ["terms", "jacobian"])
def test_mms_second_order(tol_mode):
    errs = []
    for n_el in (32, 64, 128):
        dd, bc, X0, (psi, n, p) = _d1(n_el)
        sol = dd.solve(X0, bc, tol_mode=tol_mode)
        assert sol.report.converged
        errs.append(max(np.max(np.abs(sol.psi - psi)), np.max(np.abs(sol.n - n))))
    orders = np.log2(np.array(errs[:-1]) / np.array(errs[1:]))
    assert np.all((orders > 1.8) & (orders < 2.2)), orders


def test_equilibrium_diode_has_zero_current():
    delta = 1e-3
    x = np.linspace(-20.0, 20.0, 257)
    c = cv_average_piecewise(x, [0.0], [-1.0, 1.0])
    pl, pr = float(neutral_potential(-1.0, delta)), float(neutral_potential(1.0, delta))
    dd = DriftDiffusion1D(x, c, DDParams(delta=delta, recombination=False))
    bc = [
        (pl, delta * np.exp(pl), delta * np.exp(-pl)),
        (pr, delta * np.exp(pr), delta * np.exp(-pr)),
    ]
    X0 = np.empty(3 * len(x))
    psi0 = np.where(x < 0, pl, pr)
    X0[0::3], X0[1::3], X0[2::3] = psi0, delta * np.exp(psi0), delta * np.exp(-psi0)
    sol = dd.solve(X0, bc, tol_mode="jacobian")
    assert sol.report.converged
    assert np.max(np.abs(sol.jn)) < 1e-12 and np.max(np.abs(sol.jp)) < 1e-12
