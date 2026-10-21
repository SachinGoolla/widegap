"""Unit tests for stage-5 adjoint sensitivities (small, fast cases)."""

import numpy as np

from widegap import adjoint as A
from widegap.fv import equilibrium as eq
from widegap.fv.equilibrium import cv_average_piecewise, neutral_potential


def test_banded_transpose_matches_dense():
    rng = np.random.default_rng(0)
    ab = rng.standard_normal((3, 7))
    ab[0, 0] = ab[2, -1] = 0.0
    dense = np.diag(ab[1]) + np.diag(ab[0, 1:], 1) + np.diag(ab[2, :-1], -1)
    at = A.banded_transpose(ab)
    dense_t = np.diag(at[1]) + np.diag(at[0, 1:], 1) + np.diag(at[2, :-1], -1)
    assert np.array_equal(dense_t, dense.T)


def test_equilibrium_adjoint_matches_directional_fd():
    d = 1e-4
    x = np.linspace(-15.0, 15.0, 121)
    c = cv_average_piecewise(x, [0.0], [-1.0, 1.0])
    pl, pr = float(neutral_potential(-1.0, d)), float(neutral_potential(1.0, d))
    psi = A.polish_equilibrium(eq.solve(x, c, d, pl, pr).psi, x, c, d)
    _, g = A.equilibrium_integral_psi(psi, x, d)
    v = np.sin(np.pi * (x - x[0]) / (x[-1] - x[0]))
    v[0] = v[-1] = 0.0
    vol = eq.control_volumes(x)

    def J(cc):
        p2 = A.polish_equilibrium(eq.solve(x, cc, d, pl, pr, psi0=psi).psi, x, cc, d)
        return float(np.sum(vol * p2))

    e = 1e-4
    fd = (J(c + e * v) - J(c - e * v)) / (2 * e)
    assert abs(g @ v - fd) <= 1e-6 * abs(fd)
    assert g[0] == 0.0 and g[-1] == 0.0  # contacts are not design variables
