"""Vertex-centred finite-volume (box-method) solver for 1D Poisson.

Discretization (docs/benchmarks/poisson_1d.md, section 3):

* nodes x_0 = 0 < ... < x_N = 1, elements e_i = [x_i, x_{i+1}];
* element permittivity = length-weighted harmonic mean of eps over the
  element, h_e / int_e (1/eps) dx, which makes the flux exact for piecewise
  constant eps even when an interface lies inside an element;
* source lumped at nodes: f(x_i) times the control-volume length;
* Dirichlet values imposed at both end nodes; tridiagonal system solved by
  banded LU.

This module deliberately imports nothing from ``widegap.pinn``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import solve_banded

from widegap.problems.poisson1d import PoissonProblem1D, evaluate

_GAUSS_ORDER = 8


@dataclass(frozen=True)
class FVSolution:
    x: np.ndarray
    u: np.ndarray
    eps_elem: np.ndarray
    # Flux -eps u' on each element (positive in +x direction).
    flux: np.ndarray


def uniform_mesh(n_elements: int) -> np.ndarray:
    return np.linspace(0.0, 1.0, n_elements + 1)


def element_harmonic_eps(problem: PoissonProblem1D, x: np.ndarray) -> np.ndarray:
    """h_e / int_e eps^{-1} dx per element, splitting at eps discontinuities."""
    gx, gw = np.polynomial.legendre.leggauss(_GAUSS_ORDER)
    eps_e = np.empty(len(x) - 1)
    for i, (a, b) in enumerate(zip(x[:-1], x[1:], strict=True)):
        cuts = [a, *(c for c in problem.eps_breaks if a < c < b), b]
        inv_int = 0.0
        for lo, hi in zip(cuts[:-1], cuts[1:], strict=True):
            mid, half = 0.5 * (lo + hi), 0.5 * (hi - lo)
            q = mid + half * gx
            inv_int += half * np.sum(gw / evaluate(problem.eps, q))
        eps_e[i] = (b - a) / inv_int
    return eps_e


def solve(problem: PoissonProblem1D, x: np.ndarray) -> FVSolution:
    x = np.asarray(x, dtype=float)
    if x.ndim != 1 or len(x) < 3 or np.any(np.diff(x) <= 0):
        raise ValueError("x must be a strictly increasing 1D mesh with >= 3 nodes")
    h = np.diff(x)
    eps_e = element_harmonic_eps(problem, x)
    c = eps_e / h  # element conductances

    n = len(x)
    # Banded storage: row 0 super-diagonal, row 1 diagonal, row 2 sub-diagonal.
    ab = np.zeros((3, n))
    rhs = np.zeros(n)

    vol = np.zeros(n)
    vol[:-1] += 0.5 * h
    vol[1:] += 0.5 * h
    fx = evaluate(problem.f, x)

    # Interior rows: sum of outward fluxes = f_i * vol_i.
    ab[1, 1:-1] = c[:-1] + c[1:]
    ab[0, 2:] = -c[1:]  # coefficient of u_{i+1} in row i
    ab[2, :-2] = -c[:-1]  # coefficient of u_{i-1} in row i
    rhs[1:-1] = fx[1:-1] * vol[1:-1]

    # Dirichlet rows.
    ab[1, 0] = 1.0
    ab[0, 1] = 0.0
    rhs[0] = problem.g0
    ab[1, -1] = 1.0
    ab[2, -2] = 0.0
    rhs[-1] = problem.g1

    u = solve_banded((1, 1), ab, rhs)
    flux = -c * np.diff(u)
    return FVSolution(x=x, u=u, eps_elem=eps_e, flux=flux)
