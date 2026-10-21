"""Nonlinear equilibrium Poisson equation, Boltzmann statistics (stage 2).

Scaled form (docs/equations.md, docs/benchmarks/stage2_equilibrium_pn.md):

    -psi'' = delta * (exp(-psi) - exp(psi)) + C(x) + s(x)

with psi in units of V_T = k_B T / q, x in units of the extrinsic Debye length
L_D = sqrt(eps V_T / (q N0)), and densities in units of N0; delta = n_i / N0.
``s`` is an optional manufactured source (zero for physical problems).

Discretization: vertex-centred finite volume (as in ``widegap.fv.poisson1d``)
with the nonlinear charge lumped at nodes. Newton's method with the analytic
Jacobian; steps whose max-norm exceeds ``max_step`` are scaled down to it
(documented limiter; every limited step is counted). Imports nothing from
``widegap.pinn``.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.linalg import solve_banded


@dataclass
class NewtonReport:
    converged: bool
    iterations: int
    limited_steps: int
    residual_history: list[float] = field(default_factory=list)
    limited_flags: list[bool] = field(default_factory=list)


@dataclass
class EquilibriumSolution:
    x: np.ndarray
    psi: np.ndarray
    n: np.ndarray  # scaled electron density n / N0
    p: np.ndarray  # scaled hole density p / N0
    report: NewtonReport


def neutral_potential(c: np.ndarray, delta: float) -> np.ndarray:
    """Charge-neutral potential: delta (e^-psi - e^psi) + C = 0."""
    return np.arcsinh(np.asarray(c, dtype=float) / (2.0 * delta))


def control_volumes(x: np.ndarray) -> np.ndarray:
    h = np.diff(x)
    vol = np.zeros_like(x)
    vol[:-1] += 0.5 * h
    vol[1:] += 0.5 * h
    return vol


def cv_average_piecewise(x: np.ndarray, breaks: list[float], values: list[float]) -> np.ndarray:
    """Exact control-volume average of a piecewise-constant function.

    ``values[k]`` holds on the k-th interval delimited by sorted ``breaks``
    (len(values) == len(breaks) + 1).
    """
    lo = np.r_[x[0], 0.5 * (x[:-1] + x[1:])]
    hi = np.r_[0.5 * (x[:-1] + x[1:]), x[-1]]
    edges = [-np.inf, *breaks, np.inf]
    acc = np.zeros_like(x)
    for k, v in enumerate(values):
        a, b = edges[k], edges[k + 1]
        acc += v * np.clip(np.minimum(hi, b) - np.maximum(lo, a), 0.0, None)
    return acc / (hi - lo)


def residual(psi, x, c_cv, delta, s_cv):
    """Finite-volume balance F_i (flux form, multiplied by control volume)."""
    h = np.diff(x)
    vol = control_volumes(x)
    flux = -np.diff(psi) / h  # -psi' on elements
    f = np.zeros_like(psi)
    f[1:-1] = (flux[1:] - flux[:-1]) - vol[1:-1] * (
        delta * (np.exp(-psi[1:-1]) - np.exp(psi[1:-1])) + c_cv[1:-1] + s_cv[1:-1]
    )
    return f


def solve(
    x: np.ndarray,
    c_cv: np.ndarray,
    delta: float,
    psi_left: float,
    psi_right: float,
    s_cv: np.ndarray | None = None,
    psi0: np.ndarray | None = None,
    tol: float = 1e-12,
    max_iter: int = 40,
    max_step: float = 2.0,
) -> EquilibriumSolution:
    x = np.asarray(x, dtype=float)
    if x.ndim != 1 or len(x) < 3 or np.any(np.diff(x) <= 0):
        raise ValueError("x must be a strictly increasing 1D mesh with >= 3 nodes")
    c_cv = np.asarray(c_cv, dtype=float)
    s_cv = np.zeros_like(x) if s_cv is None else np.asarray(s_cv, dtype=float)
    psi = neutral_potential(c_cv, delta) if psi0 is None else np.array(psi0, dtype=float)
    psi[0], psi[-1] = psi_left, psi_right

    h = np.diff(x)
    g = 1.0 / h
    vol = control_volumes(x)
    rep = NewtonReport(converged=False, iterations=0, limited_steps=0)

    for _ in range(max_iter):
        f = residual(psi, x, c_cv, delta, s_cv)
        rnorm = float(np.max(np.abs(f)))
        rep.residual_history.append(rnorm)
        if rnorm <= tol:
            rep.converged = True
            break
        ab = np.zeros((3, len(x)))
        ab[1, 1:-1] = g[:-1] + g[1:] + vol[1:-1] * delta * (np.exp(-psi[1:-1]) + np.exp(psi[1:-1]))
        ab[0, 2:] = -g[1:]
        ab[2, :-2] = -g[:-1]
        ab[1, 0] = ab[1, -1] = 1.0
        dpsi = solve_banded((1, 1), ab, -f)
        step = float(np.max(np.abs(dpsi)))
        limited = step > max_step
        if limited:
            dpsi *= max_step / step
            rep.limited_steps += 1
        rep.limited_flags.append(limited)
        psi = psi + dpsi
        rep.iterations += 1

    return EquilibriumSolution(
        x=x, psi=psi, n=delta * np.exp(psi), p=delta * np.exp(-psi), report=rep
    )
