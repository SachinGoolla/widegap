"""Manufactured-solution specifications for the 1D Poisson benchmark.

The benchmark is defined in ``docs/benchmarks/poisson_1d.md``. This module is
the only code shared by the reference solver and the PINN. It holds the problem
data (permittivity, source, exact solution, Dirichlet values) and nothing
about how either method discretizes or trains.

Every function takes an array ``x`` and an array namespace ``xp`` (``numpy``
or ``torch``) so that the same formula is evaluated by both backends.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Callable

import numpy as np

ArrayFn = Callable[[Any, Any], Any]


@dataclass(frozen=True)
class PoissonProblem1D:
    """-(eps u')' = f on (0, 1), u(0) = g0, u(1) = g1."""

    name: str
    eps: ArrayFn
    f: ArrayFn
    u_exact: ArrayFn
    g0: float
    g1: float
    # Points in (0, 1) where eps is discontinuous; empty for smooth eps.
    eps_breaks: tuple[float, ...] = field(default=())
    gating: bool = True


def _ones_like(x, xp):
    return xp.ones_like(x)


def p1_smooth() -> PoissonProblem1D:
    pi = math.pi
    return PoissonProblem1D(
        name="P1",
        eps=_ones_like,
        f=lambda x, xp: pi**2 * xp.sin(pi * x),
        u_exact=lambda x, xp: xp.sin(pi * x) + x,
        g0=0.0,
        g1=1.0,
    )


def p2_variable_eps() -> PoissonProblem1D:
    pi = math.pi
    return PoissonProblem1D(
        name="P2",
        eps=lambda x, xp: 1.0 + x**2,
        f=lambda x, xp: -2.0 * x * (pi * xp.cos(pi * x) + 1.0)
        + (1.0 + x**2) * pi**2 * xp.sin(pi * x),
        u_exact=lambda x, xp: xp.sin(pi * x) + x,
        g0=0.0,
        g1=1.0,
    )


def p3_interface(kappa: float = 2.5, xs: float = 0.5) -> PoissonProblem1D:
    a1 = 1.0 / (xs + (1.0 - xs) / kappa)
    a2 = a1 / kappa

    def eps(x, xp):
        return xp.where(x < xs, xp.ones_like(x), kappa * xp.ones_like(x))

    def u(x, xp):
        return xp.where(x <= xs, a1 * x, a1 * xs + a2 * (x - xs))

    return PoissonProblem1D(
        name="P3",
        eps=eps,
        f=lambda x, xp: xp.zeros_like(x),
        u_exact=u,
        g0=0.0,
        g1=1.0,
        eps_breaks=(xs,),
    )


def p4_steep_layer(delta: float = 0.02) -> PoissonProblem1D:
    def s(x):
        return (x - 0.5) / delta

    def f(x, xp):
        t = xp.tanh(s(x))
        return 2.0 / delta**2 * (1.0 - t**2) * t

    return PoissonProblem1D(
        name="P4",
        eps=_ones_like,
        f=f,
        u_exact=lambda x, xp: xp.tanh(s(x)),
        g0=math.tanh(-0.5 / delta),
        g1=math.tanh(0.5 / delta),
        gating=False,
    )


def linear_charge_free() -> PoissonProblem1D:
    """u = 0.3 + 0.4 x with eps = 1, f = 0 (used by acceptance check R5)."""
    return PoissonProblem1D(
        name="linear",
        eps=_ones_like,
        f=lambda x, xp: xp.zeros_like(x),
        u_exact=lambda x, xp: 0.3 + 0.4 * x,
        g0=0.3,
        g1=0.7,
    )


CASES: dict[str, Callable[[], PoissonProblem1D]] = {
    "P1": p1_smooth,
    "P2": p2_variable_eps,
    "P3": p3_interface,
    "P4": p4_steep_layer,
}


def evaluate(fn: ArrayFn, x: np.ndarray) -> np.ndarray:
    """Evaluate a problem function on a numpy array."""
    return np.asarray(fn(np.asarray(x, dtype=float), np), dtype=float)
