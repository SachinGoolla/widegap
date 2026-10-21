"""Manufactured solution D1 for the stage-3 drift-diffusion benchmark.

Sources are derived symbolically (sympy) from the exact fields so that the
residual definitions in docs/benchmarks/stage3_drift_diffusion.md hold
exactly for the continuous problem.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np


@lru_cache(maxsize=1)
def _d1_lambdas(delta: float, mu_n: float, mu_p: float, tau_n: float, tau_p: float):
    import sympy as sp

    x = sp.symbols("x")
    psi = 2 * x + sp.Rational(3, 10) * sp.sin(2 * sp.pi * x)
    n = sp.exp(sp.sin(sp.pi * x) / 2)
    p = 1 / (1 + x)
    r = (n * p - delta**2) / (tau_p * (n + delta) + tau_n * (p + delta))
    jn = mu_n * (sp.diff(n, x) - n * sp.diff(psi, x))
    jp = -mu_p * (sp.diff(p, x) + p * sp.diff(psi, x))
    s_psi = -sp.diff(psi, x, 2) - (p - n)
    s_n = sp.diff(jn, x) - r
    s_p = sp.diff(jp, x) + r
    return tuple(sp.lambdify(x, e, "numpy") for e in (psi, n, p, s_psi, s_n, s_p))


D1_PARAMS = dict(delta=0.1, mu_n=1.0, mu_p=1.0, tau_n=1.0, tau_p=1.0)


def d1_fields(x: np.ndarray):
    """Exact (psi, n, p) and sources (s_psi, s_n, s_p) at x for case D1."""
    f = _d1_lambdas(**D1_PARAMS)
    x = np.asarray(x, dtype=float)
    return tuple(np.broadcast_to(fi(x), x.shape).astype(float) for fi in f)
