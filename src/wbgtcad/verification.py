"""Error norms and convergence-order estimates used by tests and benchmarks."""

from __future__ import annotations

import numpy as np


def max_error(u_h: np.ndarray, u: np.ndarray) -> float:
    return float(np.max(np.abs(np.asarray(u_h) - np.asarray(u))))


def _trapz_norm(x: np.ndarray, v: np.ndarray) -> float:
    return float(np.sqrt(np.trapezoid(v**2, x)))


def rel_l2_error(x: np.ndarray, u_h: np.ndarray, u: np.ndarray) -> float:
    """||u_h - u||_2 / ||u||_2 with trapezoidal weights on nodes ``x``."""
    x, u_h, u = (np.asarray(a, dtype=float) for a in (x, u_h, u))
    return _trapz_norm(x, u_h - u) / _trapz_norm(x, u)


def observed_orders(errors: list[float], ratio: float = 2.0) -> list[float]:
    """log(e_k / e_{k+1}) / log(ratio) for successive refinements."""
    e = np.asarray(errors, dtype=float)
    return list(np.log(e[:-1] / e[1:]) / np.log(ratio))
