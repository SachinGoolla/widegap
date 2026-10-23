"""Guarded sparse solve.

SuperLU (scipy ``spsolve``) must never see NaN/Inf: with non-finite input it can call
BLAS with illegal arguments (OpenBLAS prints "On entry to DTRSV parameter number 6 had an
illegal value") and has been observed to corrupt memory in multiprocessing runs.
``safe_spsolve`` returns an all-NaN vector instead, so that every Newton loop's
"non-finite step" test rejects the step / reports non-convergence. For finite inputs it
calls ``spsolve`` with exactly the same arguments (bit-identical results).
"""

from __future__ import annotations

import numpy as np
from scipy.sparse.linalg import spsolve

#: number of solves skipped because the matrix data or right-hand side was non-finite
nonfinite_calls = 0


def safe_spsolve(A, b):
    """``spsolve(A, b)``, or an all-NaN vector of shape ``(A.shape[1],)`` (``b``'s trailing
    shape for a matrix rhs) if ``A.data`` or ``b`` contains NaN/Inf."""
    global nonfinite_calls
    data = A.data if hasattr(A, "data") else np.asarray(A)
    b_arr = np.asarray(b.toarray() if hasattr(b, "toarray") else b)
    if not (np.isfinite(data).all() and np.isfinite(b_arr).all()):
        nonfinite_calls += 1
        shape = (A.shape[1],) + tuple(b_arr.shape[1:])
        return np.full(shape, np.nan)
    return spsolve(A, b)
