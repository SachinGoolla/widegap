"""Stage 23: bookkeeping of the artifact budget (claimed vs real, blocking voltage, citations).

Spec: docs/benchmarks/stage23_artifact_budget.md (pre-registered; AUTHORITATIVE).  Pure functions
only (no solves): the runner ``benchmarks/stage23/run.py`` feeds them recorded and recomputed values.
"""

from __future__ import annotations

import json
import math
from pathlib import Path


def artifact_fraction(claimed: float, real: float) -> float:
    """``1 - real / claimed`` (spec metric); NaN if ``claimed`` is zero."""
    if claimed == 0:
        return float("nan")
    return float(1.0 - real / claimed)


def blocking_voltage(bv: float, v_onset: float | None, v_reached: float) -> dict:
    """``min(BV, V_pt)`` with ``V_pt`` the 1D-DD onset voltage.

    ``v_onset`` None means the DD march ended (first non-converged target) before the onset: then
    ``V_pt >= v_reached`` and the blocking voltage is only bounded below by ``min(BV, v_reached)``
    ("lower bound only"), unless ``v_reached >= BV``, where it is exactly ``BV``."""
    if v_onset is not None:
        return {"blocking_V": float(min(bv, v_onset)), "lower_bound_only": False,
                "limited_by": "punch-through" if v_onset < bv else "avalanche"}  # fmt: skip
    if v_reached >= bv:
        return {"blocking_V": float(bv), "lower_bound_only": False, "limited_by": "avalanche"}
    return {"blocking_V": float(v_reached), "lower_bound_only": True, "limited_by": "undetermined"}


def rel_change(a: float | None, b: float | None) -> float | None:
    """|b / a - 1| (None if either is missing); both None -> 0 (equal ceilings are 'equal')."""
    if a is None and b is None:
        return 0.0
    if a is None or b is None:
        return None
    return float(abs(b / a - 1.0))


class Cites:
    """Fetch recorded values from result files, remembering (file, key path) for the table."""

    def __init__(self, root: Path, loader=None):
        self.root = Path(root)
        self.used: list[dict] = []
        self._cache: dict = {}
        self._loader = loader or (lambda p: json.loads(Path(p).read_text()))

    def get(self, rel_file: str, *keys):
        """Value at ``keys`` of ``rel_file`` (relative to ``root``); raises KeyError if missing."""
        if rel_file not in self._cache:
            self._cache[rel_file] = self._loader(self.root / rel_file)
        v = self._cache[rel_file]
        for k in keys:
            try:
                v = v[k]
            except (KeyError, IndexError, TypeError) as exc:
                raise KeyError(f"{rel_file}: missing key path {list(keys)} at {k!r}") from exc
        self.used.append({"file": rel_file, "key": "/".join(str(k) for k in keys)})
        return v

    def ref(self, rel_file: str, *keys) -> dict:
        """Like :meth:`get`, but returns ``{"value", "file", "key"}``."""
        return {"value": self.get(rel_file, *keys), "file": rel_file,
                "key": "/".join(str(k) for k in keys)}  # fmt: skip


def entry(
    eid: str,
    label: str,
    quantity: str,
    claimed: float,
    real: float,
    unit: str,
    sources: dict,
    note: str = "",
    **extra,
) -> dict:
    """One row of the artifact-budget table."""
    row = {"id": eid, "label": label, "quantity": quantity, "unit": unit,
           "claimed": float(claimed), "real": float(real),
           "artifact_fraction": artifact_fraction(claimed, real),
           "sources": sources, "note": note}  # fmt: skip
    row.update(extra)
    if not math.isfinite(row["artifact_fraction"]):
        row["artifact_fraction"] = None
    return row
