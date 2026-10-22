"""Material-parameter registry with enforced provenance (stage 7).

Rules (docs/benchmarks/stage7_sic_materials.md, section 1): a parameter is
usable only if its full text was read, its DOI metadata was checked with its registration agency,
it has a locator, and it was independently re-verified by someone other than
its extractor. Unusable parameters raise ``ProvenanceError``; there is no
fallback value.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

REQUIRED_FIELDS = (
    "name",
    "value",
    "unit",
    "param_class",
    "doi",
    "doi_metadata_match",
    "access",
    "locator",
    "extracted_by",
    "verified_by",
    "verified_value_match",
)
CLASSES = ("measured", "fit-to-data", "first-principles", "convention", "owner-approved-convention")


class ProvenanceError(ValueError):
    """Raised when a parameter does not satisfy the provenance rules."""


@dataclass(frozen=True)
class Parameter:
    name: str
    value: float
    unit: str
    param_class: str
    doi: str
    locator: str
    temperature_K: float | None
    orientation: str | None
    validity: str | None
    record: dict

    @property
    def citation(self) -> str:
        return self.record.get("citation", self.doi)


def check_entry(e: dict) -> None:
    if e.get("param_class") == "owner-approved-convention":
        # Documented exception (stage-7 amendment 2026-10-21): no read source;
        # requires an explicit, dated owner approval recorded in the entry.
        for f in ("name", "value", "unit", "approved_by", "approved_on", "rationale"):
            if f not in e or e[f] in (None, ""):
                raise ProvenanceError(f"{e.get('name', '?')}: owner-approved entry missing {f}")
        return
    missing = [f for f in REQUIRED_FIELDS if f not in e or e[f] in (None, "")]
    if missing:
        raise ProvenanceError(f"{e.get('name', '?')}: missing fields {missing}")
    if e["param_class"] not in CLASSES:
        raise ProvenanceError(f"{e['name']}: unknown class {e['param_class']!r}")
    if e["access"] != "full-text-read":
        raise ProvenanceError(f"{e['name']}: full text not read (access={e['access']!r})")
    if e["doi_metadata_match"] is not True:
        raise ProvenanceError(
            f"{e['name']}: DOI metadata not verified against its registration agency"
        )
    if e["verified_by"] == e["extracted_by"]:
        raise ProvenanceError(f"{e['name']}: not independently re-verified")
    if e["verified_value_match"] is not True:
        raise ProvenanceError(f"{e['name']}: re-verifier did not confirm the value")


class Registry:
    def __init__(self, entries: list[dict]):
        self._raw = {e["name"]: e for e in entries}

    @classmethod
    def from_json(cls, path: str | Path) -> Registry:
        data = json.loads(Path(path).read_text())
        return cls(data["parameters"])

    def names(self) -> list[str]:
        return sorted(self._raw)

    def usable(self, name: str) -> bool:
        try:
            self.get(name)
        except (ProvenanceError, KeyError):
            return False
        return True

    def get(self, name: str) -> Parameter:
        if name not in self._raw:
            raise KeyError(f"parameter {name!r} not in registry")
        e = self._raw[name]
        check_entry(e)
        return Parameter(
            name=e["name"],
            value=float(e["value"]),
            unit=e["unit"],
            param_class=e["param_class"],
            doi=e["doi"],
            locator=e["locator"],
            temperature_K=e.get("temperature_K"),
            orientation=e.get("orientation"),
            validity=e.get("validity"),
            record=e,
        )

    def value(self, name: str) -> float:
        return self.get(name).value

    def missing(self, required: list[str]) -> list[str]:
        return [r for r in required if not self.usable(r)]
