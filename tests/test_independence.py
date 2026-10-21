"""The reference solver and the PINN must not share discretization code."""

from pathlib import Path

SRC = Path(__file__).resolve().parents[1] / "src" / "wbgtcad"


def _imports(path: Path) -> str:
    return "\n".join(
        line for line in path.read_text().splitlines() if line.lstrip().startswith(("import", "from"))
    )


def test_fv_does_not_import_pinn():
    for f in (SRC / "fv").glob("*.py"):
        assert "pinn" not in _imports(f), f


def test_pinn_does_not_import_fv():
    for f in (SRC / "pinn").glob("*.py"):
        assert "fv" not in _imports(f), f
