"""The reference solver and the PINN must not share discretization code."""

from pathlib import Path

SRC = Path(__file__).resolve().parents[1] / "src" / "widegap"


def _imports(path: Path) -> str:
    return "\n".join(
        line
        for line in path.read_text().splitlines()
        if line.lstrip().startswith(("import", "from"))
    )


def test_fv_does_not_import_pinn():
    for f in (SRC / "fv").glob("*.py"):
        assert "pinn" not in _imports(f), f


def test_pinn_does_not_import_fv():
    for f in (SRC / "pinn").glob("*.py"):
        assert "fv" not in _imports(f), f


def test_benchmark_summary_is_json_serializable():
    """Regression: numpy bool_ criteria crashed summary.json writing."""
    import importlib.util
    import json

    spec = importlib.util.spec_from_file_location(
        "bench_run", Path(__file__).resolve().parents[1] / "benchmarks" / "poisson1d" / "run.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    fv = mod.run_fv()
    pinn = [{"case": "P1", "e_l2_rel_exact": 1e-6} for _ in mod.SEEDS]
    json.dumps(mod.summarize(fv, pinn))
