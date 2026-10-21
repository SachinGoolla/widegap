"""Regenerate the README benchmark tables from raw results JSON.

Rewrites the text between ``<!-- results:poisson1d:begin -->`` and
``<!-- results:poisson1d:end -->`` in README.md, so that every number in the
README comes from results/poisson1d/*.json.

Run:  python benchmarks/poisson1d/report.py [--results DIR] [--readme FILE]
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
BEGIN, END = "<!-- results:poisson1d:begin -->", "<!-- results:poisson1d:end -->"


def _e(v: float) -> str:
    return f"{v:.1e}"


def render(res: Path) -> str:
    fv = json.loads((res / "fv.json").read_text())
    man = json.loads((res / "manifest.json").read_text())
    summ = json.loads((res / "summary.json").read_text())
    pinn = sorted(
        (json.loads(p.read_text()) for p in res.glob("pinn_*_seed*.json")),
        key=lambda r: (r["case"], r["seed"]),
    )
    lv = fv["levels"]
    out = [
        "### Finite-volume reference solver (uniform meshes)",
        "",
        "Max-norm nodal error against the exact solution (`results/poisson1d/fv.json`):",
        "",
        "| Case | " + " | ".join(f"N={n}" for n in lv) + " | Observed orders |",
        "|---|" + "---|" * (len(lv) + 1),
    ]
    for c in ("P1", "P2", "P4"):
        d = fv["cases"][c]
        orders = ", ".join(f"{p:.2f}" for p in d["orders_inf"])
        out.append(f"| {c} | " + " | ".join(_e(v) for v in d["e_inf"]) + f" | {orders} |")
    p3 = fv["cases"]["P3"]
    t512 = fv["cases"]["P1"]["solve_seconds"][-1]
    out += [
        "",
        f"- P3 (dielectric interface): max error {_e(p3['64'])} (N=64, interface on a node), "
        f"{_e(p3['63'])} (N=63, interface inside an element).",
        f"- Linear exact solution: max error {_e(fv['cases']['linear'])}.",
        f"- R5b discrete flux balance (P1, N=512): relative imbalance {_e(fv['R5b_rel_imbalance'])}.",
        "- R5c end-element boundary-flux error, N=32..512: "
        + ", ".join(_e(v) for v in fv["R5c_boundary_flux_err"])
        + "; observed orders "
        + ", ".join(f"{p:.2f}" for p in fv["R5c_orders"])
        + ".",
        f"- FV solve time at N=512: {t512 * 1e3:.1f} ms.",
        "",
        "### PINN (3x32 tanh, float64, exact Dirichlet ansatz, Adam then L-BFGS)",
        "",
        "Every seed that was run is listed. Relative L2 errors are evaluated on a fixed 2001-point grid.",
        "",
        "| Case | Seed | rel. L2 vs exact | rel. L2 vs FV (N=512) | final residual MSE | train time (s) | peak RSS (MiB) |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in pinn:
        out.append(
            f"| {r['case']} | {r['seed']} | {_e(r['e_l2_rel_exact'])} | {_e(r['e_l2_rel_vs_fv512'])} "
            f"| {_e(r['final_residual_mse'])} | {r['train_seconds']:.1f} | {r['peak_rss_mib']:.0f} |"
        )
    crit = ", ".join(f"{k} {'pass' if v else 'FAIL'}" for k, v in summ["criteria"].items())
    out += [
        "",
        f"**Gate:** {'passed' if summ['gate_passed'] else 'FAILED'} ({crit}). "
        f"Pre-registered configuration: {summ['pre_registered_configuration']}.",
        "",
        f"Provenance: commit `{man['git_commit'][:7]}` (clean tree: {not man['git_dirty']}), "
        f"{man['cpu_model']}, {man['torch_threads']} torch threads, CPU only, "
        f"torch {man['packages']['torch']}, NumPy {man['packages']['numpy']}, "
        f"SciPy {man['packages']['scipy']}, run at {man['timestamp_utc'][:19]} UTC.",
    ]
    return "\n".join(out)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--results", type=Path, default=REPO / "results" / "poisson1d")
    ap.add_argument("--readme", type=Path, default=REPO / "README.md")
    args = ap.parse_args()
    text = args.readme.read_text()
    if BEGIN not in text or END not in text:
        raise SystemExit(f"markers {BEGIN} / {END} not found in {args.readme}")
    block = f"{BEGIN}\n{render(args.results)}\n{END}"
    args.readme.write_text(
        re.sub(re.escape(BEGIN) + r".*?" + re.escape(END), lambda _: block, text, flags=re.S)
    )


if __name__ == "__main__":
    main()
