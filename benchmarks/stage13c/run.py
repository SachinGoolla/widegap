"""Stage-13c benchmark (docs/benchmarks/stage12c_offstate_fix.md section 4): the stage-13b re-run
with the corrected off-state model; U1-U7 (gating), R1-R2.

EVERYTHING except the off-state quasi-Fermi model is the stage-13b run, executed unchanged: this
script imports benchmarks/stage13b/run.py, replaces only its ``CornerAdjoint`` by
``widegap.offstate.CorrectedCornerAdjoint`` (phi_p = 0, phi_n = V_D w, w the discrete harmonic
extension; see widegap.offstate and the stage-12c spec), and renames the output (default
results/stage13c) and the figures (stage13c_*.png; repo figures/ only for the default output,
else <out>/figures).  Same rounded device (r = 0.1 um), K = 10 cuts, tau = 5 V, 49 dense lines,
theta bounds, Pareto levels, SLSQP and thresholds U1-U7 as stage 13b.

After the run, <out>/comparison.json and summary.json["stage13c_vs_13b"] record the comparison
with results/stage13b (dense true-minimum BV on h and h/2, corner peak field, limiting
locations).
"""

# ruff: noqa: E402
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import importlib.util
import json
import sys
from pathlib import Path

from widegap.offstate import CorrectedCornerAdjoint

REPO = Path(__file__).resolve().parents[2]
DEFAULT_OUT = REPO / "results" / "stage13c"
B13 = REPO / "results" / "stage13b"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


M = _load("stage13b_run", REPO / "benchmarks/stage13b/run.py")  # fixed name: fork-pool pickling
M.CornerAdjoint = CorrectedCornerAdjoint  # the only physics change (harmonic w)
M.DEFAULT_OUT = DEFAULT_OUT


def _relabel(fn, pos_png, pos_label):
    """Wrap a 13b plotting function: stage13b_ -> stage13c_ in the file name, 13b -> 13c in the label."""

    def wrapper(*a, **k):
        a = list(a)
        p = Path(a[pos_png])
        a[pos_png] = p.with_name(p.name.replace("stage13b_", "stage13c_"))
        a[pos_label] = str(a[pos_label]).replace("stage 13b", "stage 13c")
        return fn(*a, **k)

    wrapper.__name__ = fn.__name__
    return wrapper


M.plot_per_cut = _relabel(M.plot_per_cut, 3, 4)
M.plot_profile = _relabel(M.plot_profile, 3, 4)
M.S13.plot_maps = _relabel(M.S13.plot_maps, 3, 4)


def compare(b13_summary, b13_results, c_summary, c_results):
    """13c vs 13b: dense true-minimum BV (h, h/2), corner peak fields, limiting locations."""
    v13, v13c = b13_summary["values"], c_summary["values"]

    def pair(f):
        return {"stage13b": f(v13, b13_results), "stage13c": f(v13c, c_results)}

    return {
        "criteria": {"stage13b": b13_summary["criteria"], "stage13c": c_summary["criteria"]},
        "U5_dense_true_min_BV_V": {
            "initial_h_h2": pair(lambda v, r: v["U5"]["min_bv_initial_h_h2"]),
            "final_h_h2": pair(lambda v, r: v["U5"]["min_bv_final_h_h2"]),
            "improvement_h": pair(lambda v, r: v["U5"]["improvement_h"]),
            "improvement_h2": pair(lambda v, r: v["U5"]["improvement_h2"]),
            "x_argmin_initial_um": pair(lambda v, r: v["U5"]["x_argmin_initial_um"]),
            "x_argmin_final_um": pair(lambda v, r: v["U5"]["x_argmin_final_um"]),
        },
        "U1_corner_peak_field": {
            "E_peak_MV_cm_h_h2_h4": pair(lambda v, r: v["U1"]["E_peak_MV_cm_h_h2_h4"]),
            "v_volt": pair(lambda v, r: v["U1"]["v_volt"]),
            "order": pair(lambda v, r: v["U1"]["order"]),
            "rel_dist_richardson": pair(lambda v, r: v["U1"]["rel_dist"]),
            "locations_um": pair(lambda v, r: v["U1"]["locations_um"]),
        },
        "U6_global_peak_at_initial_BV": pair(lambda v, r: v["U6"]),
        "objective": {
            "bv_obj_initial_V": pair(lambda v, r: v["bv_obj_initial_V"]),
            "bv_obj_final_V": pair(lambda v, r: v["bv_obj_final_V"]),
            "limiting_cut_um_initial": pair(lambda v, r: r["initial"]["limiting_cut_x_um"]),
            "limiting_cut_um_final": pair(lambda v, r: r["final"]["limiting_cut_x_um"]),
            "R_initial_final": pair(lambda v, r: [v["R_initial"], v["R_final"]]),
        },
        "stage13_optimum_on_rounded_device": pair(lambda v, r: v["stage13_comparison"]),
        "blind_spot": pair(lambda v, r: v["blind_spot"]),
        "theta_final": pair(lambda v, r: r["final"]["theta"]),
    }


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--maxiter", type=int, default=M.MAX_ITER)
    ap.add_argument("--workers", type=int, default=min(10, os.cpu_count() or 1))
    args = ap.parse_args(argv)
    out = args.out
    M.main(["--out", str(out), "--maxiter", str(args.maxiter), "--workers", str(args.workers)])
    res = json.loads((out / "results.json").read_text())
    res["config"]["offstate_model"] = (
        "stage 12c: phi_p = 0, phi_n = V_D w, w = discrete harmonic extension (gw-weighted Laplace "
        "on the semiconductor graph; 0 on n-type nodes of the source box, 1 on other n-type nodes); "
        "classification of the theta = 0 device"
    )
    res["config"]["offstate_node_classification"] = (
        "theta = 0 device; n-type nodes in the source box enlarged by 2.5 r have w = 0, other "
        "n-type nodes w = 1 (stage 13b rule), the remaining nodes carry the harmonic w "
        "(stage 12c) instead of phi_n = 0"
    )
    (out / "results.json").write_text(json.dumps(res, indent=2, default=float) + "\n")
    summ = json.loads((out / "summary.json").read_text())
    if (B13 / "summary.json").exists() and (B13 / "results.json").exists():
        cmp_ = compare(json.loads((B13 / "summary.json").read_text()),
                       json.loads((B13 / "results.json").read_text()), summ, res)  # fmt: skip
        summ["stage13c_vs_13b"] = cmp_
        (out / "comparison.json").write_text(json.dumps(cmp_, indent=2, default=float) + "\n")
        (out / "summary.json").write_text(json.dumps(summ, indent=2, default=float) + "\n")


if __name__ == "__main__":
    main()
