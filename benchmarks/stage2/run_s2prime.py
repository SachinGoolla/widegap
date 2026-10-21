"""Stage-2b benchmark: re-evaluate stage 2 with criterion S2' (docs/benchmarks/stage2b_newton_criterion.md).

Writes results/stage2b/{results,summary,manifest}.json.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

from widegap.fv.equilibrium import (
    banded_matvec,
    control_volumes,
    cv_average_piecewise,
    jacobian_banded,
    neutral_potential,
    residual,
    solve,
)
from widegap.manifest import build_manifest

REPO = Path(__file__).resolve().parents[2]
EPS_MACH = np.finfo(float).eps

_spec = importlib.util.spec_from_file_location("stage2_run", Path(__file__).with_name("run.py"))
stage2 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(stage2)


def _e1_problem(n):
    delta = 1e-3

    def exact(x):
        return 3.0 * np.tanh(4.0 * (x - 0.5))

    def source(x):
        t = np.tanh(4.0 * (x - 0.5))
        return 96.0 * (1.0 - t**2) * t - delta * (np.exp(-exact(x)) - np.exp(exact(x)))

    x = np.linspace(0.0, 1.0, n + 1)
    return dict(
        x=x,
        c_cv=np.zeros_like(x),
        delta=delta,
        psi_left=exact(0.0),
        psi_right=exact(1.0),
        s_cv=source(x),
    )


def _junction_problem(n, a, b, delta=1e-8, half_width=40.0):
    x = np.linspace(-half_width, half_width, n + 1)
    return dict(
        x=x,
        c_cv=cv_average_piecewise(x, [0.0], [-a, b]),
        delta=delta,
        psi_left=float(neutral_potential(-a, delta)),
        psi_right=float(neutral_potential(b, delta)),
        s_cv=np.zeros_like(x),
    )


def jacobian_check(prob, seed=0) -> dict:
    """S2'a: analytic J v vs central differences at initial, midway and final states."""
    full = solve(**prob)
    mid = solve(**prob, max_iter=max(1, full.report.iterations // 2))
    psi0 = neutral_potential(prob["c_cv"], prob["delta"])
    psi0[0], psi0[-1] = prob["psi_left"], prob["psi_right"]
    rng = np.random.default_rng(seed)
    args = (prob["x"], prob["c_cv"], prob["delta"], prob["s_cv"])
    worst, eps = 0.0, 1e-6
    for psi in (psi0, mid.psi, full.psi):
        ab = jacobian_banded(psi, prob["x"], prob["delta"])
        for _ in range(5):
            v = rng.standard_normal(len(psi))
            v[0] = v[-1] = 0.0
            jv = banded_matvec(ab, v)[1:-1]
            fd = (residual(psi + eps * v, *args) - residual(psi - eps * v, *args))[1:-1] / (2 * eps)
            worst = max(worst, float(np.max(np.abs(jv - fd)) / np.max(np.abs(jv))))
    return {"max_rel_err": worst, "states": 3, "directions": 5}


def rate_check(prob) -> dict:
    """S2'c/d: quadratic ratios on unlimited steps whose output residual exceeds the floor."""
    sol = solve(**prob)
    x, psi, d = prob["x"], sol.psi, prob["delta"]
    flux = np.abs(np.diff(psi) / np.diff(x))
    rho = d * (np.exp(-psi) - np.exp(psi)) + prob["c_cv"] + prob["s_cv"]
    terms = np.zeros_like(psi)
    terms[1:-1] = flux[1:] + flux[:-1] + control_volumes(x)[1:-1] * np.abs(rho[1:-1])
    floor = 1e3 * EPS_MACH * float(np.max(terms))
    r, flags = sol.report.residual_history, sol.report.limited_flags
    q = [r[j + 1] / r[j] ** 2 for j in range(len(flags)) if not flags[j] and r[j + 1] > floor]
    return {
        "converged": sol.report.converged,
        "iterations": sol.report.iterations,
        "limited_steps": sol.report.limited_steps,
        "floor": floor,
        "qualifying_ratios": q,
        "pass": sol.report.converged and len(q) >= 2 and all(v <= 1e3 for v in q),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=REPO / "results" / "stage2b")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "manifest.json").write_text(
        json.dumps(build_manifest(REPO, sys.argv, [0]), indent=2) + "\n"
    )

    e1 = _e1_problem(256)
    e2 = _junction_problem(1024, 1.0, 1.0)
    e3 = _junction_problem(1024, 100.0, 1.0)

    # Unchanged criteria S1, S3-S5 recomputed with the stage-2 functions.
    levels = [64, 128, 256, 512, 1024]
    e1_runs = [stage2._e1(n) for n in levels]
    orders = stage2.observed_orders([e for _, e, _ in e1_runs])
    _, m2 = stage2._junction(1024, 1.0, 1.0)

    conv = {
        f"E1_N{n}": sol.report.converged and sol.report.iterations <= 40
        for n, (sol, _, _) in zip(levels, e1_runs, strict=True)
    }
    for name, prob in (("E2", e2), ("E3", e3)):
        s = solve(**prob)
        conv[name] = s.report.converged and s.report.iterations <= 40

    res = {
        "S1_orders": orders,
        "S2a_jacobian": {"E1": jacobian_check(e1), "E2": jacobian_check(e2)},
        "S2b_converged": conv,
        "S2c_rate_E2": rate_check(e2),
        "S2d_rate_E1": rate_check(e1),
        "S2d_rate_E3": rate_check(e3),
        "E2": m2,
    }
    crit = {
        "S1": all(1.9 <= p <= 2.1 for p in orders[-3:]),
        "S2'a": all(v["max_rel_err"] <= 1e-6 for v in res["S2a_jacobian"].values()),
        "S2'b": all(conv.values()),
        "S2'c": res["S2c_rate_E2"]["pass"],
        "S3": m2["bulk_neutrality_max_rel"] <= 1e-6,
        "S4": m2["total_charge_rel"] <= 1e-8,
        "S5": m2["peak_field_rel_diff"] <= 0.10,
    }
    crit = {k: bool(v) for k, v in crit.items()}
    summary = {
        "criteria": crit,
        "gate_passed": all(crit.values()),
        "reported_S2'd": {"E1": res["S2d_rate_E1"]["pass"], "E3": res["S2d_rate_E3"]["pass"]},
    }
    (args.out / "results.json").write_text(json.dumps(res, indent=2) + "\n")
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
