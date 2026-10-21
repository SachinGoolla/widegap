"""Stage-2 benchmark (docs/benchmarks/stage2_equilibrium_pn.md): criteria S1-S6.

Writes results/stage2/{results,summary,manifest}.json. S7 (PINN) is not run:
its configuration must first be fixed by an amendment (see the benchmark doc).
"""

from __future__ import annotations

import argparse
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

from widegap.fv.equilibrium import cv_average_piecewise, neutral_potential, solve
from widegap.manifest import build_manifest
from widegap.verification import observed_orders

REPO = Path(__file__).resolve().parents[2]


def _e1(n: int):
    delta = 1e-3

    def exact(x):
        return 3.0 * np.tanh(4.0 * (x - 0.5))

    def source(x):
        t = np.tanh(4.0 * (x - 0.5))
        return 96.0 * (1.0 - t**2) * t - delta * (np.exp(-exact(x)) - np.exp(exact(x)))

    x = np.linspace(0.0, 1.0, n + 1)
    t0 = time.perf_counter()
    sol = solve(x, np.zeros_like(x), delta, exact(0.0), exact(1.0), s_cv=source(x))
    return sol, float(np.max(np.abs(sol.psi - exact(x)))), time.perf_counter() - t0


def _junction(n: int, a: float, b: float, delta: float = 1e-8, half_width: float = 40.0):
    """Abrupt junction at x=0: N_A = a (x<0), N_D = b (x>0), scaled by N0."""
    x = np.linspace(-half_width, half_width, n + 1)
    c = cv_average_piecewise(x, [0.0], [-a, b])
    pl, pr = neutral_potential(-a, delta), neutral_potential(b, delta)
    t0 = time.perf_counter()
    sol = solve(x, c, delta, pl, pr)
    elapsed = time.perf_counter() - t0
    vol = np.diff(x)
    rho = sol.p - sol.n + c  # scaled space charge at nodes
    bulk = np.abs(x) >= 20.0
    neutral = np.max(np.abs(rho[bulk]) / np.abs(c[bulk]))
    trap = np.sum(0.5 * (rho[1:] + rho[:-1]) * vol)
    trap_abs = np.sum(0.5 * (np.abs(c[1:]) + np.abs(c[:-1])) * vol)
    vbi = float(pr - pl)
    e_num = float(np.max(np.abs(np.diff(sol.psi) / vol)))
    e_dep = float(np.sqrt(2.0 * a * b * (vbi - 2.0) / (a + b)))
    return sol, {
        "Vbi_thermal": vbi,
        "bulk_neutrality_max_rel": float(neutral),
        "total_charge_rel": float(abs(trap) / trap_abs),
        "peak_field_numeric": e_num,
        "peak_field_depletion": e_dep,
        "peak_field_rel_diff": abs(e_num - e_dep) / e_dep,
        "solve_seconds": elapsed,
    }


def _newton(rep) -> dict:
    """S2 quantities. Ratio for unlimited iteration j: r_{j+1} / r_j^2."""
    r, flags = rep.residual_history, rep.limited_flags
    ratios = [r[j + 1] / r[j] ** 2 for j in range(len(flags)) if not flags[j]]
    return {
        "converged": rep.converged,
        "iterations": rep.iterations,
        "limited_steps": rep.limited_steps,
        "final_residual": r[-1],
        "residual_history": r,
        "limited_flags": flags,
        "unlimited_ratios": ratios,
        "final3_unlimited_ratios": ratios[-3:],
        "S2_evaluable": len(ratios) >= 3,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=REPO / "results" / "stage2")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "manifest.json").write_text(
        json.dumps(build_manifest(REPO, sys.argv, []), indent=2) + "\n"
    )

    levels = [64, 128, 256, 512, 1024]
    e1 = [_e1(n) for n in levels]
    errs = [e for _, e, _ in e1]
    res = {
        "E1": {
            "levels": levels,
            "e_inf": errs,
            "orders": observed_orders(errs),
            "solve_seconds": [t for _, _, t in e1],
            "newton_N256": _newton(e1[levels.index(256)][0].report),
        }
    }
    sol2, m2 = _junction(1024, 1.0, 1.0)
    res["E2"] = {**m2, "N": 1024, "newton": _newton(sol2.report)}
    sol3, m3 = _junction(1024, 100.0, 1.0)
    res["E3"] = {**m3, "N": 1024, "N_A/N_D": 100.0, "newton": _newton(sol3.report)}
    res["peak_rss_mib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0

    def s2(nw):
        return (
            nw["converged"]
            and nw["S2_evaluable"]
            and all(q <= 1e3 for q in nw["final3_unlimited_ratios"])
        )

    crit = {
        "S1": all(1.9 <= p <= 2.1 for p in res["E1"]["orders"][-3:]),
        "S2_E1": s2(res["E1"]["newton_N256"]),
        "S2_E2": s2(res["E2"]["newton"]),
        "S3": res["E2"]["bulk_neutrality_max_rel"] <= 1e-6,
        "S4": res["E2"]["total_charge_rel"] <= 1e-8,
        "S5": res["E2"]["peak_field_rel_diff"] <= 0.10,
    }
    crit = {k: bool(v) for k, v in crit.items()}
    summary = {
        "criteria": crit,
        "gate_passed": all(crit.values()),
        "reported_S6_E3": {k: res["E3"][k] for k in m3},
        "S7_pinn": "not run (configuration amendment required first)",
    }
    (args.out / "results.json").write_text(json.dumps(res, indent=2) + "\n")
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
