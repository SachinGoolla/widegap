"""Stage-3 benchmark (docs/benchmarks/stage3_drift_diffusion.md, incl. amendments).

Writes results/stage3/{results,summary,manifest}.json. If results/stage3/devsim.json
exists (benchmarks/stage3/devsim_crosscheck.py), D7 is reported from it.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

from widegap.fv.drift_diffusion import DDParams, DriftDiffusion1D
from widegap.fv.equilibrium import cv_average_piecewise, neutral_potential
from widegap.fv.equilibrium import solve as solve_equilibrium
from widegap.manifest import build_manifest
from widegap.problems.drift_diffusion import D1_PARAMS, d1_fields
from widegap.verification import observed_orders

REPO = Path(__file__).resolve().parents[2]
L, NDIODE, DELTA = 40.0, 1024, 1e-3
VBI = float(2.0 * np.arcsinh(1.0 / (2.0 * DELTA)))


def _pack(psi, n, p):
    X = np.empty(3 * len(psi))
    X[0::3], X[1::3], X[2::3] = psi, n, p
    return X


def _jac_check(dd, X, bc, seed=0):
    rng = np.random.default_rng(seed)
    J = dd.jacobian(X)
    worst = 0.0
    for _ in range(5):
        v = rng.standard_normal(len(X))
        v[[0, 1, 2, -3, -2, -1]] = 0.0
        e = 1e-7
        fd = (dd.residual(X + e * v, bc) - dd.residual(X - e * v, bc)) / (2 * e)
        jv = J @ v
        worst = max(worst, float(np.max(np.abs(jv - fd)) / np.max(np.abs(jv))))
    return worst


def _newton_summary(rep):
    return {
        "converged": rep.converged,
        "iterations": rep.iterations,
        "limited_steps": rep.limited_steps,
        "final_residuals": rep.residual_history[-1],
        "tolerance": rep.tolerance,
        "residual_history": rep.residual_history,
    }


def run_d1():
    levels = [32, 64, 128, 256, 512]
    errs = {"psi": [], "n": [], "p": []}
    newton, jac = {}, {}
    for N in levels:
        x = np.linspace(0.0, 1.0, N + 1)
        psi, n, p, s_psi, s_n, s_p = d1_fields(x)
        dd = DriftDiffusion1D(
            x, np.zeros_like(x), DDParams(**D1_PARAMS), {"psi": s_psi, "n": s_n, "p": s_p}
        )
        bc = [(psi[0], n[0], p[0]), (psi[-1], n[-1], p[-1])]
        X0 = _pack(psi[0] + (psi[-1] - psi[0]) * x, np.ones_like(x), 0.7 * np.ones_like(x))
        t0 = time.perf_counter()
        sol = dd.solve(X0, bc)
        el = time.perf_counter() - t0
        for k, a, b in (("psi", sol.psi, psi), ("n", sol.n, n), ("p", sol.p, p)):
            errs[k].append(float(np.max(np.abs(a - b))))
        newton[N] = {**_newton_summary(sol.report), "seconds": el}
        if N == 128:
            jac["D1_solution"] = _jac_check(dd, _pack(sol.psi, sol.n, sol.p), bc)
            jac["D1_initial"] = _jac_check(dd, X0, bc)
        newton[N]["positive"] = bool(np.all(sol.n > 0) and np.all(sol.p > 0))
    orders = {k: observed_orders(v) for k, v in errs.items()}
    return {"levels": levels, "e_inf": errs, "orders": orders, "newton": newton}, jac


def _diode(params: DDParams, biases_fwd, biases_rev):
    x = np.linspace(-L, L, NDIODE + 1)
    c = cv_average_piecewise(x, [0.0], [-1.0, 1.0])
    pl, pr = float(neutral_potential(-1.0, DELTA)), float(neutral_potential(1.0, DELTA))
    eq = solve_equilibrium(x, c, DELTA, pl, pr)
    dd = DriftDiffusion1D(x, c, params)
    nl, pl_ = DELTA * np.exp(pl), DELTA * np.exp(-pl)
    nr, pr_ = DELTA * np.exp(pr), DELTA * np.exp(-pr)

    def bc(v):
        return [(pl + v, nl, pl_), (pr, nr, pr_)]

    X_eq = _pack(eq.psi, eq.n, eq.p)
    out = {}
    jac = None
    for branch in (biases_fwd, biases_rev):
        X = X_eq.copy()
        for v in branch:
            t0 = time.perf_counter()
            sol = dd.solve(X, bc(v))
            el = time.perf_counter() - t0
            X = _pack(sol.psi, sol.n, sol.p)
            jm = float(np.mean(sol.j))
            out[v] = {
                "j_mean": jm,
                "jn_max_abs": float(np.max(np.abs(sol.jn))),
                "jp_max_abs": float(np.max(np.abs(sol.jp))),
                "j_rel_spread": float(np.max(np.abs(sol.j - jm)) / abs(jm)) if jm else None,
                "positive": bool(np.all(sol.n > 0) and np.all(sol.p > 0)),
                "seconds": el,
                **_newton_summary(sol.report),
            }
            if v == 0:
                out[v]["psi_vs_stage2_max"] = float(np.max(np.abs(sol.psi - eq.psi)))
            if v == 4 and jac is None:
                jac = {
                    "D3_V4_solution": _jac_check(dd, X, bc(v)),
                    "D3_V4_initial": _jac_check(dd, X_eq, bc(v)),
                }
    return out, jac


def j_short(v):
    w = np.sqrt(2.0 * 2.0 * max(VBI - v - 2.0, 0.0))
    xp = xn = w / 2.0
    return DELTA**2 * np.expm1(v) * (1.0 / (L - xp) + 1.0 / (L - xn))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=REPO / "results" / "stage3")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "manifest.json").write_text(
        json.dumps(build_manifest(REPO, sys.argv, [0]), indent=2) + "\n"
    )

    d1, jac1 = run_d1()
    fwd = list(range(0, 11))
    rev = [0, -1, -2, -3, -4, -5]
    d3, jac3 = _diode(DDParams(delta=DELTA, recombination=False), fwd, rev)
    d4, _ = _diode(DDParams(delta=DELTA, tau_n=1e3, tau_p=1e3), list(range(0, 6)), [0])
    jac = {**jac1, **jac3}

    ideality = {
        f"{v}-{v + 1}": float(np.log(d3[v + 1]["j_mean"] / d3[v]["j_mean"])) for v in range(4, 8)
    }
    ratio_short = {v: d3[v]["j_mean"] / j_short(v) for v in range(2, 9)}
    ratio_rev = {v: d3[v]["j_mean"] / j_short(v) for v in (-2, -5)}
    ideality_d4 = {
        f"{v}-{v + 1}": float(np.log(d4[v + 1]["j_mean"] / d4[v]["j_mean"])) for v in range(1, 4)
    }

    all_runs = [*d1["newton"].values(), *d3.values(), *d4.values()]
    crit = {
        "D1": all(all(1.8 <= q <= 2.2 for q in o[-3:]) for o in d1["orders"].values()),
        "D2": all(v <= 1e-6 for v in jac.values()),
        "D3a": d3[0]["psi_vs_stage2_max"] <= 1e-10
        and max(d3[0]["jn_max_abs"], d3[0]["jp_max_abs"]) <= 1e-12,
        "D3b": all(d3[v]["j_rel_spread"] <= 1e-6 for v in d3 if abs(v) >= 3),
        "D3c": all(0.97 <= q <= 1.03 for q in ideality.values()),
        "D3d": all(0.8 <= q <= 1.2 for q in ratio_short.values()),
        "D5": all(r["positive"] for r in all_runs),
        "D6": all(r["converged"] and r["iterations"] <= 60 for r in all_runs),
    }
    crit = {k: bool(v) for k, v in crit.items()}

    d7 = None
    devsim_file = args.out / "devsim.json"
    if devsim_file.exists():
        ds = json.loads(devsim_file.read_text())
        rows = {int(r["k"]): r for r in ds.get("results", ds.get("sweep", []))}
        d7 = {
            k: d3[k]["j_mean"] / rows[k]["j_scaled"]
            for k in range(2, 9)
            if k in rows and rows[k].get("j_scaled")
        }

    res = {
        "D1": d1,
        "jacobian_rel_err": jac,
        "D3": {str(k): v for k, v in d3.items()},
        "D3c_ideality": ideality,
        "D3d_ratio_to_short_diode": {str(k): v for k, v in ratio_short.items()},
        "D3e_reverse_ratio": {str(k): v for k, v in ratio_rev.items()},
        "D4": {str(k): v for k, v in d4.items()},
        "D4_ideality": ideality_d4,
        "D7_ratio_to_devsim": d7,
        "Vbi": VBI,
    }
    summary = {
        "criteria": crit,
        "gate_passed": all(crit.values()),
        "reported": {
            "D3e": res["D3e_reverse_ratio"],
            "D4_ideality": ideality_d4,
            "D7": d7 if d7 is not None else "devsim.json not present",
        },
    }
    (args.out / "results.json").write_text(json.dumps(res, indent=2, default=float) + "\n")
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
