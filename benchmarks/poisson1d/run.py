"""Run the pre-registered 1D Poisson benchmark (docs/benchmarks/poisson_1d.md).

Writes raw JSON to --out:
  fv.json                 reference-solver errors and observed orders (R1-R5c)
  pinn_<case>_seed<k>.json one file per PINN run (N1-N3)
  summary.json            criterion-by-criterion pass/fail
  manifest.json           reproducibility manifest

Each PINN run executes in its own subprocess so that peak resident memory
(ru_maxrss) is attributable to that run alone.
"""

from __future__ import annotations

import argparse
import json
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from wbgtcad.fv.poisson1d import solve, uniform_mesh
from wbgtcad.manifest import build_manifest
from wbgtcad.problems.poisson1d import CASES, evaluate, linear_charge_free
from wbgtcad.verification import max_error, observed_orders, rel_l2_error

REPO = Path(__file__).resolve().parents[2]
LEVELS = [32, 64, 128, 256, 512]
SEEDS = [0, 1, 2, 3, 4]
PINN_CASES = ["P1", "P2", "P4"]
EVAL_POINTS = 2001


def _dump(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=2) + "\n")


def run_fv() -> dict:
    out: dict = {"levels": LEVELS, "cases": {}}
    for name in ("P1", "P2", "P4"):
        prob = CASES[name]()
        e_inf, e_l2, t = [], [], []
        for n in LEVELS:
            t0 = time.perf_counter()
            sol = solve(prob, uniform_mesh(n))
            t.append(time.perf_counter() - t0)
            ue = evaluate(prob.u_exact, sol.x)
            e_inf.append(max_error(sol.u, ue))
            e_l2.append(rel_l2_error(sol.x, sol.u, ue))
        out["cases"][name] = {
            "e_inf": e_inf,
            "e_l2_rel": e_l2,
            "orders_inf": observed_orders(e_inf),
            "solve_seconds": t,
        }
    p3 = CASES["P3"]()
    out["cases"]["P3"] = {
        str(n): max_error((s := solve(p3, uniform_mesh(n))).u, evaluate(p3.u_exact, s.x))
        for n in (63, 64)
    }
    lin = linear_charge_free()
    s = solve(lin, uniform_mesh(64))
    out["cases"]["linear"] = max_error(s.u, evaluate(lin.u_exact, s.x))

    # R5b discrete conservation and R5c boundary-flux error, P1.
    p1 = CASES["P1"]()
    flux_err = []
    for n in LEVELS:
        x = uniform_mesh(n)
        s = solve(p1, x)
        h = np.diff(x)
        vol = np.zeros_like(x)
        vol[:-1] += 0.5 * h
        vol[1:] += 0.5 * h
        fv = evaluate(p1.f, x)[1:-1] * vol[1:-1]
        if n == LEVELS[-1]:
            out["R5b_rel_imbalance"] = abs(s.flux[-1] - s.flux[0] - fv.sum()) / np.abs(fv).sum()
        # exact flux -u'(x) with u = sin(pi x) + x
        exact0, exact1 = -(np.pi + 1.0), -(-np.pi + 1.0)
        flux_err.append(max(abs(s.flux[0] - exact0), abs(s.flux[-1] - exact1)))
    out["R5c_boundary_flux_err"] = flux_err
    out["R5c_orders"] = observed_orders(flux_err)
    return out


def run_single_pinn(case: str, seed: int, adam_steps: int | None, lbfgs: int | None) -> dict:
    from wbgtcad.pinn.poisson1d import PINNConfig, train

    kw = {"seed": seed}
    if adam_steps is not None:
        kw["adam_steps"] = adam_steps
    if lbfgs is not None:
        kw["lbfgs_max_iter"] = lbfgs
    prob = CASES[case]()
    res = train(prob, PINNConfig(**kw))
    xe = np.linspace(0.0, 1.0, EVAL_POINTS)
    ue = evaluate(prob.u_exact, xe)
    up = res.predict(prob, xe)
    ref = solve(prob, uniform_mesh(512))
    up_nodes = res.predict(prob, ref.x)
    return {
        "case": case,
        "seed": seed,
        "config": res.config,
        "e_l2_rel_exact": rel_l2_error(xe, up, ue),
        "e_inf_exact": max_error(up, ue),
        "e_l2_rel_vs_fv512": rel_l2_error(ref.x, up_nodes, ref.u),
        "final_residual_mse": res.final_loss,
        "train_seconds": res.train_seconds,
        "loss_history_every_100": res.loss_history[::100],
        "n_loss_evals": len(res.loss_history),
        "peak_rss_mib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0,
    }


def summarize(fv: dict, pinn: list[dict]) -> dict:
    c = fv["cases"]
    crit = {
        "R1_P1": all(1.9 <= p <= 2.1 for p in c["P1"]["orders_inf"][-3:]),
        "R1_P2": all(1.9 <= p <= 2.1 for p in c["P2"]["orders_inf"][-3:]),
        "R2": c["P1"]["e_inf"][-1] <= 1e-5,
        "R3": c["P2"]["e_inf"][-1] <= 5e-5,
        "R4": all(v <= 1e-12 for v in c["P3"].values()),
        "R5_linear": c["linear"] <= 1e-12,
        "R5b": fv["R5b_rel_imbalance"] <= 1e-10,
    }
    p1 = [r for r in pinn if r["case"] == "P1"]
    crit["N1"] = len(p1) == len(SEEDS) and all(r["e_l2_rel_exact"] <= 1e-3 for r in p1)
    gating = crit
    return {"criteria": crit, "gate_passed": all(gating.values())}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=REPO / "results" / "poisson1d")
    ap.add_argument("--cases", nargs="+", default=PINN_CASES)
    ap.add_argument("--seeds", nargs="+", type=int, default=SEEDS)
    ap.add_argument("--adam-steps", type=int, default=None, help="override (non-benchmark runs only)")
    ap.add_argument("--lbfgs-iter", type=int, default=None, help="override (non-benchmark runs only)")
    ap.add_argument("--single", nargs=2, metavar=("CASE", "SEED"), help=argparse.SUPPRESS)
    args = ap.parse_args()

    if args.single:
        rec = run_single_pinn(args.single[0], int(args.single[1]), args.adam_steps, args.lbfgs_iter)
        print(json.dumps(rec))
        return

    args.out.mkdir(parents=True, exist_ok=True)
    overridden = args.adam_steps is not None or args.lbfgs_iter is not None
    manifest = build_manifest(REPO, sys.argv, args.seeds)
    manifest["pre_registered_configuration"] = not overridden
    _dump(args.out / "manifest.json", manifest)

    fv = run_fv()
    _dump(args.out / "fv.json", fv)
    print(f"FV done: P1 e_inf(N=512) = {fv['cases']['P1']['e_inf'][-1]:.3e}", flush=True)

    pinn = []
    for case in args.cases:
        for seed in args.seeds:
            cmd = [sys.executable, __file__, "--single", case, str(seed)]
            if args.adam_steps is not None:
                cmd += ["--adam-steps", str(args.adam_steps)]
            if args.lbfgs_iter is not None:
                cmd += ["--lbfgs-iter", str(args.lbfgs_iter)]
            proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
            rec = json.loads(proc.stdout.strip().splitlines()[-1])
            _dump(args.out / f"pinn_{case}_seed{seed}.json", rec)
            pinn.append(rec)
            print(
                f"PINN {case} seed {seed}: rel L2 = {rec['e_l2_rel_exact']:.3e}, "
                f"{rec['train_seconds']:.1f} s, peak RSS {rec['peak_rss_mib']:.0f} MiB",
                flush=True,
            )

    summary = summarize(fv, pinn)
    summary["pre_registered_configuration"] = not overridden
    _dump(args.out / "summary.json", summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
