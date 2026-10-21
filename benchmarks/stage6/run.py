"""Stage-6 benchmark (docs/benchmarks/stage6_optimization.md).

Writes results/stage6/{results,summary,manifest}.json and
figures/stage6_optimized_profile.{png,svg}.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from widegap.fv import reverse_bias as rb
from widegap.manifest import build_manifest
from widegap.optimize import PTDiodeDesign, evaluate, feasible, optimize

REPO = Path(__file__).resolve().parents[2]
E_CRIT, KAPPA = 15.0, 200.0


def grad_check(des, theta, seeds=(0, 1, 2), eps=1e-4):
    _, g, _ = evaluate(des, theta, E_CRIT, KAPPA)
    out = []
    for s in seeds:
        u = np.random.default_rng(s).standard_normal(des.n_modes)
        u /= np.linalg.norm(u)
        fp = evaluate(des, theta + eps * u, E_CRIT, KAPPA, grad=False)[0]
        fm = evaluate(des, theta - eps * u, E_CRIT, KAPPA, grad=False)[0]
        fd = (fp - fm) / (2 * eps)
        out.append({"adjoint": float(g @ u), "fd": fd, "rel_err": abs(g @ u - fd) / abs(fd)})
    return out


def exact_bv(theta, h):
    """Independent exact (non-smoothed) BV: new design/mesh, fresh stage-4 sweep + Brent."""
    des = PTDiodeDesign(h=h)
    c, _ = des.c_cv(np.asarray(theta))
    d = rb.Diode1D(des.x, c, des.delta)
    sw = rb.sweep(d, range(0, 151))
    bv, psi, _ = rb.breakdown_voltage(d, sw, E_CRIT, xtol=1e-12)
    return float(bv), d, psi


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=REPO / "results" / "stage6")
    ap.add_argument("--criteria", choices=["v1", "6b"], default="v1")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "manifest.json").write_text(
        json.dumps(build_manifest(REPO, sys.argv, [0, 1, 2]), indent=2) + "\n"
    )

    des = PTDiodeDesign()
    th0 = np.zeros(des.n_modes)
    theta, rec = optimize(des, th0, max_iter=60, e_crit=E_CRIT, kappa=KAPPA)

    gc_init = grad_check(des, th0)
    gc_final = grad_check(des, theta)

    # Nested-family check: linear gradings theta = (t1, 0, 0, 0, 0)
    grid = []
    for t1 in np.linspace(-2.0, 2.0, 41):
        th = np.zeros(des.n_modes)
        th[0] = t1
        ok = feasible(des, th)[0]
        bvs, zero_bias_breakdown = None, False
        if ok:
            try:
                bvs = evaluate(des, th, E_CRIT, KAPPA, grad=False)[0]
            except ValueError:  # B_s >= 1 already at V_R = 0, i.e. BV_s <= 0
                zero_bias_breakdown = True
        grid.append(
            {
                "theta1": float(t1),
                "feasible": ok,
                "bv_s": bvs,
                "breaks_down_at_zero_bias": zero_bias_breakdown,
            }
        )
    best_lin = max(g["bv_s"] for g in grid if g["bv_s"] is not None)

    # Independent exact-BV re-verification on two meshes
    reverify = {}
    for h in (0.025, 0.0125):
        b0, _, _ = exact_bv(th0, h)
        b1, d1, p1 = exact_bv(theta, h)
        reverify[str(h)] = {"bv_initial": b0, "bv_final": b1, "improvement": b1 - b0}

    acc = rec.bv_s
    crit = {
        "S6a": all(r["rel_err"] <= 1e-5 for r in gc_init + gc_final)
        if args.criteria == "v1"
        else all(r["rel_err"] <= 1e-5 for r in gc_init)
        and all(abs(r["adjoint"] - r["fd"]) <= 1e-6 and abs(r["fd"]) <= 1e-4 for r in gc_final),
        "S6b": all(b > a for a, b in zip(acc[:-1], acc[1:], strict=True)),
        "S6c": max(rec.r_rel_err) <= 1e-12 and min(rec.b_min) >= 0.05 and max(rec.b_max) <= 20.0,
        "S6d": rec.stop_reason in ("gradient", "stagnation", "line search") and len(acc) - 1 <= 60,
        "S6e": all(v["improvement"] > 0 for v in reverify.values())
        and reverify["0.0125"]["improvement"] >= 0.5 * reverify["0.025"]["improvement"],
        "S6f": acc[-1] >= best_lin - 1e-6 * acc[-1],
    }
    crit = {k: bool(v) for k, v in crit.items()}

    res = {
        "theta_final": theta.tolist(),
        "record": rec.__dict__,
        "grad_check_initial": gc_init,
        "grad_check_final": gc_final,
        "linear_grid": grid,
        "best_linear_bv_s": best_lin,
        "reverify_exact_bv": reverify,
        "final_grad_norm": rec.grad_norm[-1],
    }
    summary = {
        "criteria_set": args.criteria,
        "criteria": crit,
        "gate_passed": all(crit.values()),
        "bv_s_initial": acc[0],
        "bv_s_final": acc[-1],
        "best_linear_bv_s": best_lin,
        "linear_grid_zero_bias_breakdowns": sum(g["breaks_down_at_zero_bias"] for g in grid),
        "iterations": len(acc) - 1,
        "evaluations": rec.evaluations,
        "stop_reason": rec.stop_reason,
        "reverify_exact_bv": reverify,
        "worst_grad_rel_err": max(r["rel_err"] for r in gc_init + gc_final),
    }
    (args.out / "results.json").write_text(json.dumps(res, indent=2, default=float) + "\n")
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    plot(des, th0, theta, args.out)


def plot(des, th0, theta, out):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    man = json.loads((out / "manifest.json").read_text())
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    xs = np.linspace(0, des.W, 400)
    for th, lab, col in ((th0, "initial (uniform)", "#0072B2"), (theta, "optimized", "#D55E00")):
        ax[0].plot(xs, des.b(th, xs)[0], color=col, label=lab)
        bv, d, psi = exact_bv(th, des.h)
        xm = 0.5 * (d.x[1:] + d.x[:-1])
        ax[1].plot(xm, np.abs(d.element_field(psi)), color=col, label=f"{lab}: BV = {bv:.1f} V_T")
    ax[0].set_yscale("log")
    ax[0].set_xlabel("x in drift region  [L_D]")
    ax[0].set_ylabel("drift doping b(x) / N0")
    ax[0].set_title(r"Doping at fixed $\int dx/b$")
    ax[0].legend(fontsize=8)
    ax[1].axhline(E_CRIT, color="0.5", ls="--", lw=0.8)
    ax[1].set_xlim(-1, des.W + 1)
    ax[1].set_xlabel("x  [L_D]")
    ax[1].set_ylabel("|E| at field-based BV  [V_T/L_D]")
    ax[1].set_title("Field profile at B = 1 (proxy, no avalanche)")
    ax[1].legend(fontsize=8)
    fig.text(
        0.5,
        0.005,
        f"widegap | commit {man['git_commit'][:7]} | {man['cpu_model']} | CPU only",
        ha="center",
        fontsize=7,
        color="0.4",
    )
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    (REPO / "figures").mkdir(exist_ok=True)
    for ext in ("png", "svg"):
        fig.savefig(REPO / "figures" / f"stage6_optimized_profile.{ext}", dpi=200)


if __name__ == "__main__":
    main()
