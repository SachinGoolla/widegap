"""Stage-15d benchmark (docs/benchmarks/stage15_caveat_closure.md, section 15d): L1 (gating), L2, R3.

Local-field mesh adequacy of the stage-13d field measures E_ox, W1, W2, W3 at V_blk = 400 V on the
meshes h, h/2, h/4 (levels 0, 1, 2) for the 13d initial (theta = 0) and final designs (theta from
results/stage13d/results.json).  Re-evaluation only: no optimization, no 13d definition changed.

* R3 (diagnosis): position of the hard-max edge of every measure on every mesh and design, its
  distance to the fixed features, and the pinned classification (definitions: widegap.localfield).
* Replacement measure Ebar_k: area-weighted mean of |E| over the fixed square of side L = 0.1 um
  centred at the hard-max position of the 13d INITIAL design on h (the centre is fixed for all
  designs and meshes), clipped to the oxide (E_ox) or the SiC (W1-W3).
* L1 (gating): each Ebar_k of the 13d initial design on h, h/2, h/4 has observed order >= 0.8 and
  the finest value is within 2 % of its Richardson value.
* L2 (reported): D6 re-evaluated with Ebar_k: does any measure increase from the initial to the
  final design on h and h/2 (relative tolerance 1e-6, as 13d D6)?

Also recorded (context, not criteria): the hard maxima and the smooth maxima S_kappa on the three
meshes, the number of edges of each measure and the smooth-max bias bound ln(N)/kappa.

Writes <out>/{results,summary,manifest}.json (default results/stage15d) and figures
stage15d_*.png (repo figures/ for the default output, else <out>/figures).
"""

# ruff: noqa: E402
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import importlib.util
import json
import multiprocessing as mp
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

from widegap import localfield as lf
from widegap.adjoint_sic_oxide import KAPPA_PER_MV_CM, MEASURES, V_BLK_V, WINDOWS_UM
from widegap.manifest import build_manifest

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[2]
DEFAULT_OUT = REPO / "results" / "stage15d"
R13D = REPO / "results" / "stage13d"
LEVELS = (0, 1, 2)  # h, h/2, h/4
D6_RTOL = 1e-6  # as stage 13d D6
REGION = {"E_ox": "oxide", "W1": "sic", "W2": "sic", "W3": "sic"}


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


M13D = _load("stage13d_run", REPO / "benchmarks/stage13d/run.py")  # fixed name: fork pickling
TH0, DESIGN = M13D.TH0, M13D.DESIGN


# ============================================================ job (module level, forked)
def job_state(args):
    """One design on one mesh at V_blk: hard-max records, smooth maxima and (if ``centres`` given or
    ``centres_from_self``) the area-mean measures."""
    level, theta, tag, centres, centres_from_self = args
    t0 = time.perf_counter()
    ev = M13D.make_ev13d(level, with_cuts=False)
    o = ev.fields(np.asarray(theta, float), grad=False, v_volt=V_BLK_V)
    psi = o["X"][0::3]
    um = float(ev.sc.l_d) * 1e4
    geo = lf.Geometry.from_device(ev.sc, DESIGN, ev.prm, M13D.M.R_UM)
    feats = lf.Features(geo)
    recs, means = {}, {}
    for k in MEASURES:
        recs[k] = lf.hard_max_record(ev.measures[k], psi, ev.mesh, um, feats, WINDOWS_UM.get(k))
        recs[k]["smooth_max_MV_cm"] = o[k]["S"]
        recs[k]["n_edges"] = int(len(ev.measures[k].a))
        recs[k]["smooth_bias_bound_MV_cm"] = float(np.log(len(ev.measures[k].a)) / KAPPA_PER_MV_CM)
    cen = (
        {k: (recs[k]["x_um"], recs[k]["y_um"]) for k in MEASURES} if centres_from_self else centres
    )
    if cen is not None:
        for k in MEASURES:
            means[k] = lf.mean_field(ev.mesh, psi, ev._c, um, cen[k], REGION[k])
    return {"tag": tag, "level": level, "n_nodes": int(ev.mesh.n_nodes),
            "mesh_cells": [int(ev.mesh.nx), int(ev.mesh.ny)], "records": recs, "means": means,
            "centres_um": {k: list(v) for k, v in cen.items()} if cen is not None else None,
            "geometry_um": {f: getattr(geo, f) for f in geo.__dataclass_fields__},
            "contours": {"n_source": feats.src_pts[:: max(1, len(feats.src_pts) // 300)].tolist(),
                         "p_body": feats.body_pts[:: max(1, len(feats.body_pts) // 400)].tolist()}
            if level == 0 else None,
            "stats": dict(ev.stats), "seconds": time.perf_counter() - t0}  # fmt: skip


# ============================================================ analysis
def analyse(states):
    """R3 / L1 / L2 from ``states[(tag, level)]`` (dicts returned by :func:`job_state`)."""
    out = {}
    # ---- R3 diagnosis ----
    r3 = {}
    for tag in ("initial", "final"):
        r3[tag] = {}
        for k in MEASURES:
            recs = [states[(tag, lev)]["records"][k] for lev in LEVELS]
            r3[tag][k] = {
                "per_mesh": [
                    {
                        kk: r[kk]
                        for kk in (
                            "x_um",
                            "y_um",
                            "edge_dir",
                            "cell_um",
                            "hard_max_MV_cm",
                            "smooth_max_MV_cm",
                            "nearest_feature",
                            "nearest_distance_um",
                            "nearest_distance_cells",
                            "near_features",
                            "distance_um",
                        )
                    }
                    for r in recs
                ],
                **lf.classify_pinned(recs),
            }
    out["R3"] = r3
    # ---- L1 on the initial design ----
    l1 = {k: lf.l1_assess([states[("initial", lev)]["means"][k]["mean_MV_cm"] for lev in LEVELS])
          for k in MEASURES}  # fmt: skip
    for k in MEASURES:
        l1[k]["area_fraction"] = states[("initial", 0)]["means"][k]["area_fraction"]
        l1[k]["centre_um"] = states[("initial", 0)]["centres_um"][k]
    out["L1"] = {"per_measure": l1, "pass": bool(all(v["pass"] for v in l1.values())),
                 "failed": [k for k, v in l1.items() if not v["pass"]]}  # fmt: skip
    # ---- L2 ----
    l2 = {}
    for lev in (0, 1):
        l2[str(lev)] = {}
        for k in MEASURES:
            i = states[("initial", lev)]["means"][k]["mean_MV_cm"]
            f = states[("final", lev)]["means"][k]["mean_MV_cm"]
            l2[str(lev)][k] = {"initial": i, "final": f, "rel_change": f / i - 1.0,
                               "increases": bool(f > i * (1.0 + D6_RTOL))}  # fmt: skip
    out["L2"] = {"per_mesh_level": l2,
                 "any_increase": bool(any(v["increases"] for lv in l2.values() for v in lv.values())),
                 "increasing": sorted({f"{k}@L{lev}" for lev, lv in l2.items()
                                       for k, v in lv.items() if v["increases"]})}  # fmt: skip
    # ---- context: hard vs smooth vs mean convergence of the initial design ----
    ctx = {}
    for k in MEASURES:
        ctx[k] = {
            "hard_max": lf.l1_assess([states[("initial", lev)]["records"][k]["hard_max_MV_cm"] for lev in LEVELS]),
            "smooth_max": lf.l1_assess([states[("initial", lev)]["records"][k]["smooth_max_MV_cm"] for lev in LEVELS]),
            "n_edges": [states[("initial", lev)]["records"][k]["n_edges"] for lev in LEVELS],
            "smooth_bias_bound_MV_cm": [states[("initial", lev)]["records"][k]["smooth_bias_bound_MV_cm"] for lev in LEVELS],
            "smooth_minus_hard_MV_cm": [states[("initial", lev)]["records"][k]["smooth_max_MV_cm"]
                                        - states[("initial", lev)]["records"][k]["hard_max_MV_cm"] for lev in LEVELS],
        }  # fmt: skip
    out["context_initial"] = ctx
    return out


def check_against_13d(states):
    """Smooth maxima of the initial design on h, h/2, h/4 against results/stage13d/summary.json R1."""
    p = R13D / "summary.json"
    if not p.exists():
        return None
    r1 = json.loads(p.read_text())["values"]["R1"]
    rec = {}
    for tag in ("initial", "final"):
        for k in MEASURES:
            ref = r1[tag][k]["values_MV_cm_h_h2_h4"]
            got = [states[(tag, lev)]["records"][k]["smooth_max_MV_cm"] for lev in LEVELS]
            rec[f"{tag}_{k}"] = {"stage13d": ref, "here": got,
                                 "max_rel_diff": float(max(abs(a / b - 1.0) for a, b in zip(got, ref, strict=True)))}  # fmt: skip
    return {"per_entry": rec, "max_rel_diff": max(v["max_rel_diff"] for v in rec.values())}


# ============================================================ figures
def plot_all(figdir, res, states, label):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figdir = Path(figdir)
    figdir.mkdir(parents=True, exist_ok=True)
    foot = lambda fig: fig.text(0.5, 0.005, label, ha="center", fontsize=7, color="0.4")  # noqa: E731
    col = {"E_ox": "#0072B2", "W1": "#D55E00", "W2": "#009E73", "W3": "#CC79A7"}
    mk = {0: "o", 1: "s", 2: "^"}
    s0 = states[("initial", 0)]
    geo = s0["geometry_um"]
    # 1: hard-max positions on the cross-section
    fig, ax = plt.subplots(figsize=(8.6, 5.2), constrained_layout=True)
    for name, c in (("n_source", "0.35"), ("p_body", "0.65")):
        pts = np.asarray(s0["contours"][name])
        ax.plot(
            pts[:, 0],
            pts[:, 1],
            ".",
            ms=1.5,
            color=c,
            label=f"N_D = N_A contour ({name.replace('_', ' ')})",
        )
    for k, ((x0, x1), (y0, y1)) in WINDOWS_UM.items():
        ax.add_patch(
            plt.Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, ec=col[k], lw=0.8, ls="--")
        )
        ax.text(x0, y1 + 0.01, k, color=col[k], fontsize=8)
    for tag, fill in (("initial", True), ("final", False)):
        for lev in LEVELS:
            for k in MEASURES:
                r = states[(tag, lev)]["records"][k]
                ax.plot(
                    r["x_um"],
                    r["y_um"],
                    mk[lev],
                    ms=7 - lev,
                    color=col[k],
                    mfc=col[k] if fill else "none",
                    mew=1.2,
                )
    for k, c in s0["centres_um"].items():
        ax.add_patch(
            plt.Rectangle((c[0] - 0.05, c[1] - 0.05), 0.1, 0.1, fill=False, ec=col[k], lw=1.2)
        )
    for name, (x, y) in (("gate edge", (geo["gate_x1"], geo["t_ox"])), ("contact start", (geo["src_contact_x0"], 0.0)),
                         ("body contact", (geo["x_bc"], 0.0))):  # fmt: skip
        ax.plot(x, y, "k*", ms=8)
        ax.annotate(name, (x, y), fontsize=7, xytext=(3, 6), textcoords="offset points")
    ax.axhline(0.0, color="k", lw=0.5)
    ax.set_xlim(-0.1, 3.0)
    ax.set_ylim(-1.45, 0.25)
    ax.set_xlabel("x [um]")
    ax.set_ylabel("y [um]")
    ax.set_title("hard-max edge positions at 400 V (filled: initial, open: final; o h, s h/2, ^ h/4)\n"
                 "boxes: windows; small squares: fixed L = 0.1 um centres; stars: interface features", fontsize=9)  # fmt: skip
    hs = [plt.Line2D([], [], marker="o", ls="", color=col[k], label=k) for k in MEASURES]
    ax.legend(handles=hs + ax.get_legend_handles_labels()[0][:2], fontsize=7, loc="lower right")
    foot(fig)
    fig.savefig(figdir / "stage15d_hardmax_positions.png", dpi=130)
    plt.close(fig)
    # 2: convergence of hard max, smooth max and Ebar (initial design)
    fig, axs = plt.subplots(1, 4, figsize=(12.5, 3.9), constrained_layout=True)
    hrel = np.array([1.0, 0.5, 0.25])
    for ax, k in zip(axs, MEASURES, strict=True):
        hard = [states[("initial", lev)]["records"][k]["hard_max_MV_cm"] for lev in LEVELS]
        smooth = [states[("initial", lev)]["records"][k]["smooth_max_MV_cm"] for lev in LEVELS]
        mean = [states[("initial", lev)]["means"][k]["mean_MV_cm"] for lev in LEVELS]
        ax.plot(hrel, hard, "o-", color="0.2", label="hard max")
        ax.plot(hrel, smooth, "s--", color="0.55", label="smooth max (13d)")
        ax.plot(hrel, mean, "^-", color=col[k], label="area mean $\\bar E$ (L = 0.1 um)")
        l1 = res["L1"]["per_measure"][k]
        ax.set_title(f"{k}: L1 {'pass' if l1['pass'] else 'FAIL'}"
                     + (f" (p = {l1['order']:.2f})" if l1["order"] is not None else " (not monotone)"), fontsize=9)  # fmt: skip
        ax.set_xscale("log", base=2)
        ax.invert_xaxis()
        ax.set_xlabel("h / h$_0$")
        ax.grid(True, alpha=0.3)
    axs[0].set_ylabel("field [MV/cm]")
    axs[0].legend(fontsize=7)
    foot(fig)
    fig.savefig(figdir / "stage15d_convergence.png", dpi=130)
    plt.close(fig)


# ============================================================ main
def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--workers", type=int, default=min(6, os.cpu_count() or 1))
    args = ap.parse_args(argv)
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    man = build_manifest(REPO, sys.argv, [])
    (out / "manifest.json").write_text(json.dumps(man, indent=2) + "\n")
    t_start = time.perf_counter()
    th_final = np.asarray(json.loads((R13D / "results.json").read_text())["final"]["theta"], float)
    res: dict = {"config": {
        "V_blk_V": V_BLK_V, "kappa_per_MV_cm": KAPPA_PER_MV_CM, "levels": list(LEVELS),
        "windows_um": {k: [list(v[0]), list(v[1])] for k, v in WINDOWS_UM.items()},
        "side_um": lf.SIDE_UM, "near_cells": lf.NEAR_CELLS, "on_boundary_tol_um": lf.ON_BOUNDARY_TOL_UM,
        "order_min": lf.ORDER_MIN, "rich_tol": lf.RICH_TOL, "d6_rtol": D6_RTOL,
        "theta_initial": TH0.tolist(), "theta_final": th_final.tolist(),
        "centre_rule": "hard-max edge midpoint of the 13d initial design on h, per measure; fixed for all designs/meshes",
        "regions": REGION,
    }}  # fmt: skip
    states: dict = {}
    # phase 1: initial design on h defines the centres (serial, short)
    s0 = job_state((0, TH0, "initial", None, True))
    states[("initial", 0)] = s0
    print(f"initial L0 done {s0['seconds']:.0f}s", flush=True)
    centres = {k: tuple(v) for k, v in s0["centres_um"].items()}
    jobs = [(1, TH0, "initial"), (2, TH0, "initial"), (0, th_final, "final"), (1, th_final, "final"),
            (2, th_final, "final")]  # fmt: skip
    ctx = mp.get_context("fork")
    with ProcessPoolExecutor(args.workers, mp_context=ctx) as pool:
        futs = {
            pool.submit(job_state, (lev, th, tag, centres, False)): (tag, lev)
            for lev, th, tag in jobs
        }
        for f in as_completed(futs):
            r = f.result()
            states[futs[f]] = r
            print(f"done {futs[f]} {r['seconds']:.0f}s", flush=True)
    res["states"] = {f"{tag}_L{lev}": {kk: v for kk, v in s.items() if kk != "contours"}
                     for (tag, lev), s in sorted(states.items())}  # fmt: skip
    ana = analyse(states)
    ana["consistency_with_13d_R1"] = check_against_13d(states)
    res.update(ana)
    label = f"stage 15d | commit {str(man.get('git_commit'))[:8]} | {man['timestamp_utc'][:19]}Z"
    figdir = (REPO / "figures") if out.resolve() == DEFAULT_OUT.resolve() else out / "figures"
    fig_err = None
    try:
        plot_all(figdir, res, states, label)
    except Exception as exc:  # figures are optional
        fig_err = f"{type(exc).__name__}: {exc}"
    res["figures"] = {"dir": str(figdir), "error": fig_err}
    res["runtime_seconds"] = time.perf_counter() - t_start
    n_failed = sum(s["stats"]["failed_evals"] for s in states.values())
    summary = {
        "criteria": {"L1": ana["L1"]["pass"]},
        "gate_passed_15d_L1": ana["L1"]["pass"],
        "values": {
            "L1": {k: {"order": v["order"], "rel_dist_to_richardson": v["rel_dist"], "pass": v["pass"],
                       "reason": v["reason"]} for k, v in ana["L1"]["per_measure"].items()},
            "L1_failed_measures": ana["L1"]["failed"],
            "L2_any_increase": ana["L2"]["any_increase"], "L2_increasing": ana["L2"]["increasing"],
            "L2": ana["L2"]["per_mesh_level"],
            "R3_pinned": {tag: {k: {"pinned": v["pinned"], "pinning_features": v["pinning_features"]}
                                for k, v in d.items()} for tag, d in ana["R3"].items()},
            "failed_evals": n_failed,
            "consistency_with_13d_R1_max_rel_diff": (ana["consistency_with_13d_R1"] or {}).get("max_rel_diff"),
        },
        "runtime_seconds": res["runtime_seconds"],
    }  # fmt: skip
    (out / "results.json").write_text(json.dumps(res, indent=2, default=float) + "\n")
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    print(json.dumps(summary["criteria"], indent=2))
    print(f"runtime {res['runtime_seconds']:.1f} s")


if __name__ == "__main__":
    main()
