"""Stage-15b benchmark (docs/benchmarks/stage15_caveat_closure.md, section 15b): F1-F5 (gating), R1-R2.

Field-line ionization-integral re-evaluation of the avalanche BV (no optimization): the initial
design and the final designs of stages 13b, 13c, 13d, on the corrected (stage-12c) off-state model,
on h (level 0) and h/2 (level 1).  Final thetas are read from results/stage13{b,c,d}/results.json.
The vertical-cut reference (R1) is recomputed here with the SAME corrected evaluator and the stage-13b
dense set (49 lines), so that cuts and field lines see identical electrostatics (the stage-13b stored
cut BVs used the old off-state rule).

Study choices (fixed before the run; the spec leaves them free):
  * N_s = 400 contour seeds + the junction crossings of the 49 dense vertical cuts;
  * RK4 step ds = 0.00125 um (F2 halves it).  At ds = 0.01, 0.005 and 0.0025 um a few lines orbit a
    sink of the nodal field in the p-body (|E| stays above the low-field threshold within ds of the
    sink) without meeting either stopping rule, so F5 fails there; none do at 0.00125 / 0.000625 um
    in the pre-run diagnosis.  The choice is made on that termination ground only (recorded in
    results.json "step_choice_note"), not on any BV outcome;
  * march step 20 V, bracket tolerance 0.02 V, secant clip [2 %, 98 %]: those of dense_bv.

Writes <out>/{results,summary,manifest}.json (default results/stage15b) and figures
stage15b_field_lines.png, stage15b_bv_comparison.png (repo figures/ for the default output, else
<out>/figures).  --dev runs a reduced pipeline check (h only, few seeds) and is never an official
result.
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

from widegap.adjoint_sic import EvalError
from widegap.adjoint_sic_corner import DenseCuts, dense_bv
from widegap.fieldline import (
    E_FRAC,
    STEP_LIMIT,
    TOL_V,
    V_STEP_V,
    FieldLineSet,
    diode_states,
    dmos_states,
    junction_contour,
    refine_bv,
)
from widegap.fv.avalanche import avalanche_bv
from widegap.manifest import build_manifest
from widegap.offstate import CorrectedCornerAdjoint, UniformPN2D
from widegap.sic import SiC4H, TracedRegistry

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[2]
DEFAULT_OUT = REPO / "results" / "stage15b"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


M = _load(
    "stage13b_run", REPO / "benchmarks/stage13b/run.py"
)  # constants only (CornerAdjoint unused)
DESIGN = M.DESIGN

# ---- study parameters ------------------------------------------------------------------------
N_S = 400
DS_UM = 0.00125
TOL = {"F1": 1e-3, "F2": 1e-4, "F3": 2e-3, "F4": 0.03}
LEVELS = (0, 1)
DIODE_ND_CM = 1e16
N_VIZ_LINES = 70
DESIGNS = {  # tag -> (results dir, key)
    "initial": ("stage13d", "initial"),
    "13b": ("stage13b", "final"),
    "13c": ("stage13c", "final"),
    "13d": ("stage13d", "final"),
}


def load_thetas():
    th = {}
    for tag, (d, key) in DESIGNS.items():
        r = json.loads((REPO / "results" / d / "results.json").read_text())
        th[tag] = np.asarray(r[key]["theta"], float)
    return th


def make_ev(level):
    sc = SiC4H(TracedRegistry.from_json())
    return CorrectedCornerAdjoint(sc, DESIGN, level, r_um=M.R_UM, cuts_um=(), tau_v=M.TAU_V)


def um_per_unit(sc):
    return float(sc.l_d) * 1e4


# ============================================================ parallel jobs (module level)
def job_fl(args):
    """Field-line BV of one design on one mesh (+ limiting line, + plot data if ``viz``)."""
    level, theta, tag, n_s, ds_um, viz = args
    t0 = time.perf_counter()
    ev = make_ev(level)
    sc = ev.sc
    th = np.asarray(theta, float)
    ev.set_theta(th)
    u = um_per_unit(sc)
    xs = M.dense_positions(sc, ev)
    fls = FieldLineSet.build(ev.mesh, ev.nd, ev.na, ev.ch, u, cut_x_um=xs, n_s=n_s, ds_um=ds_um)
    out = {"tag": tag, "level": level, "n_s": n_s, "ds_um": ds_um, "info": fls.info,
           "n_nodes": int(ev.mesh.n_nodes), "n_seeds": int(len(fls.seeds))}  # fmt: skip
    try:
        states, at = dmos_states(ev, th, fls, v_step_v=V_STEP_V)
        r = refine_bv(states, at, fls, sc.volt_s, sc.volt_phys, tol_v=TOL_V)
        out.update(r)
        o = fls.evaluate(at(r["bv_s"]), keep_lines=True)
        out["limiting_line"] = fls.line_record(o, o["argmax"])
        out["imax_at_bv_state"] = o["imax"]
        out["step_limit_lines_at_bv_state"] = int(np.sum(o["codes"] == STEP_LIMIT))
        out["n_lines_I_gt_0.5_imax"] = int(np.sum(o["I"] > 0.5 * o["imax"]))
        if viz:
            fld = o["field"]
            e_mv = np.hypot(fld.ex, fld.ey) * float(sc.e_unit) * 1e-6
            k_all = np.argsort(fls.seeds[:, 0] * 1e3 + fls.seeds[:, 1])
            pick = k_all[np.linspace(0, len(k_all) - 1, N_VIZ_LINES).astype(int)]
            chains = junction_contour(fld.x, fld.y, _net_sic(ev, fld))
            out["viz"] = {
                "x_um": (fld.x * u).tolist(),
                "y_um": (fld.y * u).tolist(),
                "E_MV_cm": e_mv.tolist(),
                "lines_um": [(o["lines"][k] * u).tolist() for k in pick],
                "limiting_index": int(o["argmax"]),
                "limiting_line_um": (o["lines"][o["argmax"]] * u).tolist(),
                "contour_um": [(c * u).tolist() for c in chains],
                "bv": r["bv"],
            }
    except EvalError as exc:
        out.update({"converged": False, "error": str(exc), "bv": float("nan")})
    out["stats"] = fls.stats
    out["n_stagnation_stops"] = fls.stats["stops"]["stagnation"]
    out["seconds"] = time.perf_counter() - t0
    return out


def _net_sic(ev, fld):
    mesh = ev.mesh
    j = mesh.y <= 1e-12 * max(1.0, float(np.max(np.abs(mesh.y))))
    return (ev.nd - ev.na).reshape(mesh.nx + 1, mesh.ny + 1)[:, j]


def job_dense(args):
    """Dense vertical-cut BV (49 lines) of one design on one mesh, corrected off state."""
    level, theta, tag = args
    t0 = time.perf_counter()
    ev = make_ev(level)
    sc = ev.sc
    xs = M.dense_positions(sc, ev)
    d = dense_bv(ev, np.asarray(theta, float), DenseCuts(ev, xs), tol_v=0.02,
                 loose_tol_v=M.LOOSE_TOL_V.get(level, 3.0))  # fmt: skip
    bv = d["bv"]
    jm = int(np.nanargmin(bv))
    return {"tag": tag, "level": level, "x_um": xs.tolist(), "bv": bv.tolist(),
            "bv_min": float(bv[jm]), "x_argmin_um": float(xs[jm]),
            "resolved_all": bool(d["resolved_all"]), "n_unbracketed": int(d["n_unbracketed"]),
            "seconds": time.perf_counter() - t0}  # fmt: skip


def job_diode(args):
    """F1: laterally uniform p+n (12c-C5): field-line BV vs the stage-12a 1D BV on the same mesh."""
    level, n_s, ds_um = args
    t0 = time.perf_counter()
    sc = SiC4H(TracedRegistry.from_json())
    ch = sc.niwa()
    pn = UniformPN2D(sc, DIODE_ND_CM, level=level)
    b1 = float(sc.volt_phys(avalanche_bv(pn.diode, ch, "n", v_step=pn.v_est / 3.0,
                                         v_max=6.0 * pn.v_est, xtol=1e-12).bv))  # fmt: skip
    fls = FieldLineSet.build(pn.mesh, pn.nd, pn.na, ch, um_per_unit(sc), cut_x_um=[0.0],
                             n_s=n_s, ds_um=ds_um)  # fmt: skip
    out = {"level": level, "n_s": n_s, "ds_um": ds_um, "bv_1d_V": b1, "info": fls.info}
    try:
        states, at = diode_states(pn, fls, pn.v_est / 3.0, 6.0 * pn.v_est)
        r = refine_bv(states, at, fls, sc.volt_s, sc.volt_phys, tol_v=TOL_V)
        o = fls.evaluate(at(r["bv_s"]), keep_lines=True)
        out.update(r)
        out["rel_diff"] = abs(r["bv"] - b1) / b1
        out["max_dx_over_lines_um"] = float(
            max(np.ptp(p[:, 0]) for p in o["lines"]) * um_per_unit(sc))  # fmt: skip
    except EvalError as exc:
        out.update({"converged": False, "error": str(exc), "bv": float("nan"), "rel_diff": np.inf})
    out["stats"] = fls.stats
    out["seconds"] = time.perf_counter() - t0
    return out


# ============================================================ plots
def plot_all(figdir, res, label):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figdir = Path(figdir)
    figdir.mkdir(parents=True, exist_ok=True)
    foot = lambda fig: fig.text(0.5, 0.005, label, ha="center", fontsize=7, color="0.4")  # noqa: E731

    # 1: field lines over |E| for the 13d initial and final designs (h)
    panels = [(t, res["fl"].get(f"{t}_L0")) for t in ("initial", "13d")]
    panels = [(t, r["viz"]) for t, r in panels if r is not None and "viz" in r]
    if panels:
        fig, axs = plt.subplots(1, len(panels), figsize=(5.2 * len(panels), 6.2),
                                constrained_layout=True, squeeze=False)  # fmt: skip
        vmax = max(np.max(v["E_MV_cm"]) for _, v in panels)
        for ax, (tag, v) in zip(axs[0], panels, strict=True):
            x, y, e = np.array(v["x_um"]), np.array(v["y_um"]), np.array(v["E_MV_cm"])
            pm = ax.pcolormesh(x, y, e.T, shading="gouraud", cmap="viridis", vmin=0, vmax=vmax)
            for p in v["lines_um"]:
                p = np.array(p)
                ax.plot(p[:, 0], p[:, 1], "-", color="white", lw=0.6, alpha=0.8)
            for c in v["contour_um"]:
                c = np.array(c)
                ax.plot(c[:, 0], c[:, 1], "--", color="#E69F00", lw=1.0)
            pl = np.array(v["limiting_line_um"])
            ax.plot(pl[:, 0], pl[:, 1], "-", color="#D55E00", lw=2.0, label="limiting line")
            ax.set_title(f"{'initial' if tag == 'initial' else 'stage 13d final'}: V = {v['bv']:.0f} V",
                         fontsize=10)  # fmt: skip
            ax.set_xlabel("x  [um]")
            ax.set_ylabel("y  [um]")
            ax.set_aspect("equal")
            ax.legend(fontsize=8, loc="lower right")
        fig.colorbar(pm, ax=axs[0].tolist(), shrink=0.7, label="|E|  [MV/cm]")
        foot(fig)
        fig.savefig(figdir / "stage15b_field_lines.png", dpi=130)
        plt.close(fig)

    # 2: BV comparison, cut vs field line, per design and mesh
    tags = [t for t in DESIGNS if any(f"{t}_L{lv}" in res["fl"] for lv in LEVELS)]
    fig, ax = plt.subplots(figsize=(9.0, 4.8), constrained_layout=True)
    bars = [("cut", 0, "#0072B2", "//"), ("line", 0, "#D55E00", "//"),
            ("cut", 1, "#0072B2", ""), ("line", 1, "#D55E00", "")]  # fmt: skip
    present = [b for b in bars if any(f"{t}_L{b[1]}" in res["fl"] for t in tags)]
    w = 0.8 / max(len(present), 1)
    for j, (kind, lv, col, hatch) in enumerate(present):
        vals = []
        for t in tags:
            src = res["dense"] if kind == "cut" else res["fl"]
            r = src.get(f"{t}_L{lv}")
            vals.append(np.nan if r is None else (r["bv_min"] if kind == "cut" else r["bv"]))
        xs = np.arange(len(tags)) + (j - (len(present) - 1) / 2) * w
        name = f"{'vertical cut' if kind == 'cut' else 'field line'}, {'h' if lv == 0 else 'h/2'}"
        ax.bar(xs, vals, w * 0.95, color=col, hatch=hatch, edgecolor="white", label=name)
        for xx, vv in zip(xs, vals, strict=True):
            if np.isfinite(vv):
                ax.text(xx, vv, f"{vv:.0f}", ha="center", va="bottom", fontsize=7)
    ax.set_xticks(range(len(tags)), ["initial" if t == "initial" else f"stage {t}" for t in tags])
    ax.set_ylabel("avalanche BV  [V]")
    ax.set_title("BV: dense vertical cuts vs field lines (corrected off state)")
    ax.legend(fontsize=8, ncol=2)
    ax.grid(True, axis="y", alpha=0.3)
    foot(fig)
    fig.savefig(figdir / "stage15b_bv_comparison.png", dpi=130)
    plt.close(fig)


# ============================================================ main
def rel(a, b):
    return float(abs(a - b) / abs(b))


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--workers", type=int, default=min(14, os.cpu_count() or 1))
    ap.add_argument("--dev", action="store_true", help="reduced pipeline check (h only, few seeds)")
    ap.add_argument("--n-seeds", type=int, default=None, help="contour seeds N_s (default 400)")
    args = ap.parse_args(argv)
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    (out / "manifest.json").write_text(
        json.dumps(build_manifest(REPO, sys.argv, []), indent=2) + "\n"
    )
    n_s = args.n_seeds or (40 if args.dev else N_S)
    levels = (0,) if args.dev else LEVELS
    th = load_thetas()
    names = ("initial", "13d") if args.dev else tuple(DESIGNS)
    t_start = time.perf_counter()

    jobs = []  # (kind, key, fn, args)
    for t in names:
        for lv in levels:
            jobs.append(
                (
                    "fl",
                    f"{t}_L{lv}",
                    job_fl,
                    (lv, th[t], t, n_s, DS_UM, t in ("initial", "13d") and lv == 0),
                )
            )
            jobs.append(("dense", f"{t}_L{lv}", job_dense, (lv, th[t], t)))
    jobs.append(
        ("fl", "initial_L0_dsHalf", job_fl, (0, th["initial"], "initial", n_s, DS_UM / 2, False))
    )
    for t in ("initial", "13d"):
        jobs.append(("fl", f"{t}_L0_Ns2", job_fl, (0, th[t], t, 2 * n_s, DS_UM, False)))
    for lv in levels:
        jobs.append(("diode", f"L{lv}", job_diode, (lv, n_s, DS_UM)))
    if args.dev:  # one h/2 field-line run: timing and F4 plumbing only
        jobs.append(("fl", "initial_L1", job_fl, (1, th["initial"], "initial", n_s, DS_UM, False)))
    # longest first: the h/2 jobs dominate
    jobs.sort(key=lambda j: -(1 if "_L1" in j[1] or j[1] == "L1" else 0))

    res: dict = {"fl": {}, "dense": {}, "diode": {}}
    ctx = mp.get_context("fork")
    with ProcessPoolExecutor(args.workers, mp_context=ctx) as pool:
        futs = {pool.submit(fn, a): (kind, key) for kind, key, fn, a in jobs}
        for f in as_completed(futs):
            kind, key = futs[f]
            res[kind][key] = f.result()
            print(f"done {kind} {key}: {res[kind][key]['seconds']:.0f}s", flush=True)

    fl, dn = res["fl"], res["dense"]
    # ---- criteria -----------------------------------------------------------------------------
    crit: dict = {}
    f1 = {k: {"bv_fieldline_V": v["bv"], "bv_1d_V": v["bv_1d_V"], "rel_diff": v["rel_diff"],
              "max_line_dx_um": v.get("max_dx_over_lines_um")} for k, v in res["diode"].items()}  # fmt: skip
    crit["F1"] = {"values": f1, "tol": TOL["F1"],
                  "pass": bool(all(v["rel_diff"] <= TOL["F1"] for v in f1.values()))}  # fmt: skip
    a, b = fl["initial_L0"], fl["initial_L0_dsHalf"]
    r2 = rel(a["bv"], b["bv"])
    crit["F2"] = {"bv_ds_V": a["bv"], "bv_ds_half_V": b["bv"], "ds_um": [a["ds_um"], b["ds_um"]],
                  "rel_diff": r2, "tol": TOL["F2"], "pass": bool(r2 <= TOL["F2"])}  # fmt: skip
    f3 = {}
    for t in ("initial", "13d"):
        f3[t] = {"bv_Ns_V": fl[f"{t}_L0"]["bv"], "bv_2Ns_V": fl[f"{t}_L0_Ns2"]["bv"],
                 "rel_diff": rel(fl[f"{t}_L0"]["bv"], fl[f"{t}_L0_Ns2"]["bv"])}  # fmt: skip
    crit["F3"] = {"values": f3, "N_s": [n_s, 2 * n_s], "tol": TOL["F3"],
                  "pass": bool(all(v["rel_diff"] <= TOL["F3"] for v in f3.values()))}  # fmt: skip
    if "initial_L1" in fl:
        r4 = rel(fl["initial_L0"]["bv"], fl["initial_L1"]["bv"])
        crit["F4"] = {"bv_h_V": fl["initial_L0"]["bv"], "bv_h2_V": fl["initial_L1"]["bv"],
                      "rel_diff": r4, "tol": TOL["F4"], "pass": bool(r4 <= TOL["F4"])}  # fmt: skip
    else:
        crit["F4"] = {"pass": None, "note": "dev run: no h/2"}
    all_runs = [*fl.values(), *res["diode"].values()]
    crit["F5"] = {
        "n_root_finds": len(all_runs),
        "n_not_converged": int(sum(not r.get("converged", False) for r in all_runs)),
        "step_limit_lines": int(sum(r["stats"]["stops"]["step_limit"] for r in all_runs)),
        "lines_traced": int(sum(r["stats"]["n_lines_traced"] for r in all_runs)),
        "stops_total": {k: int(sum(r["stats"]["stops"][k] for r in all_runs))
                        for k in ("left_sic", "low_field", "stagnation", "step_limit")},
    }  # fmt: skip
    crit["F5"]["n_stagnation_stops"] = int(sum(r["stats"]["stops"]["stagnation"] for r in all_runs))
    crit["F5"]["info_step_limit_lines_at_bv_states"] = int(
        sum(r.get("step_limit_lines_at_bv_state", 0) for r in fl.values())
    )  # information only; the gate counts every evaluated state (as written)
    crit["F5"]["pass"] = bool(
        crit["F5"]["n_not_converged"] == 0 and crit["F5"]["step_limit_lines"] == 0
    )

    # ---- reported -----------------------------------------------------------------------------
    r1 = {}
    for t in names:
        for lv in levels:
            k = f"{t}_L{lv}"
            f, d = fl[k], dn[k]
            ll = f.get("limiting_line", {})
            r1[k] = {"bv_fieldline_V": f["bv"], "bv_cut_min_V": d["bv_min"],
                     "fieldline_minus_cut_V": f["bv"] - d["bv_min"],
                     "fieldline_over_cut_minus_1": f["bv"] / d["bv_min"] - 1.0,
                     "limiting_cut_x_um": d["x_argmin_um"], "limiting_line": ll,
                     "cut_resolved_all": d["resolved_all"]}  # fmt: skip
    r2d = {}
    for lv in levels:
        for t in names:
            if t == "initial":
                continue
            k0, k1 = f"initial_L{lv}", f"{t}_L{lv}"
            r2d[f"{t}_L{lv}"] = {
                "fieldline_final_minus_initial_V": fl[k1]["bv"] - fl[k0]["bv"],
                "cut_final_minus_initial_V": dn[k1]["bv_min"] - dn[k0]["bv_min"],
                "fieldline_ratio": fl[k1]["bv"] / fl[k0]["bv"],
            }  # fmt: skip
    gate = None if args.dev else bool(all(crit[k]["pass"] for k in ("F1", "F2", "F3", "F4", "F5")))
    res["config"] = {
        "dev": bool(args.dev), "N_s": n_s, "ds_um": DS_UM, "e_frac": E_FRAC, "march_step_V": V_STEP_V,
        "bracket_tol_V": TOL_V, "tolerances": TOL, "levels": list(levels), "designs": list(names),
        "thetas": {t: th[t].tolist() for t in names}, "n_dense_cuts": M.N_DENSE,
        "offstate_model": "stage 12c (phi_n = V_D w harmonic, phi_p = 0), theta = 0 classification; "
                          "used for ALL designs including the stage-13b final (optimized on the old rule)",
        "diode_N_D_cm3": DIODE_ND_CM, "rounding_length_um": M.R_UM,
        "step_choice_note": "ds = 0.01, 0.005, 0.0025 um leave a few lines orbiting a p-body sink of the nodal field without meeting a stopping rule (13d initial, V = 180-200 V states); none at 0.00125 / 0.000625 um (40 seeds, all marched states). Chosen on this termination ground only.",
        "F4_denominator": "h/2 value", "F2_denominator": "ds/2 value",
        "cut_crossing_seeds": "junction crossings of the 49 dense vertical cuts of stage 13b",
    }  # fmt: skip
    res["runtime_seconds"] = time.perf_counter() - t_start
    summary = {"dev": bool(args.dev), "gate_passed": gate,
               "F1": crit["F1"]["pass"], "F2": crit["F2"]["pass"], "F3": crit["F3"]["pass"],
               "F4": crit["F4"]["pass"], "F5": crit["F5"]["pass"], "n_stagnation_stops": crit["F5"]["n_stagnation_stops"], "criteria": crit,
               "R1": r1, "R2": r2d, "runtime_seconds": res["runtime_seconds"]}  # fmt: skip
    res["criteria"], res["R1"], res["R2"] = crit, r1, r2d

    def dump(name, obj):
        (out / name).write_text(json.dumps(obj, indent=2, default=float) + "\n")

    # keep the (large) plot payload out of results.json
    viz = {k: v.pop("viz") for k, v in fl.items() if "viz" in v}
    dump("results.json", res)
    dump("summary.json", summary)
    for k, v in fl.items():
        if k in viz:
            v["viz"] = viz[k]
    figdir = REPO / "figures" if out == DEFAULT_OUT else out / "figures"
    plot_all(figdir, res, f"stage 15b{' (DEV)' if args.dev else ''}  N_s={n_s}  ds={DS_UM} um")
    print(json.dumps({k: summary[k] for k in ("dev", "gate_passed", "F1", "F2", "F3", "F4", "F5")}))


if __name__ == "__main__":
    main()
