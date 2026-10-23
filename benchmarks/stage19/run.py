"""Stage-19 benchmark (docs/benchmarks/stage19_punchthrough_check.md): P1-P3 (gating), R1-R2.

Independent 1D drift-diffusion punch-through check of the study DMOSFET (stage-13d device at
theta = 0: ``OxideAdjoint`` -> ``CorrectedCornerAdjoint`` -> ``RoundedSiCDMOS``, rounding length
0.1 um, the 13d mesh levels).  Writes <out>/{results,summary,manifest}.json, <out>/stderr.log
(default results/stage19) and figures stage19_*.png (repo figures/ for the default output, else
<out>/figures).  ``--dev`` refuses an output under results/.

Choices fixed BEFORE any official run (recorded in results.json["config"]):

* device / columns: the 2D device nodes at the mesh column nearest x = 2.0 um (1.998 um) and at
  x = 1.6 um (the worst column of the 17c scoping, exactly a mesh column); vertical lines from the
  surface (source side) to the drain contact, nodes = the 2D mesh node rows y <= 0.  Mesh h =
  stage-13d level 0 column grid, h/2 = level 1; level 2 (h/4) is a supplementary low-bias
  convergence trend only (not gated).
* doping: the cell-centre functions of ``RoundedDMOSFET`` (same masks) evaluated at the column's x,
  averaged over the 1D control volumes (``widegap.punchthrough.column_doping``).
* physics: stage-9a edge core (``AnisoMobilityDD``) on a 1-column mesh: incomplete ionization
  (Goetz N donor, Koizumi Al acceptor, device levels), Caughey-Thomas velocity saturation with the
  device's doping-dependent vertical (E || c) electron mobility, no impact ionization, no
  recombination (R = 0), constant interface-free bulk.  Contacts: top = n+ source (psi_neutral, phi
  = 0), bottom = n+ drain (psi_neutral + V_D, phi = V_D).
* configurations: ``floating`` (primary, gated: no body contact) and ``tied`` (supplementary:
  hole-only contact phi_p = 0 at the node nearest mid-body depth).
* continuation: one ``continue_bias`` march per run through a FIXED grid, 1 V steps to 150 V then
  5 V steps to 400 V, starting from the equilibrium state; stops at the first target that does not
  converge (reached range = last converged target).  Convergence is the unchanged stage-8.1
  Newton test (every free row |F|/S <= 1e2 eps and last step <= 1e-8).
* J_n: electron flux on the edge just below the potential minimum of the barrier window (surface
  to d_body + 0.5 um), in A/cm^2 (physical unit; all log tests use ln|J_n| in A/cm^2).
* barrier: psi(surface node) - min psi over the window [V_T].
* P1a: |ln J_h - ln J_h/2| / |ln J_h/2| <= 1e-3 at V_D = 10 V, both columns (the absolute
  |dln J| is reported as well, and the same ratio in scaled units).
* P1b: |barrier_1D(10 V) - barrier_12c(10 V)| <= 1 V_T for every (column, mesh); the stage-12c
  value is the zero-current off-state (phi_p = 0, body tied) of the same 2D column from
  ``solve_at_voltage``; 5/20/40 V are tabulated.  Gated on the floating configuration; the tied
  one is reported.
* V_pt: smallest V_D with |J_n| > 1e3 |J_n(10 V)| (log-linear interpolation on the grid).
* P3: |V_pt(h/2) - V_pt(h)| / V_pt(h/2) <= 5 % at both columns (floating).
* P2: floating, both columns, h and h/2: every target up to 150 V converged, and no "On entry to"
  line in stderr.log.
* R1/R2: avalanche BV = results/stage12c/results.json jobs_raw["mode|harmonic"]["dense_bv_min_V"]
  (dense minimum, level 0) and results/stage13d/results.json initial.bv_obj (soft-min) /
  dense.initial_L0/L1.bv_min; verdict "confirmed" iff the finest-mesh floating V_pt of the worst
  column is below the dense-minimum avalanche BV.
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
import resource
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

from widegap import linalg
from widegap.manifest import build_manifest
from widegap.punchthrough import (
    Column1D,
    failure_record,
    first_crossing,
    interp_at,
    punch_voltage,
    zero_current_barrier,
)

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[2]
DEFAULT_OUT = REPO / "results" / "stage19"
BV_12C = REPO / "results" / "stage12c" / "results.json"
BV_13D = REPO / "results" / "stage13d" / "results.json"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


M13D = _load("stage13d_run", REPO / "benchmarks/stage13d/run.py")  # fixed name: fork-pool pickling
TH0 = M13D.TH0

COLUMNS_UM = (2.0, 1.6)
LEVELS = (0, 1)  # h, h/2 (gated)
LEVEL_SUPP = 2  # h/4, low bias only
BODIES = ("floating", "tied")
GRID_V = [float(v) for v in range(1, 151)] + [float(v) for v in range(155, 401, 5)]
GRID_SUPP_V = [float(v) for v in range(1, 41)]
V_REF = 10.0  # V_pt reference and P1 bias
V_ZC = (5.0, 10.0, 20.0, 40.0)  # stage-12c zero-current comparison biases
V_REACH = 150.0  # P2
TOL = {"P1a_rel_log": 1e-3, "P1b_VT": 1.0, "P3_rel": 0.05, "J_factor": 1.0e3}
J_LEVELS = (1e-3, 1.0, 1e3)  # A/cm^2, supplementary onset voltages
BARRIER_LOW_VT = 10.0  # comparison metric with the 17c scoping (not V_pt)
SCOPING_ESTIMATE_V = (100.0, 150.0)


def rss_mib():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


# ============================================================ one run (fork pool)
def job_run(args):
    """One (column, level, body) march; returns the per-voltage series and the failure record."""
    x_um, level, body, grid, zc = args
    t0 = time.perf_counter()
    nf0 = linalg.nonfinite_calls
    ev = M13D.make_ev13d(level, with_cuts=False)
    col = Column1D(ev, x_um, TH0, body)
    eq = col.equilibrium()
    series = []
    fail = None
    v_last = 0.0
    if not eq.report.converged:
        fail = {"stage": "equilibrium", "target_V": 0.0}
    else:
        res = col.march(eq.X, 0.0, grid)
        for v, sol, its, solves in res:
            if not sol.report.converged:
                fail = failure_record(col, sol, v)
                break
            o = col.observe(sol)
            o.update({"V": float(v), "iterations": int(its), "solves": int(solves)})
            series.append(o)
            v_last = float(v)
    out = {
        "key": f"x{x_um:g}|L{level}|{body}",
        "x_um_requested": x_um,
        "x_um_column": col.x_um,
        "level": level,
        "body": body,
        "n_nodes": col.n,
        "body_tie_node_depth_um": col.body_node_um,
        "t_um": col.t_um.tolist(),
        "nd_cm3": (col.nd * col.sc.n0).tolist(),
        "na_cm3": (col.na * col.sc.n0).tolist(),
        "equilibrium": col.observe(eq),
        "series": series,
        "V_last_converged": v_last,
        "failure": fail,
        "reached_target": fail is None,
        "stats": dict(col.stats),
        "nonfinite_solves": linalg.nonfinite_calls - nf0,
    }
    if zc:
        out["zero_current_12c"] = zero_current_barrier(col, TH0, V_ZC)
    out["seconds"] = time.perf_counter() - t0
    return out


# ============================================================ analysis
def arr(run, key):
    return np.array([s[key] for s in run["series"]], dtype=float)


def at_v(run, v, key):
    return interp_at(arr(run, "V"), arr(run, key), v)


def curve_metrics(run):
    """V_pt, barrier at V_pt, onset voltages for one run."""
    v = arr(run, "V")
    j = np.abs(arr(run, "Jn_A_cm2"))
    bar = arr(run, "barrier_VT")
    pt = punch_voltage(v, j, V_REF, TOL["J_factor"])
    out = {"V_pt": pt["V_pt"], "J_ref_10V": pt["J_ref"], "J_threshold": pt["threshold"]}
    if pt["V_pt"] is not None:
        out["barrier_at_V_pt_VT"] = interp_at(v, bar, pt["V_pt"])
        out["J_at_V_pt_A_cm2"] = float(
            np.exp(interp_at(v, np.log(np.maximum(j, 1e-300)), pt["V_pt"]))
        )
        out["barrier_drop_at_V_pt_VT"] = float(interp_at(v, bar, V_REF) - out["barrier_at_V_pt_VT"])
    ln = np.log(np.maximum(j, 1e-300))
    out["onset_V"] = {f"{lv:g}": first_crossing(v, ln, np.log(lv), True) for lv in J_LEVELS}
    out["V_barrier_below_10VT"] = first_crossing(v, bar, BARRIER_LOW_VT, above=False)
    out["V_reached"] = float(v[-1]) if len(v) else 0.0
    return out


def evaluate(runs, bv):
    """P1-P3 and R1/R2 from the finished runs."""
    cols = sorted({r["x_um_requested"] for r in runs.values()}, reverse=True)
    key = lambda x, lv, b: f"x{x:g}|L{lv}|{b}"  # noqa: E731
    vals = {}
    # ---------------- P1a
    p1a = {}
    for body in BODIES:
        rows = []
        for x in cols:
            h, h2 = runs[key(x, 0, body)], runs[key(x, 1, body)]
            jh, jh2 = at_v(h, V_REF, "Jn_A_cm2"), at_v(h2, V_REF, "Jn_A_cm2")
            row = {"x_um": x, "J_h_A_cm2": jh, "J_h2_A_cm2": jh2}
            if jh and jh2:
                lh, lh2 = np.log(abs(jh)), np.log(abs(jh2))
                row.update(
                    abs_dlnJ=float(abs(lh - lh2)),
                    rel_log=float(abs(lh - lh2) / abs(lh2)),
                    rel_J=float(abs(jh / jh2 - 1.0)),
                )
                js = _j_scale(runs)
                ls_h, ls_h2 = np.log(abs(jh) / js), np.log(abs(jh2) / js)
                row["rel_log_scaled_units"] = float(abs(ls_h - ls_h2) / abs(ls_h2))
            r4 = runs.get(key(x, LEVEL_SUPP, body))
            if r4 is not None:
                j4 = at_v(r4, V_REF, "Jn_A_cm2")
                row["J_h4_A_cm2"] = j4
                row["abs_dlnJ_h2_h4"] = float(abs(np.log(abs(jh2) / abs(j4)))) if j4 else None
            rows.append(row)
        p1a[body] = rows
    # ---------------- P1b
    p1b = {}
    for body in BODIES:
        rows = []
        for x in cols:
            for lv in LEVELS:
                r = runs[key(x, lv, body)]
                zc = runs[key(x, lv, "floating")]["zero_current_12c"]
                ent = {"x_um": x, "level": lv, "by_V": {}}
                for v in V_ZC:
                    b1 = at_v(r, v, "barrier_VT")
                    b2 = zc[v] if v in zc else zc[float(v)]
                    ent["by_V"][f"{v:g}"] = {
                        "barrier_1D": b1,
                        "barrier_12c": b2,
                        "diff": None if b1 is None else float(b1 - b2),
                    }
                ent["phi_p_body_V_at_10V"] = at_v(r, V_REF, "phi_p_body_V")
                d10 = ent["by_V"][f"{V_REF:g}"]["diff"]
                ent["diff_10V"] = d10
                ent["ok"] = bool(d10 is not None and abs(d10) <= TOL["P1b_VT"])
                rows.append(ent)
        p1b[body] = rows
    ok_a = all(
        r.get("rel_log") is not None and r["rel_log"] <= TOL["P1a_rel_log"] for r in p1a["floating"]
    )
    ok_b = all(e["ok"] for e in p1b["floating"])
    vals["P1"] = {"P1a": p1a, "P1b": p1b, "P1a_ok_floating": bool(ok_a), "P1b_ok_floating": bool(ok_b),
                  "tied_P1a_ok": bool(all(r.get("rel_log") is not None and r["rel_log"] <= TOL["P1a_rel_log"]
                                          for r in p1a["tied"])),
                  "tied_P1b_ok": bool(all(e["ok"] for e in p1b["tied"]))}  # fmt: skip
    # ---------------- P2
    p2 = {}
    for body in BODIES:
        p2[body] = [
            {
                "x_um": x,
                "level": lv,
                "V_last_converged": runs[key(x, lv, body)]["V_last_converged"],
                "reached_150V": bool(runs[key(x, lv, body)]["V_last_converged"] >= V_REACH),
                "failure": runs[key(x, lv, body)]["failure"],
            }
            for x in cols
            for lv in LEVELS
        ]
    p2_ok = all(e["reached_150V"] for e in p2["floating"]) and bv["on_entry_to_lines"] == 0
    vals["P2"] = {"by_body": p2, "on_entry_to_lines": bv["on_entry_to_lines"], "ok": bool(p2_ok),
                  "grid": "1 V steps to 150 V, 5 V steps to 400 V, one continue_bias march"}  # fmt: skip
    # ---------------- curves, V_pt, P3
    cm = {k: curve_metrics(r) for k, r in runs.items() if r["level"] in LEVELS}
    p3 = {}
    for body in BODIES:
        rows = []
        for x in cols:
            a, b = cm[key(x, 0, body)]["V_pt"], cm[key(x, 1, body)]["V_pt"]
            rel = None if (a is None or b is None) else float(abs(b - a) / b)
            rows.append({"x_um": x, "V_pt_h": a, "V_pt_h2": b, "rel_change": rel,
                         "ok": bool(rel is not None and rel <= TOL["P3_rel"])})  # fmt: skip
        p3[body] = rows
    vals["P3"] = {"by_body": p3, "ok": bool(all(e["ok"] for e in p3["floating"]))}
    vals["curve_metrics"] = cm
    # ---------------- R1 / R2
    fin = {
        body: {x: cm[key(x, 1, body)]["V_pt"] for x in cols} for body in BODIES
    }  # finest gated mesh
    worst = {}
    for body in BODIES:
        cand = [(v, x) for x, v in fin[body].items() if v is not None]
        worst[body] = min(cand) if cand else None
    bv_min = bv["avalanche_dense_min_BV_V"]
    verdict = {}
    for body in BODIES:
        if worst[body] is None:
            verdict[body] = "undetermined (V_pt not reached within the converged range)"
        elif worst[body][0] < bv_min:
            verdict[body] = "confirmed (DD V_pt below the avalanche BV)"
        else:
            verdict[body] = "refuted (DD V_pt not below the avalanche BV)"
    vals["R1"] = {
        "avalanche_BV": bv,
        "scoping_estimate_V": list(SCOPING_ESTIMATE_V),
        "analytic_1D_estimate_V": 150.0,
        "V_pt_finest_mesh": {b: {f"{x:g}": v for x, v in fin[b].items()} for b in BODIES},
        "worst_column_V_pt": {
            b: None if worst[b] is None else {"V_pt": worst[b][0], "x_um": worst[b][1]}
            for b in BODIES
        },
        "note_V_pt_definition": (
            "V_pt as defined (1e3 x J_n(10 V)) marks a barrier drop of about ln(1e3) = 7 V_T from "
            "a ~1e-45 A/cm2 level; it is not an onset of appreciable current.  onset_V (J_n = 1e-3,"
            " 1, 1e3 A/cm2) and V_barrier_below_10VT (the 17c scoping metric) are reported in "
            "curve_metrics."
        ),
    }
    vals["R2"] = {"verdict_by_body": verdict, "verdict": verdict["floating"]}
    return vals


def _j_scale(runs):
    r = next(iter(runs.values()))
    return r.get("j_scale_A_cm2", 1.0)


# ============================================================ figures
def plot_all(figdir: Path, res: dict, label: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figdir.mkdir(parents=True, exist_ok=True)
    runs = res["runs"]
    colors = {"floating": ["#1f77b4", "#9ecae1"], "tied": ["#d62728", "#fc9272"]}
    ls = {2.0: "-", 1.6: "--"}
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    for r in runs.values():
        if r["level"] not in LEVELS or not r["series"]:
            continue
        v, j, b = arr(r, "V"), np.abs(arr(r, "Jn_A_cm2")), arr(r, "barrier_VT")
        lab = f"x={r['x_um_column']:.2f} um, {r['body']}, L{r['level']}"
        kw = {"color": colors[r["body"]][r["level"]], "ls": ls[r["x_um_requested"]], "lw": 1.4}
        ax[0].semilogy(v, np.maximum(j, 1e-60), label=lab, **kw)
        ax[1].plot(v, b, label=lab, **kw)
    for r in runs.values():
        zc = r.get("zero_current_12c")
        if zc and r["body"] == "floating":
            ax[1].plot(list(map(float, zc)), list(zc.values()), "k.", ms=4)
    ax[0].set_xlabel("V_D [V]")
    ax[0].set_ylabel("|J_n| at the barrier [A/cm$^2$]")
    ax[1].set_xlabel("V_D [V]")
    ax[1].set_ylabel("barrier $\\Delta\\psi$ [$V_T$] (dots: stage-12c zero current)")
    ax[1].axhline(BARRIER_LOW_VT, color="gray", lw=0.6)
    ax[0].legend(fontsize=6)
    ax[0].set_title("electron current density")
    ax[1].set_title("source-side barrier")
    fig.suptitle(label, fontsize=8)
    fig.tight_layout()
    fig.savefig(figdir / "stage19_curves.png", dpi=130)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for r in runs.values():
        if r["level"] == 0 and r["series"]:
            ax.semilogy(r["t_um"], np.array(r["nd_cm3"]) + 1.0, color="C0", lw=1.0, alpha=0.5)
            ax.semilogy(r["t_um"], np.array(r["na_cm3"]) + 1.0, color="C3", lw=1.0, alpha=0.5)
    ax.set_xlabel("depth below the surface [um]")
    ax.set_ylabel("N_D (blue), N_A (red) + 1 [cm$^{-3}$]")
    ax.set_title("column doping (x = 2.0 and 1.6 um, level 0)")
    fig.tight_layout()
    fig.savefig(figdir / "stage19_doping.png", dpi=130)
    plt.close(fig)


# ============================================================ main
def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--workers", type=int, default=min(10, os.cpu_count() or 1))
    ap.add_argument("--dev", action="store_true", help="DEV run (refused under results/)")
    args = ap.parse_args(argv)
    out = args.out.resolve()
    if args.dev and (REPO / "results") in (out, *out.parents):
        raise SystemExit("refusing --dev with --out under results/")
    out.mkdir(parents=True, exist_ok=True)
    log_path = out / "stderr.log"
    log_fd = os.open(log_path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
    saved_fd = os.dup(2)
    os.dup2(log_fd, 2)
    try:
        run(args, out, log_path)
    finally:
        sys.stdout.flush()
        sys.stderr.flush()
        os.dup2(saved_fd, 2)
        os.close(log_fd)
        os.close(saved_fd)
        text = log_path.read_text(errors="replace")
        if text:
            sys.stderr.write(text)


def run(args, out: Path, log_path: Path) -> None:
    man = build_manifest(REPO, sys.argv, [])
    man["dev"] = bool(args.dev)
    (out / "manifest.json").write_text(json.dumps(man, indent=2) + "\n")
    t_start = time.perf_counter()
    d12c = json.loads(BV_12C.read_text())
    d13d = json.loads(BV_13D.read_text())
    bv = {
        "avalanche_dense_min_BV_V": d12c["jobs_raw"]["mode|harmonic"]["dense_bv_min_V"],
        "source_dense_min": "results/stage12c/results.json jobs_raw['mode|harmonic']['dense_bv_min_V']"
        " (rounded device r=0.1 um, theta=0, level 0, 49 dense lines)",
        "stage13d_initial_bv_obj_V": d13d["initial"]["bv_obj"],
        "stage13d_initial_dense_bv_min_L0_V": d13d["dense"]["initial_L0"]["bv_min"],
        "stage13d_initial_dense_bv_min_L1_V": d13d["dense"]["initial_L1"]["bv_min"],
        "stage13d_initial_limiting_cut_um": d13d["initial"]["limiting_cut_x_um"],
        "source_13d": "results/stage13d/results.json initial.bv_obj, dense.initial_L0/L1.bv_min",
    }
    jobs = []
    for x in COLUMNS_UM:
        for lv in LEVELS:
            for body in BODIES:
                jobs.append((x, lv, body, GRID_V, body == "floating"))
        jobs.append((x, LEVEL_SUPP, "floating", GRID_SUPP_V, False))
    ctx = mp.get_context("fork")
    runs: dict = {}
    t0 = time.perf_counter()
    with ProcessPoolExecutor(args.workers, mp_context=ctx) as pool:
        futs = {pool.submit(job_run, j): j for j in jobs}
        for f in as_completed(futs):
            r = f.result()
            runs[r["key"]] = r
            print(f"done {r['key']} V_last={r['V_last_converged']:.0f} {r['seconds']:.0f}s",
                  flush=True)  # fmt: skip
    t_par = time.perf_counter() - t0
    sys.stderr.flush()
    stderr = log_path.read_text(errors="replace") if log_path.exists() else ""
    bv["on_entry_to_lines"] = len([ln for ln in stderr.splitlines() if "On entry to" in ln])
    ev_in = runs
    # J scale for the scaled-unit variant of P1a
    sc = M13D.make_ev13d(0, with_cuts=False).sc
    from widegap.punchthrough import Q_E

    js = Q_E * sc.n0 * sc.mu0 * sc.vt / sc.l_d
    for r in ev_in.values():
        r["j_scale_A_cm2"] = float(js)
    vals = evaluate(ev_in, bv)
    crit = {
        "P1": bool(vals["P1"]["P1a_ok_floating"] and vals["P1"]["P1b_ok_floating"]),
        "P2": bool(vals["P2"]["ok"]),
        "P3": bool(vals["P3"]["ok"]),
    }
    res = {
        "config": {
            "device": "stage-13d OxideAdjoint (CorrectedCornerAdjoint -> RoundedSiCDMOS), theta = 0,"
            " r = 0.1 um, DESIGN of benchmarks/stage13d/run.py",
            "columns_um": list(COLUMNS_UM),
            "levels_gated": list(LEVELS),
            "level_supplementary": LEVEL_SUPP,
            "bodies": list(BODIES),
            "grid_V": [GRID_V[0], GRID_V[149], GRID_V[-1]],
            "tolerances": TOL,
            "v_ref_V": V_REF,
            "v_zero_current_V": list(V_ZC),
            "j_scale_A_cm2": float(js),
            "workers": args.workers,
        },
        "bv": bv,
        "runs": runs,
        "evaluation": vals,
        "parallel_seconds": t_par,
        "peak_rss_MiB": rss_mib(),
    }
    res["runtime_seconds"] = time.perf_counter() - t_start
    figdir = (REPO / "figures") if out == DEFAULT_OUT.resolve() else out / "figures"
    label = f"stage 19 | commit {str(man.get('git_commit'))[:8]} | {man['timestamp_utc'][:19]}Z"
    try:
        plot_all(figdir, res, label)
    except Exception as exc:
        res["figures_error"] = f"{type(exc).__name__}: {exc}"
    (out / "results.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    short = {
        k: {kk: vv for kk, vv in r.items() if kk not in ("series", "t_um", "nd_cm3", "na_cm3")}
        for k, r in runs.items()
    }
    summary = {
        "criteria": crit,
        "gate_passed": bool(all(crit.values())),
        "dev": bool(args.dev),
        "values": {**vals, "runs_short": short},
        "runtime_seconds": res["runtime_seconds"],
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    print(json.dumps(crit, indent=2))
    print(f"R2 verdict: {vals['R2']['verdict']}")
    print(f"runtime {res['runtime_seconds']:.1f} s")


if __name__ == "__main__":
    main()
