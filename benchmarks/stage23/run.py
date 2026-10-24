"""Stage-23 benchmark (docs/benchmarks/stage23_artifact_budget.md): T1-T3 (gating), R1.

Artifact budget: claimed vs real gains/values under each model correction, same device, mesh and
design.  Writes <out>/{results,summary,manifest}.json, <out>/stderr.log (default results/stage23)
and figures/stage23_artifact_budget.png (repo figures/ for the default output, else
<out>/figures).  ``--dev`` refuses an output under results/.

Choices fixed BEFORE any official run (recorded in results.json["config"]):

* AB1: the stage-13 device (sharp stage-12 12b design, ``SiCAdjoint``, mesh level 0, the stage-12
  ORIGINAL off-state model = ``SiCDMOS.offstate_bc`` with the theta = 0 phi_n classification).
  Planar-cut BV exactly as the stage-13 ``reverify`` (cold evaluator, ``bv(grad=False,
  warm=False)``); dense true minimum = ``dense_bv`` on the same evaluator over the 49 stage-13b
  lines x = 0..x_max (tol 0.02 V near the minimum, 1 V elsewhere), the same model.  The x = x_bc
  column is added as a 50th line, report-only (it must reproduce the planar value to the dense
  tolerance).  Dense minimum = minimum over the 49 standard lines.
* AB2: no solves; recorded stage-17 / stage-17b values (dense BV, h).  (b): the off-state BV does
  not depend on the on-state surface classification, so the 17-final BV gain is unchanged; only
  R_on,ch changes (recorded reclassification, h/2).
* AB3: the 0.8 um geometry (stage-20 ``make_ev(0.8, level)``, 12c off-state).  Dense true-min BV
  at h: the recorded values (stage 13c/13d/17b) are used in the table; they are also recomputed
  here (same ``dense_bv`` call as stage 20, theta = the recorded theta) as a cross-check.
  V_pt,es on h and h/2 with the stage-20 ``OffMarch.vpt_es`` (10 V_T level, source columns,
  bisection 0.5 V).  1D tied-body DD (``Column1D``, level h/2, stage-20 DD grid) at the column that
  attains the minimum barrier at V_pt,es on h/2; onset = first J_n >= 1 A/cm^2.  Blocking =
  min(BV, onset).  If the DD march fails before the onset, the entry is "lower bound only"
  (T2) with blocking = min(BV, V reached).
* AB4: stage-16 ``job_ron`` (bulk Hall and Noguchi channel, d_ch = 20 nm) for the 13d initial and
  final designs.  stage 16 recorded levels 0 and 2 only, so level 1 (h/2) is a NEW solve with the
  stage-16 code path; level 0 is re-run and compared with the record as a reproduction check.
* AB5: recorded stage-22 values (1.4 um geometry; the E run is the stage-21 result).
* Metric: artifact fraction 1 - real/claimed, applied to the claimed quantity (gain for AB1, AB2,
  AB5; absolute value for AB3, AB4).
* T1: recomputed planar-cut gain AND both planar BVs vs the stage-13 record, relative <= 1e-6.
* T2: every DD onset reached or the entry flagged "lower bound only"; for every AB3 design
  |V_pt,es(h/2)/V_pt,es(h) - 1| <= 5 % (both beyond the 1200 V ceiling count as equal).
* T3: every job completed without error; every dense BV resolved with no failed evaluation
  (lines that never bracket lie above 2x the first-crossing voltage, so cannot change the
  minimum; their number is reported); the DD equilibrium converged; stage-16 bulk/Picard solves
  converged; no "On entry to" line in stderr.log.  Rejected step-halving Newton trials and the
  DD's first non-converged target after its onset are REPORTED, not gated (precedent: stage 17
  J7, stage 20 G1, stage 22 Q3).
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

from widegap.adjoint_sic import SiCAdjoint
from widegap.adjoint_sic_corner import DenseCuts, dense_bv
from widegap.artifactbudget import Cites, blocking_voltage, entry, rel_change
from widegap.bodyredesign import OffMarch, column_index_matrix, source_columns
from widegap.manifest import build_manifest
from widegap.punchthrough import Column1D, first_crossing
from widegap.sic import SiC4H, TracedRegistry

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[2]
DEFAULT_OUT = REPO / "results" / "stage23"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


M13 = _load("stage13_run", REPO / "benchmarks/stage13/run.py")  # fixed names: fork-pool pickling
M20 = _load("stage20_run", REPO / "benchmarks/stage20/run.py")
M16 = _load("stage16_run", REPO / "benchmarks/stage16/run.py")

TH0 = [0.0] * 5
RES = "results/{}/results.json"
TOL = {"T1_rel": 1e-6, "T2_vpt_es_rel": 0.05}
J_ONSET = 1.0  # A/cm^2
DESIGNS = ("initial", "13", "13c", "13d", "17b")
DESIGN_LABEL = {"initial": "initial", "13": "13 final", "13c": "13c final", "13d": "13d final",
                "17b": "17b final"}  # fmt: skip
THETA_SRC = {  # (results file stage, key path) of each design's theta
    "13": ("stage13", ("final", "theta")),
    "13c": ("stage13c", ("final", "theta")),
    "13d": ("stage13d", ("final", "theta")),
    "17b": ("stage17b", ("final", "theta")),
}
BV_SRC = {  # recorded dense true-min BV on h (12c model, 0.8 um geometry)
    "initial": ("stage13d", ("dense", "initial_L0", "bv_min")),
    "13": ("stage13c", ("dense", "stage13_final_L0", "bv_min")),
    "13c": ("stage13c", ("dense", "final_L0", "bv_min")),
    "13d": ("stage13d", ("dense", "final_L0", "bv_min")),
    "17b": ("stage17b", ("dense", "final_L0", "bv_min")),
}
N_DENSE = M20.N_DENSE
V_CEIL = M20.V_CEIL


def rss_mib():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def theta_of(cites: Cites, name):
    if name == "initial":
        return list(TH0)
    st, keys = THETA_SRC[name]
    return [float(t) for t in cites.get(RES.format(st), *keys)]


# ============================================================ jobs (module level: fork pool)
def job_ab1(args):
    """Stage-13 device, old off-state model: planar-cut BV and dense true minimum (level 0)."""
    name, theta = args
    t0 = time.perf_counter()
    sc = SiC4H(TracedRegistry.from_json())
    ev = SiCAdjoint(sc, M13.DESIGN, 0)
    th = np.asarray(theta, float)
    ev.bv_cache = None
    ev.X_on = None  # cold, exactly as the stage-13 reverify
    ev.X_eq = None
    planar = float(ev.bv(th, grad=False, warm=False)["BV"])
    x_max_um = float(sc.length_cm(ev.prm.x_max)) * 1e4
    x_bc_um = float(sc.length_cm(ev.prm.x_bc)) * 1e4
    xs = np.linspace(0.0, x_max_um, N_DENSE)
    xs_all = np.r_[xs, x_bc_um]
    d = dense_bv(ev, th, DenseCuts(ev, xs_all), tol_v=0.02, loose_tol_v=M20.LOOSE_TOL_V[0])
    bv = np.asarray(d["bv"], float)
    bd = bv[:N_DENSE]
    jm = int(np.nanargmin(bd))
    return {
        "kind": "ab1", "name": name, "theta": list(map(float, theta)),
        "planar_bv": planar, "cut_x_um": float(sc.length_cm(ev.col_x)) * 1e4,
        "dense_bv_min": float(bd[jm]), "x_argmin_um": float(xs[jm]),
        "dense_bv_at_x_bc": float(bv[-1]), "x_bc_um": x_bc_um,
        "bv": bd.tolist(), "x_um": xs.tolist(),
        "n_unbracketed": int(np.isnan(bd).sum()), "resolved_all": bool(d["resolved_all"]),
        "stats": dict(ev.stats), "flips_max": int(ev.stats["flips_max"]),
        "seconds": time.perf_counter() - t0,
    }  # fmt: skip


def job_dense(args):
    """Recompute of the dense true-min BV on the 0.8 um geometry (12c model), level 0."""
    name, theta = args
    t0 = time.perf_counter()
    ev = M20.make_ev(0.8, 0)
    sc = ev.sc
    xs = np.linspace(0.0, float(sc.length_cm(ev.prm.x_max)) * 1e4, N_DENSE)
    th = np.asarray(theta, float)
    d = dense_bv(ev, th, DenseCuts(ev, xs), tol_v=0.02, loose_tol_v=M20.LOOSE_TOL_V[0])
    bv = np.asarray(d["bv"], float)
    jm = int(np.nanargmin(bv))
    return {
        "kind": "dense", "name": name, "bv_min": float(bv[jm]), "x_argmin_um": float(xs[jm]),
        "n_unbracketed": int(np.isnan(bv).sum()), "resolved_all": bool(d["resolved_all"]),
        "stats": dict(ev.stats), "seconds": time.perf_counter() - t0,
    }  # fmt: skip


def dd_tied(ev, x_um, theta):
    col = Column1D(ev, x_um, np.asarray(theta, float), "tied")
    eq = col.equilibrium()
    if not eq.report.converged:
        return {"error": "equilibrium did not converge", "equilibrium_converged": False}
    res = col.march(eq.X, 0.0, M20.DD_GRID_V)
    v, j, bar, fail = [], [], [], None
    for vv, sol, _, _ in res:
        if not sol.report.converged:
            fail = float(vv)
            break
        o = col.observe(sol)
        v.append(float(vv))
        j.append(abs(o["Jn_A_cm2"]))
        bar.append(o["barrier_VT"])
    ln = np.log(np.maximum(np.array(j), 1e-300))
    return {
        "equilibrium_converged": True, "x_um_column": col.x_um, "V": v, "Jn_A_cm2": j,
        "barrier_VT": bar, "V_reached": v[-1] if v else 0.0, "first_failed_target_V": fail,
        "onset_V": first_crossing(v, ln, np.log(J_ONSET), True) if len(v) else None,
        "V_barrier_below_10VT": first_crossing(v, bar, M20.BARRIER_LEVEL_VT, above=False),
    }  # fmt: skip


def job_vpt(args):
    """V_pt,es (zero-current barrier, stage-20 rule) and, on level 1, the tied-body DD onset."""
    name, level, theta = args
    t0 = time.perf_counter()
    ev = M20.make_ev(0.8, level)
    sc = ev.sc
    th = np.asarray(theta, float)
    ev.set_theta(th)
    cols = source_columns(ev)
    idx = column_index_matrix(ev, cols)
    xs_um = [float(sc.length_cm(ev.dev.mesh.x[i])) * 1e4 for i in cols]
    r = OffMarch(ev, th).vpt_es(idx, M20.BARRIER_LEVEL_VT, v_max_v=V_CEIL)
    out = {
        "kind": "vpt", "name": name, "level": level, "V_pt_es": r["V_pt_es"],
        "V_bracket": r.get("V_bracket"), "n_columns": len(cols), "n_solves": r["n_solves"],
        "stats": dict(ev.stats),
    }  # fmt: skip
    if r["V_pt_es"] is not None:
        out["worst_x_um"] = xs_um[r["worst_column_row"]]
    else:
        out["worst_x_um"] = xs_um[r["curve"][-1][2]]
    if level == 1:
        out["dd"] = dd_tied(ev, out["worst_x_um"], th)
    out["seconds"] = time.perf_counter() - t0
    return out


def job_ron16(args):
    """Stage-16 R_on (bulk and channel, d_ch = 20 nm) of the 13d initial / final design."""
    case, level = args
    rec = M16.job_ron((case, level, (20,), M16.registry_json(False)))
    rec["name"] = case
    return rec


def run_job(spec):
    fn = {"ab1": job_ab1, "dense": job_dense, "vpt": job_vpt, "ron16": job_ron16}[spec[0]]
    try:
        return fn(spec[1:])
    except Exception as exc:
        import traceback

        return {"kind": spec[0], "name": spec[1] if isinstance(spec[1], str) else None,
                "level": spec[2] if len(spec) > 2 and isinstance(spec[2], int) else None,
                "error": repr(exc), "traceback": traceback.format_exc()}  # fmt: skip


# ============================================================ analysis
def evaluate(res: dict, c: Cites, n_blas: int):
    errors = [f"{k}: {v['error']}" for k, v in res.items() if v.get("error")]
    table, detail = [], {}

    # ---------------- AB1 + T1
    a1 = {n: res.get(("ab1", n)) for n in ("initial", "13")}
    s13 = "results/stage13/summary.json"
    rec_init = c.ref(s13, "values", "bv_initial_V")
    rec_fin = c.ref(s13, "values", "bv_final_V")
    t1 = {"recorded": {"bv_initial": rec_init, "bv_final": rec_fin}}
    ab1_ok = all(a1[n] and not a1[n].get("error") for n in a1)
    if ab1_ok:
        g_claim = a1["13"]["planar_bv"] - a1["initial"]["planar_bv"]
        g_rec = rec_fin["value"] - rec_init["value"]
        errs = {
            "bv_initial": abs(a1["initial"]["planar_bv"] / rec_init["value"] - 1.0),
            "bv_final": abs(a1["13"]["planar_bv"] / rec_fin["value"] - 1.0),
            "gain": abs(g_claim / g_rec - 1.0),
        }
        t1.update(rel_errors=errs, tol=TOL["T1_rel"], ok=bool(max(errs.values()) <= TOL["T1_rel"]),
                  recomputed={"bv_initial": a1["initial"]["planar_bv"],
                              "bv_final": a1["13"]["planar_bv"], "gain": g_claim})  # fmt: skip
        g_real = a1["13"]["dense_bv_min"] - a1["initial"]["dense_bv_min"]
        table.append(
            entry(
                "AB1",
                "stage 13 (planar cut -> dense true minimum), old off-state, sharp device, h",
                "BV gain initial -> final",
                g_claim,
                g_real,
                "V",
                {
                    "claimed": "recomputed planar-cut BV of the stage-13 initial/final (T1: equals "
                    + f"{rec_init['file']} {rec_init['key']}, {rec_fin['key']})",
                    "real": "recomputed dense true-minimum BV (49 lines) of the same two designs, same "
                    "model; theta: " + RES.format("stage13") + " final/theta",
                },
                f"planar BV {a1['initial']['planar_bv']:.3f} -> {a1['13']['planar_bv']:.3f} V; dense "
                f"min {a1['initial']['dense_bv_min']:.3f} (x={a1['initial']['x_argmin_um']:.2f} um) -> "
                f"{a1['13']['dense_bv_min']:.3f} V (x={a1['13']['x_argmin_um']:.2f} um)",
                claimed_endpoints=[a1["initial"]["planar_bv"], a1["13"]["planar_bv"]],
                real_endpoints=[a1["initial"]["dense_bv_min"], a1["13"]["dense_bv_min"]],
            )
        )
        detail["AB1"] = {n: {k: v for k, v in a1[n].items() if k not in ("bv", "x_um")}
                         for n in a1}  # fmt: skip
    else:
        t1.update(ok=False, error="AB1 job failed")
    # ---------------- AB2 (recorded)
    s17, s17b = RES.format("stage17"), RES.format("stage17b")
    bv17 = c.ref(s17, "R1", "bv_true_min_h", "17")
    bv0 = c.ref(s17, "R1", "bv_true_min_h", "initial")
    gain17 = c.ref(s17, "R1", "bv_gain_h_17")
    bv17b = c.ref(s17b, "R4", "bv_true_min_h", "17b")
    bv0b = c.ref(s17b, "R4", "bv_true_min_h", "initial")
    reclass = c.ref(s17b, "R4", "reclass_rel_change_h2")
    rfro = c.ref(s17b, "R4", "frozen_R_s17_final", "2")
    rrec = c.ref(s17b, "R4", "reclassified_R_s17_final", "2", "R")
    g17 = bv17["value"] - bv0["value"]
    g17b = bv17b["value"] - bv0b["value"]
    table.append(
        entry(
            "AB2a",
            "stage 17 (frozen surface classification), reading (a): 17b gain",
            "BV gain over initial (dense, h)",
            g17,
            g17b,
            "V",
            {"claimed": [gain17, bv17, bv0], "real": [bv17b, bv0b]},
            "17b keeps the surface classification valid (constraint); R_on,ch = R0 in both",
            claimed_endpoints=[bv0["value"], bv17["value"]],
            real_endpoints=[bv0b["value"], bv17b["value"]],
        )
    )
    table.append(
        entry(
            "AB2b",
            "stage 17, reading (b): 17 final under dynamic reclassification",
            "BV gain over initial (dense, h)",
            g17,
            g17,
            "V",
            {"claimed": [gain17], "real": [bv17, bv0, reclass, rfro, rrec]},
            f"BV (off-state) is independent of the on-state surface classification: gain unchanged; "
            f"R_on,ch (h/2) is {reclass['value'] * 100:+.2f} % vs R0, so the design is feasible and "
            "the frozen constraint was slack. The sign of the artifact is model-dependent (this "
            "reading: no loss of gain, extra R_on margin; reading (a): -38.7 % gain)",
            claimed_endpoints=[bv0["value"], bv17["value"]],
            real_endpoints=[bv0["value"], bv17["value"]],
            R_on_rel_change_reclassified_h2=reclass["value"],
            real_not_unique=True,
        )
    )
    # ---------------- AB3
    ab3, t2_rows, t2_ok = {}, {}, True
    for n in DESIGNS:
        st, keys = BV_SRC[n]
        bv_ref = c.ref(RES.format(st), *keys)
        v0, v1 = res.get(("vpt", n, 0)), res.get(("vpt", n, 1))
        dn = res.get(("dense", n))
        if not (v0 and v1 and not v0.get("error") and not v1.get("error")):
            t2_ok = False
            continue
        dd = v1.get("dd", {})
        onset, reached = dd.get("onset_V"), dd.get("V_reached", 0.0)
        blk = blocking_voltage(bv_ref["value"], onset, reached)
        rc = rel_change(v0["V_pt_es"], v1["V_pt_es"])
        mesh_ok = bool(rc is not None and rc <= TOL["T2_vpt_es_rel"])
        t2_rows[n] = {"V_pt_es_h": v0["V_pt_es"], "V_pt_es_h2": v1["V_pt_es"], "rel_change": rc,
                      "mesh_ok": mesh_ok, "dd_onset_V": onset, "dd_V_reached": reached,
                      "dd_first_failed_target_V": dd.get("first_failed_target_V"),
                      "lower_bound_only": blk["lower_bound_only"]}  # fmt: skip
        t2_ok = t2_ok and mesh_ok
        recomputed = None if not dn or dn.get("error") else dn["bv_min"]
        ab3[n] = {
            "bv_recorded": bv_ref, "bv_recomputed": recomputed,
            "bv_recomputed_rel_diff": None if recomputed is None else
            abs(recomputed / bv_ref["value"] - 1.0),
            "V_pt_es_h": v0["V_pt_es"], "V_pt_es_h2": v1["V_pt_es"],
            "worst_x_um_h2": v1["worst_x_um"], "worst_x_um_h": v0["worst_x_um"],
            "dd_x_um_column": dd.get("x_um_column"), "dd_onset_V": onset, "dd_V_reached": reached,
            "dd_first_failed_target_V": dd.get("first_failed_target_V"),
            "dd_V_barrier_below_10VT": dd.get("V_barrier_below_10VT"), **blk,
        }  # fmt: skip
        th_src = (f"theta: {RES.format(THETA_SRC[n][0])} {'/'.join(THETA_SRC[n][1])}"
                  if n != "initial" else "theta = 0")  # fmt: skip
        vpt_txt = "n/a" if v1["V_pt_es"] is None else f"{v1['V_pt_es']:.1f}"
        table.append(
            entry(
                f"AB3[{n}]",
                f"{DESIGN_LABEL[n]} on the 0.8 um geometry (12c model, h)",
                "blocking voltage min(BV, V_pt) vs avalanche BV",
                bv_ref["value"],
                blk["blocking_V"],
                "V",
                {
                    "claimed": bv_ref,
                    "real": f"min(BV, DD onset J_n=1 A/cm2, level h/2, tied body, "
                    f"worst column x={v1['worst_x_um']:.3f} um); {th_src}",
                },
                f"V_pt,es(h/2)={vpt_txt} V (h: "
                f"{'n/a' if v0['V_pt_es'] is None else format(v0['V_pt_es'], '.1f')}); DD onset "
                f"{'not reached (V reached ' + format(reached, '.0f') + ' V)' if onset is None else format(onset, '.2f') + ' V'}"
                f"; limited by {blk['limited_by']}"
                + ("; LOWER BOUND ONLY" if blk["lower_bound_only"] else ""),
                V_pt_es_h=v0["V_pt_es"],
                V_pt_es_h2=v1["V_pt_es"],
                dd_onset_V=onset,
                lower_bound_only=blk["lower_bound_only"],
            )
        )
    t2_dd_flagged = all(r["dd_onset_V"] is not None or r["lower_bound_only"]
                        for r in t2_rows.values())  # fmt: skip
    # ---------------- AB4
    s16 = RES.format("stage16")
    ab4, ab4_repro, ron_ok = {}, {}, True
    for case, bulk_key in (("13d_initial", ("dense", "initial_L1", "R")),
                           ("13d_final", ("dense", "final_L1", "R"))):  # fmt: skip
        r0, r1 = res.get(("ron16", case, 0)), res.get(("ron16", case, 1))
        if not (r0 and r1 and not r0.get("error") and not r1.get("error")):
            ron_ok = False
            continue
        rec0 = c.ref(s16, "runs", f"ron|{case}|0", "channel", "20", "R_on_sp_ohm_cm2")
        rec2 = c.ref(s16, "runs", f"ron|{case}|2", "channel", "20", "R_on_sp_ohm_cm2")
        recb0 = c.ref(s16, "runs", f"ron|{case}|0", "bulk", "R_on_sp_ohm_cm2")
        recb2 = c.ref(s16, "runs", f"ron|{case}|2", "bulk", "R_on_sp_ohm_cm2")
        bulk13d = c.ref(RES.format("stage13d"), *bulk_key)
        ch0 = r0["channel"]["20"]
        ch1 = r1["channel"]["20"]
        okc = bool(r0["bulk_converged"] and r1["bulk_converged"] and ch0.get("picard", {}).get("converged")
                   and ch1.get("picard", {}).get("converged"))  # fmt: skip
        ron_ok = ron_ok and okc
        claimed, real = r1["bulk"]["R_on_sp_ohm_cm2"], ch1["R_on_sp_ohm_cm2"]
        ab4_repro[case] = {
            "L0_channel_rel_diff_vs_record": abs(ch0["R_on_sp_ohm_cm2"] / rec0["value"] - 1.0),
            "L0_bulk_rel_diff_vs_record": abs(r0["bulk"]["R_on_sp_ohm_cm2"] / recb0["value"] - 1.0),
            "L1_bulk_rel_diff_vs_stage13d_dense_R": abs(claimed / bulk13d["value"] - 1.0),
            "recorded_L0_channel": rec0, "recorded_L2_channel": rec2,
            "recorded_L0_bulk": recb0, "recorded_L2_bulk": recb2,
            "stage13d_bulk_L1": bulk13d, "converged": okc,
        }  # fmt: skip
        ab4[case] = {"bulk_h2": claimed, "channel_h2": real, "ratio": real / claimed}
        table.append(
            entry(
                f"AB4[{case}]",
                f"{case} (bulk Hall -> Noguchi channel mobility), h/2",
                "R_on,sp at V_GS=15 V, V_DS=0.1 V",
                claimed,
                real,
                "ohm cm2",
                {
                    "claimed": "stage-16 job_ron bulk at level 1 (new solve; cross-check "
                    + f"{bulk13d['file']} {bulk13d['key']})",
                    "real": "stage-16 job_ron channel d_ch=20 nm at level 1 (new solve; stage 16 recorded "
                    "levels 0 and 2 only: " + f"{rec0['file']} {rec0['key']}, level 2)",
                },
                f"channel/bulk = {real / claimed:.2f}x at h/2 (recorded: L0 {rec0['value'] / recb0['value']:.2f}x, "
                f"L2 {rec2['value'] / recb2['value']:.2f}x)",
                ratio=real / claimed,
            )
        )
    # ---------------- AB5
    s22 = RES.format("stage22")
    bvD = c.ref(s22, "p22b", "dense_bv_min_D")
    bvE = c.ref(s22, "p22b", "dense_bv_min_E")
    bv0_22 = c.ref(s22, "initial_dense_bv_min")
    nitD = c.ref(s22, "p22b", "runs", "D", "nit")
    nitE = c.ref(s22, "p22b", "runs", "E", "nit")
    gD, gE = bvD["value"] - bv0_22["value"], bvE["value"] - bv0_22["value"]
    table.append(
        entry(
            "AB5",
            "stage 22 (simplified D-gradient optimum -> exact E-gradient optimum), 1.4 um geometry",
            "BV gain over initial (dense, h)",
            gD,
            gE,
            "V",
            {"claimed": [bvD, bv0_22], "real": [bvE, bv0_22, c.ref(s22, "p22b", "E_run_source")]},
            f"BV(D) - BV(E) = {bvD['value'] - bvE['value']:+.4f} V; SLSQP iterations D/E = "
            f"{nitD['value']}/{nitE['value']} = {nitD['value'] / nitE['value']:.2f}x (E run from stage 21)",
            claimed_endpoints=[bv0_22["value"], bvD["value"]],
            real_endpoints=[bv0_22["value"], bvE["value"]],
            bv_diff_D_minus_E_V=bvD["value"] - bvE["value"],
            iterations_D=nitD["value"],
            iterations_E=nitE["value"],
            iteration_ratio_D_over_E=nitD["value"] / nitE["value"],
        )
    )
    # ---------------- T3
    unbr = {}
    t3_ok = not errors and n_blas == 0 and ron_ok
    for v in res.values():
        if v.get("error"):
            continue
        if v["kind"] in ("ab1", "dense"):
            unbr[f"{v['kind']}|{v['name']}"] = v["n_unbracketed"]
            t3_ok = t3_ok and bool(v["resolved_all"]) and v["stats"].get("failed_evals", 1) == 0
        if v["kind"] == "vpt" and v["level"] == 1:
            t3_ok = t3_ok and bool(v["dd"].get("equilibrium_converged"))
    t3 = {"ok": bool(t3_ok), "job_errors": errors, "on_entry_to_lines": n_blas,
          "n_unbracketed_lines": unbr, "stage16_solves_converged": bool(ron_ok)}  # fmt: skip
    crit = {"T1": bool(t1.get("ok")), "T2": bool(t2_ok and t2_dd_flagged and len(t2_rows) == len(DESIGNS)),
            "T3": bool(t3_ok)}  # fmt: skip
    extra = {"T1": t1, "T2": t2_rows, "T3": t3, "AB3": ab3, "AB4": ab4, "AB4_checks": ab4_repro,
             "detail": detail}  # fmt: skip
    return table, extra, crit


# ============================================================ figure
def plot_budget(figdir: Path, table: list, label: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figdir.mkdir(parents=True, exist_ok=True)
    byid = {r["id"]: r for r in table}
    c_cl, c_re = "#4c78a8", "#e45756"

    def bars(ax, ids, names, fmt, ylabel, title, logy=False):
        x = np.arange(len(ids))
        w = 0.38
        cl = [byid[i]["claimed"] for i in ids]
        re = [byid[i]["real"] for i in ids]
        ax.bar(x - w / 2, cl, w, color=c_cl, label="claimed")
        ax.bar(x + w / 2, re, w, color=c_re, label="real (after correction)")
        top = max(cl + re)
        for k in range(len(ids)):
            ax.text(x[k] - w / 2, cl[k], fmt.format(cl[k]), ha="center", va="bottom", fontsize=7)
            ax.text(x[k] + w / 2, re[k], fmt.format(re[k]), ha="center", va="bottom", fontsize=7)
        ax.set_xticks(x)
        af = [byid[i]["artifact_fraction"] for i in ids]
        ax.set_xticklabels(
            [
                f"{n}\nartifact {'n/a' if a is None else format(a * 100, '+.1f') + ' %'}"
                for n, a in zip(names, af, strict=True)
            ],
            fontsize=7,
        )
        ax.set_ylabel(ylabel, fontsize=8)
        ax.set_title(title, fontsize=9)
        if logy:
            ax.set_yscale("log")
        else:
            ax.set_ylim(0, top * 1.15)
        ax.spines[["top", "right"]].set_visible(False)

    fig, ax = plt.subplots(1, 3, figsize=(17, 5.2), gridspec_kw={"width_ratios": [4, 5, 2.4]})
    g_ids = [i for i in ("AB1", "AB2a", "AB2b", "AB5") if i in byid]
    g_names = {"AB1": "AB1\nplanar -> dense\n(stage 13)", "AB2a": "AB2a\nstage 17 vs 17b",
               "AB2b": "AB2b\nstage 17 reclassified", "AB5": "AB5\nD -> E gradient\n(1.4 um)"}  # fmt: skip
    bars(
        ax[0],
        g_ids,
        [g_names[i] for i in g_ids],
        "{:.1f}",
        "BV gain over initial [V]",
        "BV gains: claimed vs real",
    )
    ax[0].legend(fontsize=8, frameon=False, loc="upper right")
    b_ids = [f"AB3[{n}]" for n in DESIGNS if f"AB3[{n}]" in byid]
    b_names = [DESIGN_LABEL[n] + ("\n(lower bound)" if byid[f"AB3[{n}]"].get("lower_bound_only") else "")
               for n in DESIGNS if f"AB3[{n}]" in byid]  # fmt: skip
    bars(
        ax[1],
        b_ids,
        b_names,
        "{:.0f}",
        "voltage [V]",
        "AB3: avalanche BV (dense) vs blocking min(BV, V_pt)",
    )
    for k, i in enumerate(b_ids):
        v = byid[i].get("V_pt_es_h2")
        if v is not None:
            ax[1].plot(
                [k - 0.5, k + 0.5], [v, v], "k--", lw=1, label="V_pt,es (h/2)" if k == 0 else None
            )
    ax[1].legend(fontsize=8, frameon=False, loc="upper right")
    r_ids = [i for i in ("AB4[13d_initial]", "AB4[13d_final]") if i in byid]
    r_tab = {i: dict(byid[i], claimed=byid[i]["claimed"] * 1e3, real=byid[i]["real"] * 1e3)
             for i in r_ids}  # fmt: skip
    saved = dict(byid)
    byid.update(r_tab)
    bars(
        ax[2],
        r_ids,
        ["13d initial", "13d final"],
        "{:.2f}",
        "R_on,sp [m$\\Omega$ cm$^2$] (h/2)",
        "AB4: bulk Hall vs channel mobility",
    )
    byid.update(saved)
    fig.suptitle(
        "Stage 23 artifact budget: claimed vs real under one model correction at a time",
        fontsize=10,
    )
    fig.text(0.99, 0.005, label, ha="right", fontsize=6)
    fig.tight_layout(rect=(0, 0.02, 1, 0.96))
    fig.savefig(figdir / "stage23_artifact_budget.png", dpi=140)
    plt.close(fig)


# ============================================================ main
def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--workers", type=int, default=min(12, os.cpu_count() or 1))
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
    cites = Cites(REPO)
    thetas = {n: theta_of(cites, n) for n in DESIGNS}
    specs = [("ab1", n, thetas["initial" if n == "initial" else "13"]) for n in ("initial", "13")]
    for n in DESIGNS:
        specs += [("dense", n, thetas[n]), ("vpt", n, 0, thetas[n]), ("vpt", n, 1, thetas[n])]
    for case in ("13d_initial", "13d_final"):
        specs += [("ron16", case, 0), ("ron16", case, 1)]

    def prio(s):  # longest first
        if s[0] == "vpt" and s[2] == 1:
            return 0
        if s[0] == "ron16" and s[2] == 1:
            return 1
        return 2

    specs.sort(key=prio)
    res = {}
    ctx = mp.get_context("fork")
    with ProcessPoolExecutor(args.workers, mp_context=ctx) as pool:
        futs = {pool.submit(run_job, s): s for s in specs}
        for f in as_completed(futs):
            s, r = futs[f], f.result()
            if s[0] in ("ab1", "dense"):
                key = (s[0], s[1])
            else:  # vpt / ron16: (kind, name, level)
                key = (s[0], s[1], s[2])
            res[key] = r
            print(
                f"done {key} {r.get('seconds', float('nan')):.0f}s {r.get('error', '')}", flush=True
            )
    sys.stderr.flush()
    stderr = log_path.read_text(errors="replace") if log_path.exists() else ""
    n_blas = len([ln for ln in stderr.splitlines() if "On entry to" in ln])
    table, extra, crit = evaluate(res, cites, n_blas)
    runtime = time.perf_counter() - t_start
    seen = {(u["file"], u["key"]) for u in cites.used}
    results = {
        "config": {
            "designs": list(DESIGNS),
            "theta": {n: {"value": thetas[n],
                          "source": "zeros" if n == "initial" else {
                              "file": RES.format(THETA_SRC[n][0]),
                              "key": "/".join(THETA_SRC[n][1])}} for n in DESIGNS},  # fmt: skip
            "tolerances": TOL, "dd_onset_J_A_cm2": J_ONSET, "n_dense": N_DENSE,
            "dd_grid_V": [M20.DD_GRID_V[0], M20.DD_GRID_V[-1]], "v_ceiling_V": V_CEIL,
            "barrier_level_VT": M20.BARRIER_LEVEL_VT,
            "T3_interpretation": (
                "gated: no job error, every dense BV resolved with failed_evals = 0, DD "
                "equilibrium converged, stage-16 solves converged, no 'On entry to' line. "
                "Reported, not gated: rejected step-halving trials, unbracketed dense lines "
                "(above 2x the first crossing, cannot change the minimum), and the DD's first "
                "non-converged target after its onset (a T2 'lower bound only' matter if before)"),
            "new_solves_vs_recorded": (
                "AB1 all new (old model, stage-13 device); AB2, AB5 recorded only; AB3 dense BV "
                "recorded (table) and recomputed (cross-check), V_pt,es and DD new; AB4 level 1 "
                "NEW (stage 16 recorded levels 0, 2 only), level 0 re-run as reproduction check"),
        },
        "table": table, **{k: v for k, v in extra.items()},
        "citations": [{"file": f, "key": k} for f, k in sorted(seen)],
        "jobs": {"|".join(map(str, k)): v for k, v in res.items()},
        "on_entry_to_lines": n_blas, "runtime_seconds": runtime, "peak_rss_MiB": rss_mib(),
    }  # fmt: skip
    figdir = (REPO / "figures") if out == DEFAULT_OUT.resolve() else out / "figures"
    label = f"stage 23 | commit {str(man.get('git_commit'))[:8]} | {man['timestamp_utc'][:19]}Z"
    try:
        plot_budget(figdir, table, label)
    except Exception as exc:
        results["figures_error"] = f"{type(exc).__name__}: {exc}"
    (out / "results.json").write_text(json.dumps(results, indent=1, default=float) + "\n")
    summary = {"criteria": crit, "gate_passed": bool(all(crit.values())), "dev": bool(args.dev),
               "table": table, "T1": extra["T1"], "T2": extra["T2"], "T3": extra["T3"],
               "runtime_seconds": runtime}  # fmt: skip
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    print(json.dumps(crit, indent=2))
    for r in table:
        af = r["artifact_fraction"]
        print(
            f"{r['id']:<18} claimed {r['claimed']:.6g} real {r['real']:.6g} {r['unit']} "
            f"artifact {'n/a' if af is None else format(af, '+.4f')}"
        )
    print(f"runtime {runtime:.1f} s")


if __name__ == "__main__":
    main()
