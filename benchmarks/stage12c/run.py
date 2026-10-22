"""Stage-12c benchmark (docs/benchmarks/stage12c_offstate_fix.md): C1-C6 (all gating).

Corrected zero-current off-state model (widegap.offstate): phi_p = 0, phi_n = V_D w with w the
discrete harmonic extension of (0 on the n+ source box, 1 on the other n-type nodes).

Writes <out>/{results,summary,manifest}.json (default results/stage12c).

Choices fixed BEFORE the run (recorded in results.json):

* device for C1-C3, C6: the stage-13b initial rounded device (r = 0.1 um, theta = 0), mesh levels
  h, h/2, h/4 (C1) and h (C2, C3, C6).  C4 uses the sharp (unrounded) stage-12 12b device, level 0.
* C1: peak |E| over the blocking edges (``SiCDMOS.peak_field``) of the electrostatic state at
  V_D = 400 V on h, h/2, h/4.  Observed order p = log2(|e1/e2|) with e1 = v1 - v0, e2 = v2 - v1;
  Richardson v2 + e2/(2^p - 1).  A non-monotone sequence (e1 e2 < 0) is a FAIL.  If both |e1| and
  |e2| are <= 1e-6 |v2| (round-off/converged floor) the order is ill-posed and C1 passes on the
  floor (same rule as stage 12 S2').
* C2: fresh evaluators (same bracket path) for the harmonic and the normalized-graph-distance
  interpolation, level h: (a) peak field at 400 V, (b) planar-cut BV of the ROUNDED device
  (column x = x_bc = 3.0 um, ``SiCAdjoint.bv``), (c) minimum of the 49-line dense set
  (``dense_bv`` with tol 1e-4 V on the near-minimum lines, 1e-2 V on the others, so the
  bracket-refinement tolerance is far below the 1e-6 relative gate).  Relative differences <= 1e-6.
  The planar BV on the sharp device is reported as well.
* C3: n / max(N_A, N_D) <= 1e-6 on every node with 0.01 < w < 0.99 (contact nodes included; the
  worst non-contact node is reported too) for: the 400 V states on h, h/2, h/4, and every
  marched state of the dense march on h up to the initial design's dense-minimum BV (and the
  state at that BV).  Gated on all of them; for the harmonic and the distance w.
* C4: the sharp 12b device, planar cut (x_min = x_bc, Niwa, electron form, level 0), corrected
  model vs the stored stage-12 value in results/stage12b/results.json (R2_BV, level0|niwa|
  planar_peak), relative difference <= 0.5 %.  Diagnostics: reproduction of the stored value by
  the old rule, the 1D flat-quasi-Fermi column diode BV through the same column, and the electron
  density at the junction node in the old state.
* C5: laterally uniform 2D p+n (N_D = 1e16, 12a node set extruded onto 3 columns, same mesh as
  the stored 12a levels 0 and 1) vs the stored 12a Niwa BV (results/stage12b/results.json), <= 1e-6
  relative at each level.  The 1D value is also recomputed (reconstruction check); the old rule
  is reported.
* C6: analytic ``jacobian_bc`` vs central difference of ``residual_bc`` (Richardson over steps
  1e-2 and 5e-3), all 3N columns, per row error max_j |J_ij - D_ij| / max_j |J_ij| <= 1e-6, at the
  400 V state of the level-0 rounded device; the round-off floor of each row is reported.
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

from widegap.adjoint_sic_corner import DenseCuts, dense_bv, solve_at_voltage, state_at
from widegap.fv.avalanche import avalanche_bv
from widegap.manifest import build_manifest
from widegap.offstate import (
    CorrectedCornerAdjoint,
    UniformPN2D,
    corrected_class,
    jacobian_fd_rows,
    n_consistency,
)
from widegap.sic import IonizedReverseDiode, SiC4H, SiCDMOS, TracedRegistry

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[2]
DEFAULT_OUT = REPO / "results" / "stage12c"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


M13B = _load("stage13b_run", REPO / "benchmarks/stage13b/run.py")
DESIGN, TH0, R_UM, TAU_V = M13B.DESIGN, M13B.TH0, M13B.R_UM, M13B.TAU_V
MODES = ("harmonic", "distance")
V_C1 = 400.0  # V (spec C1)
TOL = {"C1_order": 0.8, "C1_rich": 0.02, "C2": 1e-6, "C3": 1e-6, "C4": 5e-3, "C5": 1e-6, "C6": 1e-6}
STORED = REPO / "results" / "stage12b" / "results.json"


def rss_mib():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def make_ev(level, mode):
    sc = SiC4H(TracedRegistry.from_json())
    return CorrectedCornerAdjoint(
        sc, DESIGN, level, r_um=R_UM, cuts_um=(), tau_v=TAU_V, w_mode=mode
    )


def peak_record(ev, X):
    sc = ev.sc
    psi = X[0::3]
    pk = ev.base.peak_field(psi)
    um = lambda v: float(sc.length_cm(v)) * 1e4  # noqa: E731
    return {"E_peak_MV_cm": float(sc.field_phys(pk[0])) * 1e-6, "x_um": um(pk[1]), "y_um": um(pk[2]),
            "dir": pk[3]}  # fmt: skip


# ============================================================ jobs (module level: fork pool)
def job_peak(args):
    """C1/C3: electrostatic state at V_D = 400 V, peak field and the n-consistency check."""
    level, mode, v = args
    t0 = time.perf_counter()
    ev = make_ev(level, mode)
    X = solve_at_voltage(ev, TH0, v)
    out = {"level": level, "mode": mode, "v_volt": v, "n_nodes": int(ev.mesh.n_nodes),
           **peak_record(ev, X), "C3": n_consistency(ev.base, X)}  # fmt: skip
    out["seconds"] = time.perf_counter() - t0
    return out


def job_mode(args):
    """C2/C3 for one interpolation: peak at 400 V, planar-cut BV (rounded device), dense minimum."""
    (mode,) = args
    t0 = time.perf_counter()
    ev = make_ev(0, mode)
    sc = ev.sc
    x400 = solve_at_voltage(ev, TH0, V_C1)
    peak = peak_record(ev, x400)
    ev.reset_caches()
    planar = ev.bv(TH0, grad=False)["BV"]
    xs = M13B.dense_positions(sc, ev)
    ev2 = make_ev(0, mode)
    dc = DenseCuts(ev2, xs)
    d = dense_bv(ev2, TH0, dc, tol_v=1e-4, loose_tol_v=1e-2)
    bv_min = d["bv_min"]
    v_bv = float(sc.volt_s(bv_min))
    c3_states = []
    for v_s, X, _ in d["states"]:
        if v_s <= v_bv * (1.0 + 1e-9):
            c3_states.append((sc.volt_phys(v_s), n_consistency(ev2.base, X)))
    c3_bv = n_consistency(ev2.base, state_at(ev2, d["states"], v_bv))
    worst = max([c for _, c in c3_states] + [c3_bv], key=lambda c: c["worst_ratio"])
    worst_v = ([v for v, c in c3_states if c is worst] + [bv_min])[0]
    beyond = [
        (sc.volt_phys(v_s), n_consistency(ev2.base, X)["worst_ratio"])
        for v_s, X, _ in d["states"] if v_s > v_bv * (1.0 + 1e-9)
    ]  # fmt: skip
    return {
        "mode": mode, "peak_400V": peak, "planar_bv_rounded_V": planar,
        "dense_bv_min_V": bv_min, "dense_x_argmin_um": float(dc.x_um[d["argmin"]]),
        "dense_bv": d["bv"].tolist(), "dense_x_um": xs.tolist(),
        "dense_resolved_all": d["resolved_all"], "dense_n_unbracketed": d["n_unbracketed"],
        "C3_n_states_checked": len(c3_states) + 1, "C3_worst": worst, "C3_worst_at_V": worst_v,
        "C3_worst_ratio_all_states_up_to_BV": worst["worst_ratio"],
        "C3_worst_ratio_beyond_BV_diagnostic": max([r for _, r in beyond], default=None),
        "seconds": time.perf_counter() - t0,
    }  # fmt: skip


def job_sharp(args):
    """C4: sharp 12b device, planar cut; corrected (both interpolations) vs old rule vs 1D column."""
    t0 = time.perf_counter()
    sc = SiC4H(TracedRegistry.from_json())
    ch = sc.niwa()
    out = {}
    for name, cls in (("old", SiCDMOS), ("harmonic", corrected_class(SiCDMOS, "harmonic")),
                      ("distance", corrected_class(SiCDMOS, "distance"))):  # fmt: skip
        d = cls(sc, DESIGN, 0)
        eq = d.equilibrium()
        r = d.offstate_bv(ch, sc.volt_s(20.0), sc.volt_s(1500.0), eq.X, "n", x_min=d.prm.x_bc)
        out[name] = {"BV_V": float(sc.volt_phys(r["bv"])), "x_col_um": float(sc.length_cm(r["x_col"])) * 1e4,
                     "consistent": r["consistent"], "peak_E_MV_cm": float(sc.field_phys(r["peak_E"])) * 1e-6}  # fmt: skip
        nodes = d.column_nodes(d.prm.x_bc)[1]
        psi, pn = r["X"][0::3][nodes], r["X"][1::3][nodes]
        with np.errstate(over="ignore"):
            ratio = np.exp(d.model.ld + psi - pn) / np.maximum(d.nd[nodes], d.na[nodes])
        k = int(np.argmax(ratio))
        out[name]["column_max_n_over_N"] = float(ratio[k])
        out[name]["column_max_n_over_N_at_y_um"] = (
            float(sc.length_cm(d.dev.mesh.node_xy[nodes[k], 1])) * 1e4
        )
        if name == "harmonic":
            out[name]["C3_at_BV"] = n_consistency(d, r["X"])
            col = nodes
    # 1D flat-quasi-Fermi column diode through the same column (no lateral fringing)
    depth = -d.dev.mesh.node_xy[col, 1]
    dio = IonizedReverseDiode(depth, d.nd[col], d.na[col], sc.delta, d.donors, d.acceptors)
    b = avalanche_bv(dio, ch, "n", v_step=float(sc.volt_s(20.0)), v_max=float(sc.volt_s(1500.0)),
                     xtol=1e-12)  # fmt: skip
    out["column_1d_diode_BV_V"] = float(sc.volt_phys(b.bv))
    out["seconds"] = time.perf_counter() - t0
    return out


def job_c5(args):
    """C5: laterally uniform 2D p+n, N_D = 1e16, level `level`."""
    level, nd_cm = args
    t0 = time.perf_counter()
    sc = SiC4H(TracedRegistry.from_json())
    ch = sc.niwa()
    pn = UniformPN2D(sc, nd_cm, level=level)
    d = pn.diode
    bv1 = sc.volt_phys(avalanche_bv(d, ch, "n", v_step=pn.v_est / 3.0, v_max=6.0 * pn.v_est,
                                    xtol=1e-12).bv)  # fmt: skip
    out = {"level": level, "N_D_cm3": nd_cm, "n_nodes_2d": int(pn.mesh.n_nodes),
           "n_nodes_1d": int(len(d.x)), "bv_1d_recomputed_V": float(bv1)}  # fmt: skip
    for key, rule, mode in (("new_harmonic", "new", "harmonic"), ("new_distance", "new", "distance"),
                            ("old_rule", "old", "harmonic")):  # fmt: skip
        try:
            out[key] = float(sc.volt_phys(pn.bv(ch, rule, mode)))
        except (RuntimeError, ValueError) as exc:
            out[key] = f"{type(exc).__name__}: {exc}"
    out["w_min_p_side"] = float(np.min(pn.w("harmonic")[pn.reg.semi_node]))
    out["seconds"] = time.perf_counter() - t0
    return out


def job_c6(args):
    """C6: analytic Jacobian vs FD at the 400 V state; also adjoint vs tangent dI/dV_D."""
    (step,) = args
    t0 = time.perf_counter()
    ev = make_ev(0, "harmonic")
    X = solve_at_voltage(ev, TH0, V_C1)
    rows, vals = ev._off_bc(float(ev.sc.volt_s(V_C1)))
    fd = jacobian_fd_rows(ev.mod, X, rows, vals, step=step)
    sel = np.ones(fd["n_rows"], dtype=bool)
    sel[rows] = False  # free rows (Dirichlet rows are identity: trivially exact)
    o = ev.bv(TH0, grad=True)
    wr = fd["worst_row"]
    return {
        "step": step, "n_rows": fd["n_rows"], "n_free_rows": int(sel.sum()),
        "max_err": float(fd["err"].max()), "max_err_free_rows": float(fd["err"][sel].max()),
        "max_err_single_step": float(fd["err_single_step"].max()),
        "max_floor": float(fd["floor"][sel].max()), "worst_row": wr,
        "worst_row_var": ["psi", "phi_n", "phi_p"][wr % 3], "worst_row_node": wr // 3,
        "worst_row_err": float(fd["err"][wr]), "worst_row_floor": float(fd["floor"][wr]),
        "median_err_free": float(np.median(fd["err"][sel])),
        "bv_planar_V": o["BV"], "di_dv_adjoint": o["di_dv_adjoint"], "di_dv_tangent": o["di_dv_tangent"],
        "di_dv_rel_diff": abs(o["di_dv_adjoint"] - o["di_dv_tangent"]) / abs(o["di_dv_tangent"]),
        "seconds": time.perf_counter() - t0,
    }  # fmt: skip


# ============================================================ analysis
def richardson3(vals, floor_rel=1e-6):
    e1, e2 = vals[1] - vals[0], vals[2] - vals[1]
    floor = floor_rel * abs(vals[2])
    if abs(e1) <= floor and abs(e2) <= floor:
        return {"order": None, "richardson": float(vals[2]), "rel_dist": float(abs(e2) / abs(vals[2])),
                "monotone": True, "floor_case": True, "pass": True}  # fmt: skip
    if e1 == 0.0 or e2 == 0.0 or e1 * e2 < 0.0:
        return {"order": None, "richardson": None, "rel_dist": None, "monotone": False,
                "floor_case": False, "pass": False}  # fmt: skip
    p = float(np.log2(abs(e1 / e2)))
    rich = vals[2] + e2 / (2.0**p - 1.0)
    rel = abs(vals[2] - rich) / abs(rich)
    return {"order": p, "richardson": float(rich), "rel_dist": float(rel), "monotone": True,
            "floor_case": False, "pass": bool(p >= TOL["C1_order"] and rel <= TOL["C1_rich"])}  # fmt: skip


def rel(a, b):
    return float(abs(a - b) / abs(b))


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--workers", type=int, default=min(10, os.cpu_count() or 1))
    args = ap.parse_args(argv)
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    man = build_manifest(REPO, sys.argv, [])
    (out / "manifest.json").write_text(json.dumps(man, indent=2) + "\n")
    t_start = time.perf_counter()
    stored = json.loads(STORED.read_text())
    bv12b = stored["12b"]["R2_BV"]["level0|niwa|planar_peak"]["BV_V"]
    s12a = [e for e in stored["12a"] if e["N_D_cm3"] == 1e16][0]["sets"]["niwa"]["levels"]
    ref12a = {lev: s12a[lev]["bv_n_V"] for lev in (0, 1)}

    ctx = mp.get_context("fork")
    res: dict = {"config": {
        "rounding_length_um": R_UM, "tau_V": TAU_V, "v_C1_volt": V_C1, "tolerances": TOL,
        "model": "phi_p = 0, phi_n = V_D w; w harmonic (gw-weighted Laplace on the semiconductor "
                 "graph; w=0 on n-type nodes of the (enlarged) source box, w=1 on other n-type "
                 "nodes) | normalized graph distance (C2)",
        "pinned_sets": "n-type source-box nodes (0) and other n-type nodes (1); p-type contact "
                       "nodes are ordinary Laplace nodes (ohmic phi_n overrides)",
        "stored_12b_planar_BV_V": bv12b, "stored_12a_niwa_1e16_V": ref12a,
        "workers": args.workers,
    }}  # fmt: skip
    jobs = {}
    t0 = time.perf_counter()
    with ProcessPoolExecutor(args.workers, mp_context=ctx) as pool:
        futs = {}
        for lev in (0, 1, 2):
            futs[pool.submit(job_peak, (lev, "harmonic", V_C1))] = ("peak", ("harmonic", lev))
        for lev in (0, 1):
            futs[pool.submit(job_peak, (lev, "distance", V_C1))] = ("peak", ("distance", lev))
        for mode in MODES:
            futs[pool.submit(job_mode, (mode,))] = ("mode", mode)
        futs[pool.submit(job_sharp, ())] = ("sharp", None)
        for lev in (0, 1):
            futs[pool.submit(job_c5, (lev, 1e16))] = ("c5", lev)
        futs[pool.submit(job_c6, (1e-2,))] = ("c6", None)
        for f in as_completed(futs):
            kind, key = futs[f]
            r = f.result()
            jobs[(kind, key)] = r
            print(f"done: {kind} {key} {r['seconds']:.0f}s", flush=True)
    t_par = time.perf_counter() - t0

    # ---- C1 ---------------------------------------------------------------------------------
    pk = [jobs[("peak", ("harmonic", lev))] for lev in (0, 1, 2)]
    vals = [p["E_peak_MV_cm"] for p in pk]
    c1 = richardson3(vals)
    c1.update({"E_peak_MV_cm_h_h2_h4": vals, "locations_um": [[p["x_um"], p["y_um"]] for p in pk],
               "n_nodes": [p["n_nodes"] for p in pk], "v_volt": V_C1})  # fmt: skip

    # ---- C2 ---------------------------------------------------------------------------------
    jh, jd = jobs[("mode", "harmonic")], jobs[("mode", "distance")]
    sh = jobs[("sharp", None)]
    pd0 = jobs[("peak", ("distance", 0))]
    c2 = {
        "peak_rel": rel(jd["peak_400V"]["E_peak_MV_cm"], jh["peak_400V"]["E_peak_MV_cm"]),
        "planar_bv_rounded_rel": rel(jd["planar_bv_rounded_V"], jh["planar_bv_rounded_V"]),
        "dense_min_rel": rel(jd["dense_bv_min_V"], jh["dense_bv_min_V"]),
        "values": {"peak_MV_cm": [jh["peak_400V"]["E_peak_MV_cm"], jd["peak_400V"]["E_peak_MV_cm"]],
                   "planar_bv_rounded_V": [jh["planar_bv_rounded_V"], jd["planar_bv_rounded_V"]],
                   "dense_min_V": [jh["dense_bv_min_V"], jd["dense_bv_min_V"]],
                   "dense_x_argmin_um": [jh["dense_x_argmin_um"], jd["dense_x_argmin_um"]]},
        "dense_profile_max_abs_diff_V": float(np.nanmax(np.abs(np.array(jh["dense_bv"]) - np.array(jd["dense_bv"])))),
        "peak_independent_job_distance_rel": rel(pd0["E_peak_MV_cm"], jh["peak_400V"]["E_peak_MV_cm"]),
        "sharp_planar_bv_rel": rel(sh["distance"]["BV_V"], sh["harmonic"]["BV_V"]),
    }  # fmt: skip
    c2["pass"] = bool(
        max(c2["peak_rel"], c2["planar_bv_rounded_rel"], c2["dense_min_rel"]) <= TOL["C2"]
    )

    # ---- C3 ---------------------------------------------------------------------------------
    c3_items = {}
    for lev in (0, 1, 2):
        c3_items[f"harmonic_400V_L{lev}"] = jobs[("peak", ("harmonic", lev))]["C3"]
    for lev in (0, 1):
        c3_items[f"distance_400V_L{lev}"] = jobs[("peak", ("distance", lev))]["C3"]
    c3_items["harmonic_march_to_BV_L0"] = jh["C3_worst"]
    c3_items["distance_march_to_BV_L0"] = jd["C3_worst"]
    c3_items["sharp_planar_BV_state"] = sh["harmonic"]["C3_at_BV"]
    worst_key = max(c3_items, key=lambda k: c3_items[k]["worst_ratio"])
    c3 = {"pass": bool(all(v["ok"] for v in c3_items.values())),
          "worst_ratio": c3_items[worst_key]["worst_ratio"], "worst_case": worst_key,
          "worst_node": c3_items[worst_key], "tol": TOL["C3"], "cases": c3_items,
          "worst_ratio_non_contact": max(v["free_nodes"]["worst_ratio"] for v in c3_items.values()),
          "n_states_in_march_checks": [jh["C3_n_states_checked"], jd["C3_n_states_checked"]],
          "march_worst_at_V": [jh["C3_worst_at_V"], jd["C3_worst_at_V"]],
          "beyond_BV_diagnostic": [jh["C3_worst_ratio_beyond_BV_diagnostic"], jd["C3_worst_ratio_beyond_BV_diagnostic"]]}  # fmt: skip

    # ---- C4 ---------------------------------------------------------------------------------
    c4 = {
        "bv_new_harmonic_V": sh["harmonic"]["BV_V"], "bv_new_distance_V": sh["distance"]["BV_V"],
        "bv_stored_V": bv12b, "bv_old_rule_recomputed_V": sh["old"]["BV_V"],
        "rel_new_vs_stored": rel(sh["harmonic"]["BV_V"], bv12b),
        "old_reproduces_stored_rel": rel(sh["old"]["BV_V"], bv12b),
        "column_1d_flat_qf_diode_BV_V": sh["column_1d_diode_BV_V"],
        "rel_old_vs_1d": rel(sh["old"]["BV_V"], sh["column_1d_diode_BV_V"]),
        "rel_new_vs_1d": rel(sh["harmonic"]["BV_V"], sh["column_1d_diode_BV_V"]),
        "column_max_n_over_N_old": sh["old"]["column_max_n_over_N"],
        "column_max_n_over_N_old_at_y_um": sh["old"]["column_max_n_over_N_at_y_um"],
        "column_max_n_over_N_new": sh["harmonic"]["column_max_n_over_N"],
        "peak_E_MV_cm_old_new": [sh["old"]["peak_E_MV_cm"], sh["harmonic"]["peak_E_MV_cm"]],
    }  # fmt: skip
    c4["pass"] = bool(c4["rel_new_vs_stored"] <= TOL["C4"])

    # ---- C5 ---------------------------------------------------------------------------------
    c5 = {"levels": {}}
    for lev in (0, 1):
        j = jobs[("c5", lev)]
        r = {**j, "stored_12a_V": ref12a[lev], "rel_new_harmonic": rel(j["new_harmonic"], ref12a[lev]),
             "rel_new_distance": rel(j["new_distance"], ref12a[lev]),
             "rel_1d_recomputed": rel(j["bv_1d_recomputed_V"], ref12a[lev])}  # fmt: skip
        if isinstance(j["old_rule"], float):
            r["rel_old_rule"] = rel(j["old_rule"], ref12a[lev])
        c5["levels"][str(lev)] = r
    c5["pass"] = bool(all(v["rel_new_harmonic"] <= TOL["C5"] for v in c5["levels"].values()))

    # ---- C6 ---------------------------------------------------------------------------------
    j6 = jobs[("c6", None)]
    c6 = {**j6, "pass": bool(j6["max_err"] <= TOL["C6"])}

    crit = {"C1": bool(c1["pass"]), "C2": c2["pass"], "C3": c3["pass"], "C4": c4["pass"],
            "C5": c5["pass"], "C6": c6["pass"]}  # fmt: skip
    res.update({"C1": c1, "C2": c2, "C3": c3, "C4": c4, "C5": c5, "C6": c6,
                "jobs_raw": {f"{k[0]}|{k[1]}": v for k, v in jobs.items()},
                "runtime_seconds": time.perf_counter() - t_start, "parallel_seconds": t_par,
                "peak_rss_MiB": rss_mib()})  # fmt: skip
    (out / "results.json").write_text(json.dumps(res, indent=2, default=float) + "\n")
    summary = {"criteria": crit, "gate_passed": all(crit.values()),
               "values": {k: res[k] for k in ("C1", "C2", "C3", "C4", "C5", "C6")},
               "runtime_seconds": res["runtime_seconds"]}  # fmt: skip
    # keep the summary light: drop per-case C3 nodes and per-job raw data
    summary["values"]["C3"] = {k: v for k, v in c3.items() if k != "cases"}
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    print(json.dumps(crit, indent=2))
    print(f"runtime {res['runtime_seconds']:.1f} s")


if __name__ == "__main__":
    main()
