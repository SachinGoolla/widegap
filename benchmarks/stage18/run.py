"""Stage-18 benchmark (docs/benchmarks/stage18_accumulation_bracket.md): A1-A5 (gating), R1-R2.

Bracket of the JFET accumulation-surface mobility in the stage-16 channel model:
(U) bulk Hall mobility (stage-16 default) and (L) the constant Zhao (2007) accumulation mobility
(registry entry ``mu_acc_n_zhao2007``, read through the registry) on lateral electron edges within
d_ch = 20 nm of the interface in the n-type surface columns under the gate.  The bracket reported
is [R_on(U), R_on(L)]; no best estimate is claimed.

Protocol per (design, level): stage-16 ``ron`` protocol (direct bias march to V_GS = 15 V, V_DS =
0.1 V, then the Picard channel loop, ``ChannelMap`` built on that design/level); R_on,sp =
V_DS x_cell / I_D.  Designs: 12b, 13d initial, 13d final, stage-17 final (results/stage17) and, if
present when the run is made, the 17b final (results/stage17b; recorded as absent otherwise).

* A1 : check_entry ok for the Zhao entry (incl. verified_by / verified_value_match).
* A2 : (U) reproduces the stage-16 channel R_on (12b, level 1, d_ch = 20 nm) to <= 1e-12 relative.
* A3 : every Picard / Newton solve converged; captured fd-2 text has no "On entry to" line.
* A4 : R_on(L) 12b level 1 vs level 2 <= 2 %.
* A5 : R_on(L) 12b level 1, d_ch 20 vs 40 nm <= 1 %.
* R1 : the bracket on levels 1 and 2 for every design.
* R2 : FD dR_on/dtheta_J (Richardson eps = 1e-3, 5e-4) at the 13d initial design, level 1, under
       (L) and (U) (and the bulk-mobility model for reference).  The channel map is frozen at the
       base design (as stage 17 freezes it at theta = 0) so the FD is not polluted by the discrete
       edge set changing.

Writes <out>/{results,summary,manifest}.json and stderr.log (default results/stage18) and figures
stage18_*.png (repo figures/ for the default output, else <out>/figures).  ``--dev`` marks the run
DEV (refused when --out lies under results/).
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

from widegap.adjoint_sic import richardson_fd
from widegap.channel_mobility import (
    ACC_NAME,
    ChannelMap,
    acc_mobility_from_registry,
    picard_solve,
)
from widegap.manifest import build_manifest
from widegap.materials.registry import ProvenanceError, Registry, check_entry
from widegap.sic import REGISTRY_PATH, TracedRegistry

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[2]
DEFAULT_OUT = REPO / "results" / "stage18"
R16 = REPO / "results" / "stage16" / "results.json"
R17 = REPO / "results" / "stage17" / "results.json"
R17B = REPO / "results" / "stage17b" / "results.json"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


S16 = _load("stage16_run", REPO / "benchmarks/stage16/run.py")  # fixed names: fork pickling
S12, S15F, TH0 = S16.S12, S16.S15F, S16.S15F.TH0
V_ON, V_LIN = S16.V_ON, S16.V_LIN

D_CH_NM = 20.0
A2_TOL, A4_TOL, A5_TOL = 1e-12, 0.02, 0.01
LEVELS = (1, 2)
EPS_FD = (1e-3, 5e-4)  # Richardson base step and its half
K_JFET = 4  # theta_J


def rss_mib():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


# ============================================================ registry / designs
def registry_entries():
    return json.loads(REGISTRY_PATH.read_text())["parameters"]


def a1_check():
    ent = {e["name"]: e for e in registry_entries()}.get(ACC_NAME)
    out = {"name": ACC_NAME, "in_sic4h_json": ent is not None}
    if ent is None:
        out.update(check_entry="MISSING", ok=False)
        return out
    try:
        check_entry(ent)
        out["check_entry"] = "ok"
    except (ProvenanceError, TypeError) as exc:
        out["check_entry"] = f"REFUSED: {exc}"
    out.update(verified_by=ent.get("verified_by"), verified_value_match=ent.get("verified_value_match"),
               value=ent.get("value"), unit=ent.get("unit"), doi=ent.get("doi"),
               locator=ent.get("locator"))  # fmt: skip
    out["ok"] = bool(out["check_entry"] == "ok" and ent.get("verified_by")
                     and ent.get("verified_value_match") is True)  # fmt: skip
    return out


def designs(s17b_path: Path):
    """name -> theta (None for 12b) for every design of R1; 17b recorded as absent if missing."""
    d = {"12b": None, "13d_initial": None, "13d_final": None}
    th17 = np.asarray(json.loads(R17.read_text())["final"]["theta"], float)
    d["s17_final"] = th17.tolist()
    info = {"s17_final": {"source": str(R17.relative_to(REPO)), "present": True}}
    if s17b_path.exists():
        r = json.loads(s17b_path.read_text())
        th = r.get("final", {}).get("theta")
        complete = "runtime_seconds" in r
        info["s17b_final"] = {"source": str(s17b_path), "present": th is not None,
                              "complete_results_file": bool(complete), "dev": r.get("dev")}  # fmt: skip
        if th is not None:
            d["s17b_final"] = [float(t) for t in th]
    else:
        info["s17b_final"] = {"source": str(s17b_path), "present": False}
    return d, info


# ============================================================ device
def make_dev(case, level, theta, entries):
    """Stage-16 device; ``theta`` (list) overrides the design (13d family only)."""
    base = {"12b": "12b", "13d_initial": "13d_initial", "13d_final": "13d_final"}.get(case)
    sc, P, dev = S16._setup(base or "13d_initial", level, entries)
    if theta is not None:
        dev.ev.set_theta(np.asarray(theta, float))
        dev.nd, dev.na = dev.ev.nd, dev.ev.na
    return sc, P, dev


def ron_solve(dev, sc, P, cmap, th=None, mode="channel"):
    """Direct march + (channel) Picard at theta ``th`` (None: as built).  Returns (R, record)."""
    if th is not None:
        dev.ev.set_theta(np.asarray(th, float))
        dev.nd, dev.na = dev.ev.nd, dev.ev.na
        cmap.mu_bulk = np.asarray(dev.model.mu_edge, dtype=float).copy()
    dev.model.mu_edge = cmap.mu_bulk.copy()
    eq = dev.equilibrium()
    sol = S16._bias_to(dev, sc, eq.X, V_ON, V_LIN)
    ok_bulk = bool(eq.report.converged and sol.report.converged)
    r_bulk = float(V_LIN * dev.x_cm / S16._id_a_per_cm(dev, sol))
    if mode == "bulk":
        return r_bulk, {"R_bulk": r_bulk, "converged": ok_bulk}
    vds_s = float(sc.volt_s(V_LIN))
    r = picard_solve(cmap, P, sol.X, dev.bc(float(sc.volt_s(V_ON)), 0.0, vds_s))
    rec = {"picard": S16._picard_record(r), "R_bulk": r_bulk, "bulk_converged": ok_bulk,
           "converged": bool(ok_bulk and r["converged"])}  # fmt: skip
    if r["converged"]:
        rec["R_on_sp_ohm_cm2"] = float(V_LIN * dev.x_cm / S16._id_a_per_cm(dev, r["sol"]))
        rec["counts"] = r["counts"]
    dev.model.mu_edge = cmap.mu_bulk.copy()
    return rec.get("R_on_sp_ohm_cm2", float("nan")), rec


def acc_info(cmap, sc):
    if cmap.mu_acc_cm2 is None:
        return {"n_acc_columns": 0, "n_acc_edges": 0}
    cols = np.flatnonzero(cmap.col_acc)
    x_um = (
        [float(sc.length_cm(cmap.mesh.x[i])) * 1e4 for i in (cols[0], cols[-1])]
        if len(cols)
        else []
    )
    return {"n_acc_columns": int(len(cols)), "n_acc_edges": int(len(cmap.acc_edges)),
            "acc_x_range_um": x_um, "mu_acc_cm2": cmap.mu_acc_cm2}  # fmt: skip


# ============================================================ jobs
def job_bracket(args):
    """[R_on(U), R_on(L)] for one design on one level (+ d_ch variants listed in ``d_ch_nm``)."""
    case, level, theta, d_ch_nm, entries = args
    t0 = time.perf_counter()
    rec = {"kind": "bracket", "case": case, "level": level, "error": None}
    try:
        sc, P, dev = make_dev(case, level, theta, entries)
        mu_acc = acc_mobility_from_registry(TracedRegistry(Registry(entries)))
        rec["nodes"] = int(dev.model.mesh.n_nodes)
        rec["mu_acc_cm2Vs"] = mu_acc
        gx = dev.prm.gate_x1
        out = {}
        for tag, dnm, acc in (("U", D_CH_NM, None), ("L", D_CH_NM, mu_acc),
                              *[(f"L_d{d:g}", d, mu_acc) for d in d_ch_nm]):  # fmt: skip
            cmap = ChannelMap(sc, dev.model, dev.nd, dev.na, dnm * 1e-3, gate_x_max=gx,
                              mu_acc_cm2=acc)  # fmt: skip
            _, r = ron_solve(dev, sc, P, cmap)
            r.update(acc_info(cmap, sc))
            r["n_channel_edges"] = int(len(cmap.edges))
            out[tag] = r
        rec["runs"] = out
        rec["R_on_bulk_ohm_cm2"] = out["U"]["R_bulk"]
    except Exception as exc:  # recorded, never hidden
        rec["error"] = f"{type(exc).__name__}: {exc}"
    rec["seconds"] = time.perf_counter() - t0
    rec["rss_mib"] = rss_mib()
    return rec


def job_r2(args):
    """Richardson FD of R_on w.r.t. theta_J at the 13d initial design (level 1) for one model."""
    model, level, entries = args  # model in {"bulk", "U", "L"}
    t0 = time.perf_counter()
    rec = {"kind": "r2", "model": model, "level": level, "error": None}
    try:
        sc, P, dev = make_dev("13d_initial", level, None, entries)
        mu_acc = acc_mobility_from_registry(TracedRegistry(Registry(entries)))
        cmap = ChannelMap(sc, dev.model, dev.nd, dev.na, D_CH_NM * 1e-3,
                          gate_x_max=dev.prm.gate_x1,
                          mu_acc_cm2=mu_acc if model == "L" else None)  # fmt: skip
        th0 = np.asarray(TH0, float)
        flags = []

        def f(th):
            r, d = ron_solve(dev, sc, P, cmap, th, "bulk" if model == "bulk" else "channel")
            flags.append(bool(d["converged"]))
            return d["R_bulk"] if model == "bulk" else r

        r0 = f(th0)
        u = np.zeros(5)
        u[K_JFET] = 1.0
        rich, f1, f2 = richardson_fd(f, th0, u, EPS_FD[0])
        rec.update({"R0": r0, "dR_dthetaJ": float(rich), "fd_eps_1e-3": float(f1),
                    "fd_eps_5e-4": float(f2), "spread": float(abs(f1 - f2)),
                    "theta0": th0.tolist(), "n_solves": len(flags), "all_converged": all(flags),
                    **acc_info(cmap, sc)})  # fmt: skip
    except Exception as exc:
        rec["error"] = f"{type(exc).__name__}: {exc}"
    rec["seconds"] = time.perf_counter() - t0
    rec["rss_mib"] = rss_mib()
    return rec


# ============================================================ figures
def figures(res, figdir: Path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figdir.mkdir(parents=True, exist_ok=True)
    rows = [r for r in res["R1_table"] if r["R_on_U_ohm_cm2"] and r["R_on_L_ohm_cm2"]]
    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    xs = np.arange(len(rows))
    ax.bar(xs, [r["R_on_U_ohm_cm2"] * 1e4 for r in rows], 0.4, label="(U) bulk Hall (optimistic)")
    ax.bar(
        xs + 0.4,
        [r["R_on_L_ohm_cm2"] * 1e4 for r in rows],
        0.4,
        label="(L) Zhao 17.5 (pessimistic)",
    )
    ax.set_xticks(xs + 0.2, [f"{r['case']}\nlevel {r['level']}" for r in rows], fontsize=7)
    ax.set_ylabel(r"$R_{on,sp}$ [$10^{-4}\ \Omega\,$cm$^2$]")
    ax.set_title("Stage 18: accumulation-mobility bracket [R_on(U), R_on(L)]")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(figdir / "stage18_bracket.png", dpi=130)
    plt.close(fig)
    r2 = res["R2"]["per_model"]
    if all(r2.get(m) for m in ("bulk", "U", "L")):
        fig, ax = plt.subplots(figsize=(5.5, 4.0))
        vals = [r2[m]["dR_dthetaJ"] * 1e4 for m in ("bulk", "U", "L")]
        ax.bar(["bulk", "(U)", "(L)"], vals, color=["0.6", "C0", "C1"])
        ax.set_ylabel(r"$dR_{on}/d\theta_J$ [$10^{-4}\ \Omega\,$cm$^2$]")
        ax.set_title("R2: JFET-level sensitivity, 13d initial, h/2")
        ax.grid(alpha=0.3, axis="y")
        fig.tight_layout()
        fig.savefig(figdir / "stage18_r2.png", dpi=130)
        plt.close(fig)


# ============================================================ main
def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--workers", type=int, default=min(4, os.cpu_count() or 1))
    ap.add_argument("--stage17b-results", type=Path, default=R17B)
    ap.add_argument("--dev", action="store_true", help="DEV run (refused under results/)")
    args = ap.parse_args(argv)
    out = args.out.resolve()
    if args.dev and (REPO / "results") in (out, *out.parents):
        raise SystemExit("refusing --dev with --out under results/")
    out.mkdir(parents=True, exist_ok=True)
    log_path = out / "stderr.log"  # capture C-level stderr (OpenBLAS writes to fd 2)
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
    a1 = a1_check()
    if not a1["ok"]:
        (out / "summary.json").write_text(json.dumps({"criteria": {"A1": False},
                                                      "gate_passed": False}, indent=2) + "\n")  # fmt: skip
        (out / "results.json").write_text(json.dumps({"A1": a1}, indent=2) + "\n")
        print("A1 failed: Zhao entry unusable; nothing run")
        return
    entries = registry_entries()
    dsg, dinfo = designs(args.stage17b_results)
    jobs = []
    for lv in LEVELS:
        for c, th in dsg.items():
            extra = (40.0,) if (c == "12b" and lv == 1) else ()
            jobs.append((job_bracket, ("bracket", c, lv), (c, lv, th, extra, entries)))
    for m in ("bulk", "U", "L"):
        jobs.append((job_r2, ("r2", m, 1), (m, 1, entries)))
    jobs.sort(key=lambda j: -j[1][2] if j[1][0] == "bracket" else -1)  # level 2 first
    runs: dict = {}
    ctx = mp.get_context("fork")
    with ProcessPoolExecutor(args.workers, mp_context=ctx) as pool:
        futs = {pool.submit(fn, a): key for fn, key, a in jobs}
        for f in as_completed(futs):
            r = f.result()
            runs["|".join(str(k) for k in futs[f])] = r
            print(
                f"done {futs[f]} {r['seconds']:.0f}s rss={r['rss_mib']:.0f}MiB err={r['error']}",
                flush=True,
            )
    # ---- criteria
    crit, vals = {"A1": bool(a1["ok"])}, {}

    def ron(case, lv, tag):
        r = runs.get(f"bracket|{case}|{lv}")
        if not r or r["error"]:
            return None
        return r["runs"].get(tag, {}).get("R_on_sp_ohm_cm2")

    s16 = json.loads(R16.read_text())["runs"]["ron|12b|1"]["channel"]["20"]["R_on_sp_ohm_cm2"]
    u1 = ron("12b", 1, "U")
    vals["A2"] = {"R_on_U_this_run": u1, "R_on_stage16": s16,
                  "rel_diff": None if u1 is None else abs(u1 / s16 - 1.0)}  # fmt: skip
    crit["A2"] = bool(u1 is not None and abs(u1 / s16 - 1.0) <= A2_TOL)
    failed = []
    for k, r in runs.items():
        if r["error"]:
            failed.append(k)
        elif r["kind"] == "bracket":
            failed += [f"{k}@{t}" for t, e in r["runs"].items() if not e["converged"]]
        elif not r["all_converged"]:
            failed.append(k)
    sys.stderr.flush()
    stderr_text = log_path.read_text(errors="replace") if log_path.exists() else ""
    entry_lines = [ln for ln in stderr_text.splitlines() if "On entry to" in ln]
    vals["A3"] = {"failed": failed, "on_entry_to_lines": len(entry_lines),
                  "stderr_bytes": len(stderr_text)}  # fmt: skip
    crit["A3"] = bool(not failed and not entry_lines and runs)
    l1, l2 = ron("12b", 1, "L"), ron("12b", 2, "L")
    vals["A4"] = {"R_on_L_level1": l1, "R_on_L_level2": l2,
                  "rel_diff": None if not (l1 and l2) else abs(l2 / l1 - 1.0)}  # fmt: skip
    crit["A4"] = bool(l1 and l2 and abs(l2 / l1 - 1.0) <= A4_TOL)
    l40 = ron("12b", 1, "L_d40")
    vals["A5"] = {"R_on_L_20nm": l1, "R_on_L_40nm": l40,
                  "rel_diff": None if not (l1 and l40) else abs(l40 / l1 - 1.0)}  # fmt: skip
    crit["A5"] = bool(l1 and l40 and abs(l40 / l1 - 1.0) <= A5_TOL)
    # ---- R1
    r1 = []
    for c in dsg:
        for lv in LEVELS:
            rr = runs.get(f"bracket|{c}|{lv}")
            if rr is None or rr["error"]:
                continue
            u, lo = rr["runs"]["U"], rr["runs"]["L"]
            row = {"case": c, "level": lv, "R_on_bulk_ohm_cm2": rr["R_on_bulk_ohm_cm2"],
                   "R_on_U_ohm_cm2": u.get("R_on_sp_ohm_cm2"),
                   "R_on_L_ohm_cm2": lo.get("R_on_sp_ohm_cm2"),
                   "n_acc_edges": lo["n_acc_edges"], "n_acc_columns": lo["n_acc_columns"]}  # fmt: skip
            if row["R_on_U_ohm_cm2"] and row["R_on_L_ohm_cm2"]:
                row["L_over_U_minus_1"] = row["R_on_L_ohm_cm2"] / row["R_on_U_ohm_cm2"] - 1.0
                row["U_over_bulk_minus_1"] = row["R_on_U_ohm_cm2"] / row["R_on_bulk_ohm_cm2"] - 1.0
            r1.append(row)
    # ---- R2
    r2 = {m: runs.get(f"r2|{m}|1") for m in ("bulk", "U", "L")}
    r2s = {"per_model": r2}
    if all(r2[m] and not r2[m]["error"] for m in r2):
        g = {m: r2[m]["dR_dthetaJ"] for m in r2}
        r2s["ratios"] = {"U_over_bulk": g["U"] / g["bulk"], "L_over_bulk": g["L"] / g["bulk"],
                         "L_over_U": g["L"] / g["U"]}  # fmt: skip
    runtime = time.perf_counter() - t_start
    res = {"config": {"V_ON": V_ON, "V_LIN": V_LIN, "d_ch_nm": D_CH_NM, "levels": list(LEVELS),
                      "tolerances": {"A2": A2_TOL, "A4": A4_TOL, "A5": A5_TOL},
                      "eps_fd": list(EPS_FD), "dev": bool(args.dev), "designs": dinfo,
                      "level_names": "1 = h/2, 2 = h/4"},
           "A1": a1, "runs": runs, "R1_table": r1, "R2": r2s, "runtime_seconds": runtime}  # fmt: skip
    summary = {"criteria": crit, "gate_passed": bool(all(crit.values())), "dev": bool(args.dev),
               "values": {**vals, "R1": r1, "R2": r2s, "A1": a1}, "designs": dinfo,
               "runtime_seconds": runtime}  # fmt: skip
    (out / "results.json").write_text(json.dumps(res, indent=2, default=float) + "\n")
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    figdir = REPO / "figures" if out == DEFAULT_OUT.resolve() else out / "figures"
    try:
        figures(res, figdir)
    except Exception as exc:
        print("figure error:", exc)
    print(json.dumps(crit, indent=2))
    print(f"runtime {runtime:.1f} s")


if __name__ == "__main__":
    main()
