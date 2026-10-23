"""Stage-18b benchmark (docs/benchmarks/stage18b_accumulation_excess.md): B1-B4 (gating), R1-R3.

Accumulation mobility applied to the gate-induced (excess) electron charge only on the JFET
surface edges of stage 18 (``ChannelMap(acc_mode="excess")``):
mu_edge = mu_bulk (1 - f) + mu_acc f, f = max(nbar - n_fb, 0)/nbar, mu_acc = registry
``mu_acc_n_zhao2007``; n_fb = N_D^+ - N_A^- at the node's neutral-bulk potential (incomplete
ionization of the model), averaged over the two edge nodes; mu_edge lagged in the Picard loop.
Bracket reported: [R_on(U), R_on(L')]; stage-18 L (all electrons) is read from results/stage18.

* B1 : mu_acc := mu_bulk edge by edge (excess mode) reproduces the stage-16 channel R_on (12b,
       level 1, d_ch = 20 nm) to <= 1e-12 relative.
* B2 : every Picard / Newton solve converged; fd-2 capture has no "On entry to" line.
* B3 : R_on(L') 12b level 1 vs 2 <= 2 %.   B4 : d_ch 20 vs 40 nm (12b, level 1) <= 1 %.
* R1 : bracket on levels 1, 2 for 12b, 13d initial / final, stage-17 final, 17b final (if present),
       with the stage-18 L value.
* R2 : per-column gate-induced sheet charge Q_ex = int_0^d_ch (n - n_fb)^+ dy (trapezoid on the
       in-layer nodes of each accumulation column) against C_ox (V_GS - V_FB)/q.  C_ox = eps_ox/t_ox
       (registry eps_r_SiO2_thermal, design t_ox, 1D).  Two V_FB definitions are reported: the
       device's gate-contact flat band relative to the neutral JFET surface (the one under which
       the model's oxide field vanishes: V_FB,loc = V_T (v_code + psi_0,col)), and the p-body
       V_FB of ``SiC4H.vfb(n_body)`` (the spec formula).
* R3 : Richardson FD dR_on/dtheta_J (eps = 1e-3, 5e-4) at the 13d initial design, level 1, under
       L' and U (n_fb recomputed at each theta; column/edge sets frozen at the base design).

Writes <out>/{results,summary,manifest}.json, stderr.log (default results/stage18b) and figures
stage18b_*.png.  ``--dev`` marks a DEV run (refused under results/).
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

from widegap.adjoint_sic import richardson_fd
from widegap.channel_mobility import ChannelMap, acc_mobility_from_registry, picard_solve
from widegap.manifest import build_manifest
from widegap.materials.registry import Registry
from widegap.sic import TracedRegistry

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[2]
DEFAULT_OUT = REPO / "results" / "stage18b"
R18 = REPO / "results" / "stage18" / "results.json"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


S18 = _load("stage18_run", REPO / "benchmarks/stage18/run.py")  # fixed names: fork pickling
S16, TH0 = S18.S16, S18.TH0
V_ON, V_LIN, D_CH_NM, LEVELS, K_JFET, EPS_FD = (S18.V_ON, S18.V_LIN, S18.D_CH_NM, S18.LEVELS,
                                                S18.K_JFET, S18.EPS_FD)  # fmt: skip
B1_TOL, B3_TOL, B4_TOL = 1e-12, 0.02, 0.01
rss_mib = S18.rss_mib


def solve18b(dev, sc, P, cmap, th=None):
    """Direct march + Picard (excess/U map).  Returns (R_on or nan, record, picard result)."""
    if th is not None:
        dev.ev.set_theta(np.asarray(th, float))
        dev.nd, dev.na = dev.ev.nd, dev.ev.na
        cmap.mu_bulk = np.asarray(dev.model.mu_edge, dtype=float).copy()
        cmap.refresh_doping(dev.nd, dev.na)
    dev.model.mu_edge = cmap.mu_bulk.copy()
    eq = dev.equilibrium()
    sol = S16._bias_to(dev, sc, eq.X, V_ON, V_LIN)
    ok_bulk = bool(eq.report.converged and sol.report.converged)
    vds_s = float(sc.volt_s(V_LIN))
    r = picard_solve(cmap, P, sol.X, dev.bc(float(sc.volt_s(V_ON)), 0.0, vds_s))
    rec = {"picard": S16._picard_record(r), "R_bulk": float(V_LIN * dev.x_cm
                                                          / S16._id_a_per_cm(dev, sol)),
           "bulk_converged": ok_bulk, "converged": bool(ok_bulk and r["converged"])}  # fmt: skip
    R = float("nan")
    if r["converged"]:
        R = float(V_LIN * dev.x_cm / S16._id_a_per_cm(dev, r["sol"]))
        rec["R_on_sp_ohm_cm2"] = R
        rec["counts"] = r["counts"]
    dev.model.mu_edge = cmap.mu_bulk.copy()
    return R, rec, r


def make_map(sc, dev, dnm, kind, mu_acc):
    """kind: U (stage-16), Lp (excess, Zhao), B1 (excess, mu_acc := mu_bulk)."""
    kw = dict(d_ch_um=dnm * 1e-3, gate_x_max=dev.prm.gate_x1)
    if kind == "U":
        return ChannelMap(sc, dev.model, dev.nd, dev.na, **kw)
    return ChannelMap(sc, dev.model, dev.nd, dev.na, mu_acc_cm2="bulk" if kind == "B1" else mu_acc,
                      acc_mode="excess", **kw)  # fmt: skip


def job_bracket(args):
    case, level, theta, variants, entries = args  # variants: ((tag, kind, d_ch_nm), ...)
    t0 = time.perf_counter()
    rec = {"kind": "bracket", "case": case, "level": level, "error": None}
    try:
        sc, P, dev = S18.make_dev(case, level, theta, entries)
        mu_acc = acc_mobility_from_registry(TracedRegistry(Registry(entries)))
        rec["nodes"] = int(dev.model.mesh.n_nodes)
        out = {}
        for tag, kind, dnm in variants:
            cmap = make_map(sc, dev, dnm, kind, mu_acc)
            _, r, _ = solve18b(dev, sc, P, cmap)
            r.update(S18.acc_info(cmap, sc))
            if kind != "U":
                r["n_fb_edge_range_cm3"] = [float(cmap.nfb_edge.min() * sc.n0),
                                            float(cmap.nfb_edge.max() * sc.n0)]  # fmt: skip
            out[tag] = r
        rec["runs"] = out
        rec["R_on_bulk_ohm_cm2"] = next(iter(out.values()))["R_bulk"]
    except Exception as exc:  # recorded, never hidden
        rec["error"] = f"{type(exc).__name__}: {exc}"
    rec["seconds"] = time.perf_counter() - t0
    rec["rss_mib"] = rss_mib()
    return rec


def job_sheet(args):
    """R2: per-column excess sheet charge vs C_ox (V_GS - V_FB)/q (L', level 1)."""
    case, level, theta, entries = args
    t0 = time.perf_counter()
    rec = {"kind": "sheet", "case": case, "level": level, "error": None}
    try:
        sc, P, dev = S18.make_dev(case, level, theta, entries)
        mu_acc = acc_mobility_from_registry(TracedRegistry(Registry(entries)))
        cmap = make_map(sc, dev, D_CH_NM, "Lp", mu_acc)
        _, r0, r = solve18b(dev, sc, P, cmap)
        rec["picard"] = r0["picard"]
        n = np.asarray(dev.model.densities(r["sol"].X)[0])
        xy = cmap.mesh.node_xy
        ix = cmap.node_col
        from widegap.sic import EPS0

        eps_ox = sc.reg.value("eps_r_SiO2_thermal") * EPS0 * 1e-2  # F/cm
        t_ox_cm = float(sc.length_cm(dev.prm.t_ox))
        c_ox = eps_ox / t_ox_cm  # F/cm^2
        q = 1.602176634e-19
        v_fb_body, v_code, _ = sc.vfb(S18.S15F.DESIGN.n_body)
        cols = []
        for i in np.flatnonzero(cmap.col_acc):
            nodes = np.flatnonzero(cmap.in_ch & (ix == i))
            nodes = nodes[np.argsort(-xy[nodes, 1])]
            depth = -xy[nodes, 1] * sc.l_d  # cm
            ex = np.maximum(n[nodes] - cmap.nfb_node[nodes], 0.0) * sc.n0
            q_ex = float(np.sum(0.5 * (ex[1:] + ex[:-1]) * np.diff(depth)))
            psi0 = float(cmap.psi0_node[nodes[0]])
            v_fb_loc = float(sc.vt * (v_code + psi0))
            ref_loc = c_ox * (V_ON - v_fb_loc) / q
            ref_body = c_ox * (V_ON - v_fb_body) / q
            cols.append({"x_um": float(sc.length_cm(cmap.mesh.x[i])) * 1e4, "Q_ex_cm2": q_ex,
                         "V_FB_local_V": v_fb_loc, "ref_local_cm2": ref_loc,
                         "ratio_local": q_ex / ref_loc, "ref_body_cm2": ref_body,
                         "ratio_body": q_ex / ref_body})  # fmt: skip
        rec.update({"columns": cols, "C_ox_F_cm2": c_ox, "V_FB_body_V": v_fb_body,
                    "t_ox_um": t_ox_cm * 1e4, "V_GS": V_ON})  # fmt: skip
    except Exception as exc:
        rec["error"] = f"{type(exc).__name__}: {exc}"
    rec["seconds"] = time.perf_counter() - t0
    rec["rss_mib"] = rss_mib()
    return rec


def job_r3(args):
    model, level, entries = args  # "U" or "Lp"
    t0 = time.perf_counter()
    rec = {"kind": "r3", "model": model, "level": level, "error": None}
    try:
        sc, P, dev = S18.make_dev("13d_initial", level, None, entries)
        mu_acc = acc_mobility_from_registry(TracedRegistry(Registry(entries)))
        cmap = make_map(sc, dev, D_CH_NM, model, mu_acc)
        th0 = np.asarray(TH0, float)
        flags = []

        def f(th):
            R, d, _ = solve18b(dev, sc, P, cmap, th)
            flags.append(bool(d["converged"]))
            return R

        r0 = f(th0)
        u = np.zeros(5)
        u[K_JFET] = 1.0
        rich, f1, f2 = richardson_fd(f, th0, u, EPS_FD[0])
        rec.update({"R0": r0, "dR_dthetaJ": float(rich), "fd_eps_1e-3": float(f1),
                    "fd_eps_5e-4": float(f2), "spread": float(abs(f1 - f2)),
                    "n_solves": len(flags), "all_converged": all(flags),
                    **S18.acc_info(cmap, sc)})  # fmt: skip
    except Exception as exc:
        rec["error"] = f"{type(exc).__name__}: {exc}"
    rec["seconds"] = time.perf_counter() - t0
    rec["rss_mib"] = rss_mib()
    return rec


def figures(res, figdir: Path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figdir.mkdir(parents=True, exist_ok=True)
    rows = [r for r in res["R1_table"] if r["R_on_U_ohm_cm2"] and r["R_on_Lp_ohm_cm2"]]
    fig, ax = plt.subplots(figsize=(9, 4.2))
    xs = np.arange(len(rows))
    ax.bar(xs - 0.27, [r["R_on_U_ohm_cm2"] * 1e4 for r in rows], 0.27, label="(U) bulk")
    ax.bar(xs, [r["R_on_Lp_ohm_cm2"] * 1e4 for r in rows], 0.27, label="(L') excess charge only")
    ax.bar(xs + 0.27, [r.get("R_on_L18_ohm_cm2") or np.nan for r in rows], 0.27,
           label="stage-18 L (all electrons)")  # fmt: skip
    ax.set_xticks(xs, [f"{r['case']}\nlevel {r['level']}" for r in rows], fontsize=7)
    ax.set_ylabel(r"$R_{on,sp}$ [$10^{-4}\ \Omega\,$cm$^2$]")
    ax.set_title("Stage 18b: bracket [R_on(U), R_on(L')]")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(figdir / "stage18b_bracket.png", dpi=130)
    plt.close(fig)
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.8))
    for ax, (k, r) in zip(axs, [(k, r) for k, r in res["runs"].items() if r["kind"] == "sheet"][:2],
                          strict=False):  # fmt: skip
        if r["error"]:
            continue
        x = [c["x_um"] for c in r["columns"]]
        ax.plot(x, [c["Q_ex_cm2"] for c in r["columns"]], "ko-", label="Q_ex (model)")
        ax.plot(x, [c["ref_local_cm2"] for c in r["columns"]], "r--", label="Cox(Vgs-VFB,loc)/q")
        ax.plot(x, [c["ref_body_cm2"] for c in r["columns"]], "b:", label="Cox(Vgs-VFB,body)/q")
        ax.set_title(k)
        ax.set_xlabel("x [um]")
        ax.set_ylabel("sheet charge [cm$^{-2}$]")
        ax.legend(fontsize=7)
        ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(figdir / "stage18b_sheet_charge.png", dpi=130)
    plt.close(fig)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--workers", type=int, default=min(4, os.cpu_count() or 1))
    ap.add_argument("--stage17b-results", type=Path, default=S18.R17B)
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
    entries = S18.registry_entries()
    dsg, dinfo = S18.designs(args.stage17b_results)
    jobs = []
    for lv in LEVELS:
        for c, th in dsg.items():
            var = [("U", "U", D_CH_NM), ("Lp", "Lp", D_CH_NM)]
            if c == "12b" and lv == 1:
                var += [("Lp_d40", "Lp", 40.0), ("B1", "B1", D_CH_NM)]
            jobs.append((job_bracket, ("bracket", c, lv), (c, lv, th, tuple(var), entries)))
    for c in ("12b", "13d_initial"):
        jobs.append((job_sheet, ("sheet", c, 1), (c, 1, dsg[c], entries)))
    for m in ("U", "Lp"):
        jobs.append((job_r3, ("r3", m, 1), (m, 1, entries)))
    jobs.sort(key=lambda j: -j[1][2] if j[1][0] == "bracket" else -1)
    runs: dict = {}
    with ProcessPoolExecutor(args.workers, mp_context=mp.get_context("fork")) as pool:
        futs = {pool.submit(fn, a): key for fn, key, a in jobs}
        for f in as_completed(futs):
            r = f.result()
            runs["|".join(str(k) for k in futs[f])] = r
            print(
                f"done {futs[f]} {r['seconds']:.0f}s rss={r['rss_mib']:.0f}MiB err={r['error']}",
                flush=True,
            )

    def ron(case, lv, tag):
        r = runs.get(f"bracket|{case}|{lv}")
        return None if not r or r["error"] else r["runs"].get(tag, {}).get("R_on_sp_ohm_cm2")

    crit, vals = {}, {}
    s16 = json.loads(S18.R16.read_text())["runs"]["ron|12b|1"]["channel"]["20"]["R_on_sp_ohm_cm2"]
    b1 = ron("12b", 1, "B1")
    vals["B1"] = {"R_on_B1_this_run": b1, "R_on_stage16": s16,
                  "rel_diff": None if b1 is None else abs(b1 / s16 - 1.0)}  # fmt: skip
    crit["B1"] = bool(b1 is not None and abs(b1 / s16 - 1.0) <= B1_TOL)
    failed = []
    for k, r in runs.items():
        if r["error"]:
            failed.append(k)
        elif r["kind"] == "bracket":
            failed += [f"{k}@{t}" for t, e in r["runs"].items() if not e["converged"]]
        elif r["kind"] == "sheet":
            if not r["picard"]["converged"]:
                failed.append(k)
        elif not r["all_converged"]:
            failed.append(k)
    sys.stderr.flush()
    text = log_path.read_text(errors="replace") if log_path.exists() else ""
    entry_lines = [ln for ln in text.splitlines() if "On entry to" in ln]
    vals["B2"] = {
        "failed": failed,
        "on_entry_to_lines": len(entry_lines),
        "stderr_bytes": len(text),
    }
    crit["B2"] = bool(not failed and not entry_lines and runs)
    l1, l2, l40 = ron("12b", 1, "Lp"), ron("12b", 2, "Lp"), ron("12b", 1, "Lp_d40")
    vals["B3"] = {"R_on_Lp_level1": l1, "R_on_Lp_level2": l2,
                  "rel_diff": None if not (l1 and l2) else abs(l2 / l1 - 1.0)}  # fmt: skip
    crit["B3"] = bool(l1 and l2 and abs(l2 / l1 - 1.0) <= B3_TOL)
    vals["B4"] = {"R_on_Lp_20nm": l1, "R_on_Lp_40nm": l40,
                  "rel_diff": None if not (l1 and l40) else abs(l40 / l1 - 1.0)}  # fmt: skip
    crit["B4"] = bool(l1 and l40 and abs(l40 / l1 - 1.0) <= B4_TOL)
    r18 = json.loads(R18.read_text())["runs"] if R18.exists() else {}
    r1 = []
    for c in dsg:
        for lv in LEVELS:
            rr = runs.get(f"bracket|{c}|{lv}")
            if rr is None or rr["error"]:
                continue
            u, lp = rr["runs"]["U"], rr["runs"]["Lp"]
            l18 = r18.get(f"bracket|{c}|{lv}", {}).get("runs", {}).get("L", {})
            row = {"case": c, "level": lv, "R_on_bulk_ohm_cm2": rr["R_on_bulk_ohm_cm2"],
                   "R_on_U_ohm_cm2": u.get("R_on_sp_ohm_cm2"),
                   "R_on_Lp_ohm_cm2": lp.get("R_on_sp_ohm_cm2"),
                   "R_on_L18_ohm_cm2": l18.get("R_on_sp_ohm_cm2"),
                   "n_acc_edges": lp["n_acc_edges"]}  # fmt: skip
            if row["R_on_U_ohm_cm2"] and row["R_on_Lp_ohm_cm2"]:
                row["Lp_over_U_minus_1"] = row["R_on_Lp_ohm_cm2"] / row["R_on_U_ohm_cm2"] - 1.0
                if row["R_on_L18_ohm_cm2"]:
                    row["L18_over_U_minus_1"] = (
                        row["R_on_L18_ohm_cm2"] / row["R_on_U_ohm_cm2"] - 1.0
                    )
            r1.append(row)
    r2 = {k: v for k, v in runs.items() if v["kind"] == "sheet"}
    r3 = {m: runs.get(f"r3|{m}|1") for m in ("U", "Lp")}
    r3s = {"per_model": r3}
    if all(r3[m] and not r3[m]["error"] for m in r3):
        r3s["Lp_over_U"] = r3["Lp"]["dR_dthetaJ"] / r3["U"]["dR_dthetaJ"]
    runtime = time.perf_counter() - t_start
    res = {"config": {"V_ON": V_ON, "V_LIN": V_LIN, "d_ch_nm": D_CH_NM, "levels": list(LEVELS),
                      "tolerances": {"B1": B1_TOL, "B3": B3_TOL, "B4": B4_TOL},
                      "dev": bool(args.dev), "designs": dinfo},
           "runs": runs, "R1_table": r1, "R2": r2, "R3": r3s, "runtime_seconds": runtime}  # fmt: skip
    summary = {"criteria": crit, "gate_passed": bool(all(crit.values())), "dev": bool(args.dev),
               "values": {**vals, "R1": r1, "R3": r3s}, "designs": dinfo,
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
