"""Stage 11: do PINNs add value next to the verified FV/adjoint stack?
(docs/benchmarks/stage11_pinn_question.md; criteria M1-M5, R1-R3).

Writes to --out (default results/stage11):
  fv.json            FV references, M1/M3/M5, timings, convergence ladders
  fv_rss_*.json      FV peak RSS (one subprocess per mesh)
  calibration.json   isolated Adam-step time (load normalisation, see below)
  pinn_A_seed<k>.json/.npz, pinn_B_seed<k>.json/.npz   one subprocess per seed
  r1.json, r2.json, summary.json, manifest.json, figures/stage11_*.png

Every PINN seed runs in its own subprocess (peak RSS attributable to it).
The PINN configuration is FIXED by the spec; the only options that change it
(--adam-steps, --lbfgs-iter) mark the run as a non-official dev run.

Parallelism: --jobs J runs J PINN subprocesses concurrently with
threads = cpu_count // J each, to fit the ~60 min budget. Wall times are then
measured under load. To keep the cost comparison fair, an isolated
Adam-step time (all threads, nothing else running) is measured first and each
run's training time is also reported scaled by (isolated step time / that run's
loaded step time); both raw and scaled numbers are in the results. With
--jobs 1 (sequential, all threads) the two coincide.

SPEC ISSUE (documented, see widegap.pinn_study): for the symmetric doping
C = -+a the solution is odd and J = 0 for every a, so Q-B scales the n side
only (C = -1 for x<0, +a for x>0). The symmetric degeneracy is measured and
reported in fv.json ("symmetric_degeneracy").
"""

from __future__ import annotations

import argparse
import json
import os
import resource
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
SEEDS = [0, 1, 2, 3, 4]
A_POINTS = np.linspace(0.5, 2.0, 9)
N_EVAL = 4001
N_REF = 2**18
N_FINEST_PAIR = (2**17, 2**18)
LADDER = [2**k for k in range(4, 17)]  # 16 .. 65536 intervals
M_LIST = [1, 10, 100, 1000, 10000]
THRESHOLDS = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5]
SELF = Path(__file__).resolve()


def _dump(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=2, default=float) + "\n")


def _rss_mib() -> float:
    """Peak RSS of this process. VmHWM, not ru_maxrss: after fork+exec ru_maxrss can
    still carry the (larger) parent's peak, which would inflate the subprocess."""
    for line in Path("/proc/self/status").read_text().splitlines():
        if line.startswith("VmHWM"):
            return float(line.split()[1]) / 1024.0
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def _vmrss_mib() -> float:
    for line in Path("/proc/self/status").read_text().splitlines():
        if line.startswith("VmRSS"):
            return float(line.split()[1]) / 1024.0
    return float("nan")


def _median_time(fn, reps: int) -> float:
    ts = []
    for _ in range(reps):
        t0 = time.perf_counter()
        fn()
        ts.append(time.perf_counter() - t0)
    return float(np.median(ts))


def _ncpu() -> int:
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        return os.cpu_count() or 1


def _thread_info() -> dict:
    info = {
        "cpu_count": os.cpu_count(),
        "usable_cpus": _ncpu(),
        "env": {
            k: os.environ.get(k)
            for k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS")
        },
    }
    if "torch" in sys.modules:
        import torch

        info["torch_threads"] = torch.get_num_threads()
        info["torch_interop_threads"] = torch.get_num_interop_threads()
    return info


# ============================================================ FV phase
def _rel_l2(x, a, b):
    from widegap.verification import rel_l2_error

    return rel_l2_error(x, a, b)


def _m1_pair(a_n: float) -> dict:
    """Two finest meshes (nested: coarse nodes are every second fine node)."""
    from widegap import pinn_study as ps

    n_c, n_f = N_FINEST_PAIR
    rc = ps.fv_J_and_grad(n_c, a_n)
    rf = ps.fv_J_and_grad(n_f, a_n)
    pc, pf = rc["psi"], rf["psi"][::2]
    x = rc["x"]
    abs_psi = float(np.sum(np.abs(rf["psi"]) * np.r_[np.diff(rf["x"]), 0.0]))  # ~ int|psi|
    out = {
        "a": a_n,
        "meshes": [n_c, n_f],
        "psi_max_rel": float(np.max(np.abs(pf - pc)) / np.max(np.abs(rf["psi"]))),
        "psi_l2_rel": _rel_l2(x, pc, pf),
        "J_coarse": rc["J"],
        "J_fine": rf["J"],
        "J_change_abs": abs(rf["J"] - rc["J"]),
        "J_change_rel_to_absJ": abs(rf["J"] - rc["J"]) / max(abs(rf["J"]), 1e-300),
        "J_change_rel_to_int_abs_psi": abs(rf["J"] - rc["J"]) / abs_psi,
        "dJda_coarse": rc["dJda"],
        "dJda_fine": rf["dJda"],
        "dJda_change_rel": abs(rf["dJda"] - rc["dJda"]) / abs(rf["dJda"]),
        "newton_converged": [rc["newton_converged"], rf["newton_converged"]],
    }
    return out, rf


def _m3(res_fine: dict, xe: np.ndarray) -> dict:
    """Interpolation error of the FV reference on the evaluation grid.

    (i) bound h^2/8 max|psi''| from nodal second differences;
    (ii) empirical: interpolate from every second node vs all nodes. The coarse
    interpolation error is ~4x the fine one, so fine ~ diff / 3.
    """
    from widegap import pinn_study as ps

    x, psi = res_fine["x"], res_fine["psi"]
    full = np.interp(xe, x, psi)
    half = np.interp(xe, x[::2], psi[::2])
    return {
        "n_eval": len(xe),
        "interp_error_bound": ps.interp_error_bound(x, psi),
        "interp_error_empirical_est": float(np.max(np.abs(full - half)) / 3.0),
    }


def _ladder_A(ref_eval: np.ndarray, xe: np.ndarray) -> list[dict]:
    from widegap import pinn_study as ps
    from widegap.verification import rel_l2_error

    rows = []
    for n in LADDER:
        sol, _ = ps.fv_junction(n)
        t = _median_time(lambda n=n: ps.fv_junction(n), 5 if n <= 16384 else 3)
        p = np.interp(xe, sol.x, sol.psi)
        rows.append(
            {
                "N": n,
                "seconds": t,
                "rel_l2": rel_l2_error(xe, p, ref_eval),
                "max_err": float(np.max(np.abs(p - ref_eval))),
                "converged": sol.report.converged,
            }
        )
    return rows


def _ladder_B(j_ref: np.ndarray, g_ref: np.ndarray) -> list[dict]:
    from widegap import pinn_study as ps

    rows = []
    for n in LADDER:
        ej, eg, ts = [], [], []
        for k, a in enumerate(A_POINTS):
            r = ps.fv_J_and_grad(n, float(a))
            reps = 5 if n <= 4096 else 2
            ts.append(_median_time(lambda n=n, a=a: ps.fv_J_and_grad(n, float(a)), reps))
            ej.append(abs(r["J"] - j_ref[k]) / abs(j_ref[k]))
            eg.append(abs(r["dJda"] - g_ref[k]) / abs(g_ref[k]))
        rows.append(
            {
                "N": n,
                "seconds_per_eval": float(np.median(ts)),
                "J_rel_err": ej,
                "dJda_rel_err": eg,
                "J_rel_err_max": max(ej),
                "dJda_rel_err_max": max(eg),
            }
        )
    return rows


def _m5(n: int = 1024) -> dict:
    from widegap import pinn_study as ps

    rows = []
    for a in A_POINTS:
        r = ps.fv_J_and_grad(n, float(a))
        _, g_lib = ps.fv_dJ_dC_library(n, float(a))
        fd = ps.richardson_fd_dJda(n, float(a))
        rows.append(
            {
                "a": float(a),
                "adjoint_dJda": r["dJda"],
                "adjoint_interior_part": r["dJda_interior"],
                "adjoint_contact_part": r["dJda_contact"],
                "fd_R2": fd["R2"],
                "fd_R1": fd["R1"][-1],
                "fd_error_proxy_R2_minus_R1": fd["error_proxy"],
                "rel_err_vs_R2": abs(r["dJda"] - fd["R2"]) / abs(fd["R2"]),
                "rel_err_vs_R1": abs(r["dJda"] - fd["R1"][-1]) / abs(fd["R1"][-1]),
                "dJdC_max_diff_vs_adjoint_module": float(np.max(np.abs(r["dJdC"] - g_lib))),
            }
        )
    return {
        "mesh": n,
        "rel_h": 0.02,
        "levels": 3,
        "rows": rows,
        "max_rel_err_vs_R2": max(r["rel_err_vs_R2"] for r in rows),
    }


def run_fv(out: Path) -> dict:
    from widegap import pinn_study as ps

    xe = ps.eval_grid(N_EVAL)
    fv: dict = {"eval_points": N_EVAL, "a_points": list(A_POINTS), "n_ref": N_REF}

    # ---- Q-A reference (a = 1) and M1/M3
    t0 = time.perf_counter()
    m1a, ref = _m1_pair(1.0)
    fv["ref_seconds_two_meshes"] = time.perf_counter() - t0
    fv["M1_A"] = m1a
    fv["M3_A"] = _m3(ref, xe)
    ref_eval = np.interp(xe, ref["x"], ref["psi"])
    np.savez(out / "ref_A.npz", xe=xe, psi=ref_eval)

    # ---- Q-B references at the 9 points, M1/M3
    m1b, m3b, psi_b, j_ref, g_ref = [], [], [], [], []
    for a in A_POINTS:
        m, rf = _m1_pair(float(a))
        m1b.append(m)
        m3b.append(_m3(rf, xe))
        psi_b.append(np.interp(xe, rf["x"], rf["psi"]))
        j_ref.append(rf["J"])
        g_ref.append(rf["dJda"])
    fv["M1_B"], fv["M3_B"] = m1b, m3b
    fv["J_ref"], fv["dJda_ref"] = j_ref, g_ref
    np.savez(
        out / "ref_B.npz",
        xe=xe,
        a=A_POINTS,
        J=np.array(j_ref),
        dJda=np.array(g_ref),
        psi=np.array(psi_b),
    )

    # ---- symmetric-scaling degeneracy (spec issue)
    sym = []
    for a in (0.5, 1.25, 2.0):
        sol, _ = ps.fv_junction(4096, a, a)
        sym.append({"a": a, "J": ps.trapezoid_J(sol.x, sol.psi)})
    fv["symmetric_degeneracy"] = sym

    # ---- M5
    fv["M5"] = _m5()

    # ---- timings and ladders
    fv["timing_ref_mesh_forward_s"] = _median_time(lambda: ps.fv_junction(N_REF), 3)
    fv["timing_ref_mesh_forward_plus_adjoint_s"] = _median_time(
        lambda: ps.fv_J_and_grad(N_REF, 1.25), 3
    )
    fv["ladder_A"] = _ladder_A(ref_eval, xe)
    fv["ladder_B"] = _ladder_B(np.array(j_ref), np.array(g_ref))
    fv["threads"] = _thread_info()
    return fv


# ============================================================ subprocess workers
def single_fv(n: int) -> dict:
    from widegap import pinn_study as ps

    base = _vmrss_mib()
    t0 = time.perf_counter()
    r = ps.fv_J_and_grad(n, 1.25)
    return {
        "N": n,
        "seconds_forward_plus_adjoint": time.perf_counter() - t0,
        "rss_after_imports_mib": base,
        "peak_rss_mib": _rss_mib(),
        "J": r["J"],
        "threads": _thread_info(),
    }


def calibrate() -> dict:
    import torch

    from widegap import pinn_study as ps

    out = {}
    for task, par in (("A", False), ("B", True)):
        net = ps.build_net((64, 64, 64, 64), 2 if par else 1)
        x = (torch.rand(2048, 1, dtype=torch.float64) * 2 - 1) * 40
        a = torch.rand(2048, 1, dtype=torch.float64) * 1.5 + 0.5 if par else None
        opt = torch.optim.Adam(net.parameters())
        for k in range(110):
            if k == 10:
                t0 = time.perf_counter()
            opt.zero_grad()
            ps.residual(net, x, a).pow(2).mean().backward()
            opt.step()
        out[f"{task}_adam_step_s"] = (time.perf_counter() - t0) / 100
    out["threads"] = _thread_info()
    return out


def single_pinn(task: str, seed: int, args) -> dict:
    import torch

    from widegap import pinn_study as ps
    from widegap.verification import rel_l2_error

    if args.threads:
        torch.set_num_threads(args.threads)
    base = _vmrss_mib()
    out_dir = Path(args.out)
    cfg_kw = {"seed": seed}
    mult = 1 if task == "A" else 2
    if args.adam_steps is not None:
        cfg_kw["adam_steps"] = args.adam_steps * mult
    if args.lbfgs_iter is not None:
        cfg_kw["lbfgs_max_iter"] = args.lbfgs_iter * mult
    base_cfg = ps.StudyConfig() if task == "A" else ps.parametric_config()
    cfg = ps.StudyConfig(**{**base_cfg.__dict__, **cfg_kw})

    rec: dict = {"task": task, "seed": seed, "config": cfg.__dict__ | {"hidden": list(cfg.hidden)}}
    if task == "A":
        ref = np.load(out_dir / "ref_A.npz")
        xe, psi_ref = ref["xe"], ref["psi"]

        def cb(net):
            p = ps.predict_psi(net, xe)
            return {
                "rel_l2": rel_l2_error(xe, p, psi_ref),
                "max_err": float(np.max(np.abs(p - psi_ref))),
            }

        res = ps.train(cfg, parametric=False, checkpoint=cb)
        p = ps.predict_psi(res.model, xe)
        rec.update(
            rel_l2=rel_l2_error(xe, p, psi_ref),
            max_err=float(np.max(np.abs(p - psi_ref))),
            max_err_rel_to_max_psi=float(np.max(np.abs(p - psi_ref)) / np.max(np.abs(psi_ref))),
            checkpoints=res.checkpoints,
        )
        np.savez(out_dir / f"pinn_A_seed{seed}.npz", xe=xe, psi=p)
    else:
        ref = np.load(out_dir / "ref_B.npz")
        xe, a_pts, j_ref, g_ref, psi_ref = (ref[k] for k in ("xe", "a", "J", "dJda", "psi"))
        res = ps.train(cfg, parametric=True)
        j, g = ps.surrogate_J_and_grad(res.model, a_pts)
        psi_l2 = [
            rel_l2_error(xe, ps.predict_psi(res.model, xe, float(a)), psi_ref[k])
            for k, a in enumerate(a_pts)
        ]
        # quadrature check of the surrogate J (1001-point trapezoid vs 8001)
        j8, _ = ps.surrogate_J_and_grad(res.model, a_pts, n_quad=8001)
        rec.update(
            J=list(j),
            dJda=list(g),
            J_rel_err=list(np.abs(j - j_ref) / np.abs(j_ref)),
            J_err_over_int_abs_psi=list(
                np.abs(j - j_ref) / (np.trapezoid(np.abs(psi_ref), xe, axis=1))
                if hasattr(np, "trapezoid")
                else np.abs(j - j_ref) / np.trapz(np.abs(psi_ref), xe, axis=1)
            ),
            dJda_rel_err=list(np.abs(g - g_ref) / np.abs(g_ref)),
            psi_rel_l2=psi_l2,
            quadrature_check_max_abs_diff_1001_vs_8001=float(np.max(np.abs(j - j8))),
        )
        # inference timing: M queries of (J, dJ/da); median of repeats
        inf = {}
        for m in (1, 10, 100):
            aq = np.linspace(0.5, 2.0, m)
            ps.surrogate_J_and_grad(res.model, aq)  # warm-up
            inf[str(m)] = _median_time(lambda aq=aq: ps.surrogate_J_and_grad(res.model, aq), 5)
        rec["inference_seconds"] = inf
    rec.update(
        train_seconds=res.train_seconds,
        adam_seconds=res.adam_seconds,
        adam_step_seconds=res.adam_seconds / cfg.adam_steps if cfg.adam_steps else None,
        final_residual_mse=res.final_loss,
        n_loss_evals=res.n_loss_evals,
        finite=res.finite,
        loss_history_every_100=res.loss_history[::100],
        rss_after_imports_mib=base,
        peak_rss_mib=_rss_mib(),
        threads=_thread_info(),
    )
    return rec


# ============================================================ assembling
def _first_below(vals, thr):
    for t, v in vals:
        if v <= thr:
            return t
    return None


def assemble_r1(fv: dict, pinn_a: list[dict], factors: dict) -> dict:
    ladder = fv["ladder_A"]
    rows = []
    for thr in THRESHOLDS:
        ok = [r for r in ladder if all(q["rel_l2"] <= thr for q in ladder if q["N"] >= r["N"])]
        f = ok[0] if ok else None
        row = {
            "target_rel_l2": thr,
            "fv_N": f["N"] if f else None,
            "fv_seconds": f["seconds"] if f else None,
            "pinn_seconds_raw": {},
            "pinn_seconds_scaled": {},
        }
        for rec in pinn_a:
            if rec.get("failed"):
                continue
            s = rec["seed"]
            t = _first_below([(c["train_seconds"], c["rel_l2"]) for c in rec["checkpoints"]], thr)
            row["pinn_seconds_raw"][s] = t
            row["pinn_seconds_scaled"][s] = None if t is None else t * factors.get(s, 1.0)
        rows.append(row)
    return {
        "per_seed": [
            {
                k: rec.get(k)
                for k in (
                    "seed",
                    "rel_l2",
                    "max_err",
                    "max_err_rel_to_max_psi",
                    "train_seconds",
                    "peak_rss_mib",
                    "final_residual_mse",
                    "finite",
                    "failed",
                    "error",
                )
            }
            for rec in pinn_a
        ],
        "cost_to_accuracy": rows,
    }


def _match_N(ladder_B, key, target):
    ok = [r for r in ladder_B if all(q[key] <= target for q in ladder_B if q["N"] >= r["N"])]
    return ok[0] if ok else None


def assemble_r2(fv: dict, pinn_b: list[dict], factors: dict) -> dict:
    ladder = fv["ladder_B"]
    per_seed = []
    for rec in pinn_b:
        if rec.get("failed"):
            per_seed.append({"seed": rec["seed"], "failed": True, "error": rec.get("error")})
            continue
        s = rec["seed"]
        ej, eg = max(rec["J_rel_err"]), max(rec["dJda_rel_err"])
        mj = _match_N(ladder, "J_rel_err_max", ej)
        mg = _match_N(ladder, "dJda_rel_err_max", eg)
        if mj is None or mg is None:
            n_match, t_fv = None, None
        else:
            n_match = max(mj["N"], mg["N"])
            t_fv = next(r["seconds_per_eval"] for r in ladder if r["N"] == n_match)
        inf = rec["inference_seconds"]
        t_q = inf["100"] / 100.0
        t_train_raw = rec["train_seconds"]
        t_train = t_train_raw * factors.get(s, 1.0)

        def inf_total(m, inf=inf, t_q=t_q):
            if m <= 1:
                return inf["1"]
            if m <= 10:
                return inf["10"] * m / 10.0
            return inf["100"] * m / 100.0

        costs = {}
        for m in M_LIST:
            costs[str(m)] = {
                "fv": None if t_fv is None else m * t_fv,
                "pinn_scaled": t_train + inf_total(m),
                "pinn_raw": t_train_raw + inf_total(m),
            }
        be = be_raw = None
        if t_fv is not None and t_fv > t_q:
            be = t_train / (t_fv - t_q)
            be_raw = t_train_raw / (t_fv - t_q)
        per_seed.append(
            {
                "seed": s,
                "J_rel_err_max": ej,
                "J_rel_err_median": float(np.median(rec["J_rel_err"])),
                "dJda_rel_err_max": eg,
                "dJda_rel_err_median": float(np.median(rec["dJda_rel_err"])),
                "J_rel_err": rec["J_rel_err"],
                "dJda_rel_err": rec["dJda_rel_err"],
                "psi_rel_l2": rec["psi_rel_l2"],
                "fv_matched_N": n_match,
                "fv_matched_N_for_J": mj["N"] if mj else None,
                "fv_matched_N_for_dJda": mg["N"] if mg else None,
                "fv_seconds_per_eval_at_matched_N": t_fv,
                "pinn_train_seconds_raw": t_train_raw,
                "pinn_train_seconds_scaled": t_train,
                "pinn_seconds_per_query": t_q,
                "pinn_inference_seconds": inf,
                "break_even_M": be,
                "break_even_M_raw_times": be_raw,
                "J_err_over_int_abs_psi_max": max(rec["J_err_over_int_abs_psi"]),
                "break_even_exists": be is not None,
                "costs_by_M": costs,
                "peak_rss_mib": rec["peak_rss_mib"],
            }
        )
    return {"per_seed": per_seed, "M_list": M_LIST}


def _fmt(x, p=2):
    return "n/a" if x is None else f"{x:.{p}e}"


def conclusion(fv: dict, r1: dict, r2: dict) -> list[str]:
    lines = []
    ok_a = [s for s in r1["per_seed"] if not s.get("failed")]
    if ok_a:
        best = min(s["rel_l2"] for s in ok_a)
        worst = max(s["rel_l2"] for s in ok_a)
        fvb = [r for r in fv["ladder_A"] if r["rel_l2"] <= best]
        fvr = min(fvb, key=lambda r: r["N"]) if fvb else None
        t_pinn = float(np.median([s["train_seconds"] for s in ok_a]))
        lines.append(
            f"Q-A: PINN rel L2 error over {len(ok_a)} seeds: best {_fmt(best)}, worst "
            f"{_fmt(worst)}. FV reaches the best PINN error with N="
            f"{fvr['N'] if fvr else 'n/a'} in {_fmt(fvr['seconds'] if fvr else None)} s versus a "
            f"median PINN training time of {t_pinn:.0f} s."
        )
    ok_b = [s for s in r2["per_seed"] if not s.get("failed")]
    if ok_b:
        bes = [s["break_even_M"] for s in ok_b if s["break_even_exists"]]
        lines.append(
            "Q-B: surrogate max rel error over 9 points, J: "
            f"{_fmt(min(s['J_rel_err_max'] for s in ok_b))}-"
            f"{_fmt(max(s['J_rel_err_max'] for s in ok_b))}, dJ/da: "
            f"{_fmt(min(s['dJda_rel_err_max'] for s in ok_b))}-"
            f"{_fmt(max(s['dJda_rel_err_max'] for s in ok_b))}."
        )
        if bes:
            lines.append(
                f"Break-even M at matched accuracy exists for {len(bes)}/{len(ok_b)} seeds "
                f"(min {_fmt(min(bes), 3)}, max {_fmt(max(bes), 3)}); "
                f"within M<=1e4: {sum(b <= 1e4 for b in bes)}/{len(ok_b)} seeds."
            )
        else:
            lines.append(
                "No break-even M found: FV at the matched accuracy is not slower per "
                "evaluation than PINN inference (or the ladder cannot match PINN accuracy)."
            )
    lines.append(
        "Scope: one 1D problem, one fixed PINN configuration (no tuning), CPU only; "
        "the conclusion does not extend to other problems."
    )
    return lines


# ============================================================ figures
def figures(out: Path, fv, r1, r2, pinn_a, pinn_b, man_text: str, save_repo: bool):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cols = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00"]
    figdirs = [out / "figures"] + ([REPO / "figures"] if save_repo else [])
    for d in figdirs:
        d.mkdir(parents=True, exist_ok=True)

    def save(fig, name):
        fig.text(0.5, 0.005, man_text, ha="center", fontsize=7, color="0.4")
        for d in figdirs:
            fig.savefig(d / name, dpi=120)
        plt.close(fig)

    ref = np.load(out / "ref_A.npz")
    xe, psi_ref = ref["xe"], ref["psi"]
    # --- Q-A solution and error
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    ax[0].plot(xe, psi_ref, "k", lw=2, label="FV (N=2^18)")
    for k, rec in enumerate(pinn_a):
        f = out / f"pinn_A_seed{rec['seed']}.npz"
        if rec.get("failed") or not f.exists():
            continue
        p = np.load(f)["psi"]
        ax[0].plot(xe, p, "--", color=cols[k % 5], lw=1, label=f"PINN seed {rec['seed']}")
        ax[1].semilogy(xe, np.abs(p - psi_ref) + 1e-16, color=cols[k % 5], lw=1)
    ax[0].set_xlabel("x  [L_D]")
    ax[0].set_ylabel(r"$\psi$  [V$_T$]")
    ax[0].legend(fontsize=7)
    ax[0].set_title("Q-A: pn junction potential")
    ax[1].set_xlabel("x  [L_D]")
    ax[1].set_ylabel(r"$|\psi_{PINN}-\psi_{FV}|$")
    ax[1].set_title("pointwise error, all seeds")
    save(fig, "stage11_qa_solution.png")

    # --- Q-A cost to accuracy
    fig, ax = plt.subplots(figsize=(6.5, 4.5), constrained_layout=True)
    lad = fv["ladder_A"]
    ax.loglog(
        [r["seconds"] for r in lad], [r["rel_l2"] for r in lad], "ks-", label="FV (mesh ladder)"
    )
    for k, rec in enumerate(pinn_a):
        if rec.get("failed"):
            continue
        ck = rec["checkpoints"]
        ax.loglog(
            [c["train_seconds"] for c in ck],
            [c["rel_l2"] for c in ck],
            "-",
            color=cols[k % 5],
            lw=1,
            label=f"PINN seed {rec['seed']}",
        )
    ax.set_xlabel("wall time  [s]")
    ax.set_ylabel("relative L2 error vs FV reference")
    ax.set_title("Q-A: cost to accuracy")
    ax.legend(fontsize=7)
    save(fig, "stage11_qa_cost_accuracy.png")

    # --- Q-B errors across a
    fig, ax = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    a = np.array(fv["a_points"])
    ax[0].plot(a, fv["dJda_ref"], "ks-", label="FV adjoint (N=2^18)")
    for k, rec in enumerate(pinn_b):
        if rec.get("failed"):
            continue
        ax[0].plot(
            a, rec["dJda"], "--o", ms=3, color=cols[k % 5], lw=1, label=f"seed {rec['seed']}"
        )
        ax[1].semilogy(a, rec["J_rel_err"], "-o", ms=3, color=cols[k % 5], lw=1)
        ax[2].semilogy(a, rec["dJda_rel_err"], "-o", ms=3, color=cols[k % 5], lw=1)
    ax[0].set_xlabel("a")
    ax[0].set_ylabel("dJ/da")
    ax[0].legend(fontsize=7)
    ax[0].set_title("dJ/da: surrogate autograd vs FV adjoint")
    ax[1].set_xlabel("a")
    ax[1].set_ylabel("rel. error of J")
    ax[1].set_title("surrogate J(a)")
    ax[2].set_xlabel("a")
    ax[2].set_ylabel("rel. error of dJ/da")
    ax[2].set_title("surrogate dJ/da")
    save(fig, "stage11_qb_surrogate_errors.png")

    # --- amortized cost
    fig, ax = plt.subplots(figsize=(6.5, 4.5), constrained_layout=True)
    ms = np.logspace(0, 6, 60)
    for k, s in enumerate(r2["per_seed"]):
        if s.get("failed"):
            continue
        if s["fv_seconds_per_eval_at_matched_N"] is not None:
            ax.loglog(
                ms,
                ms * s["fv_seconds_per_eval_at_matched_N"],
                color=cols[k % 5],
                lw=1,
                label=f"FV (matched to seed {s['seed']}, N={s['fv_matched_N']})",
            )
        ax.loglog(
            ms,
            s["pinn_train_seconds_scaled"] + ms * s["pinn_seconds_per_query"],
            "--",
            color=cols[k % 5],
            lw=1,
        )
    ax.axvline(1e4, color="0.7", lw=0.8)
    ax.set_xlabel("number of (J, dJ/da) evaluations M")
    ax.set_ylabel("total wall time  [s]")
    ax.set_title("Q-B amortized cost (solid FV, dashed PINN incl. training)")
    ax.legend(fontsize=6)
    save(fig, "stage11_qb_amortized_cost.png")


# ============================================================ main
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=REPO / "results" / "stage11")
    ap.add_argument("--seeds", nargs="+", type=int, default=SEEDS)
    ap.add_argument("--jobs", type=int, default=4, help="concurrent PINN subprocesses")
    ap.add_argument("--threads", type=int, default=None, help="torch threads per PINN process")
    ap.add_argument("--adam-steps", type=int, default=None, help="override (dev runs only)")
    ap.add_argument("--lbfgs-iter", type=int, default=None, help="override (dev runs only)")
    ap.add_argument("--single-pinn", nargs=2, metavar=("TASK", "SEED"), help=argparse.SUPPRESS)
    ap.add_argument("--single-fv", type=int, help=argparse.SUPPRESS)
    ap.add_argument("--calibrate", action="store_true", help=argparse.SUPPRESS)
    args = ap.parse_args()

    if args.single_pinn:
        print(json.dumps(single_pinn(args.single_pinn[0], int(args.single_pinn[1]), args)))
        return
    if args.single_fv:
        print(json.dumps(single_fv(args.single_fv)))
        return
    if args.calibrate:
        print(json.dumps(calibrate()))
        return

    from widegap.manifest import build_manifest

    t_start = time.perf_counter()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    overridden = args.adam_steps is not None or args.lbfgs_iter is not None
    jobs = max(1, args.jobs)
    threads = args.threads or max(1, _ncpu() // jobs)
    manifest = build_manifest(REPO, sys.argv, args.seeds)
    manifest.update(
        pre_registered_configuration=not overridden,
        pinn_jobs=jobs,
        pinn_threads_per_job=threads,
        fv_threads=_thread_info(),
    )
    _dump(out / "manifest.json", manifest)

    # ---- FV (isolated, sequential)
    t0 = time.perf_counter()
    fv = run_fv(out)
    for n in (1024, N_REF):
        p = subprocess.run(
            [sys.executable, str(SELF), "--single-fv", str(n)],
            capture_output=True,
            text=True,
            check=True,
        )
        rec = json.loads(p.stdout.strip().splitlines()[-1])
        fv.setdefault("subprocess_rss", {})[str(n)] = rec
        _dump(out / f"fv_rss_{n}.json", rec)
    fv["phase_seconds"] = time.perf_counter() - t0
    _dump(out / "fv.json", fv)
    print(
        f"FV done in {fv['phase_seconds']:.0f} s; M5 max rel err = "
        f"{fv['M5']['max_rel_err_vs_R2']:.2e}",
        flush=True,
    )

    # ---- isolated calibration of the Adam step time
    p = subprocess.run(
        [sys.executable, str(SELF), "--calibrate"], capture_output=True, text=True, check=True
    )
    cal = json.loads(p.stdout.strip().splitlines()[-1])
    _dump(out / "calibration.json", cal)

    # ---- PINN seeds, each in its own subprocess
    def run_one(task_seed):
        task, seed = task_seed
        cmd = [sys.executable, str(SELF), "--single-pinn", task, str(seed), "--out", str(out)]
        cmd += ["--threads", str(threads)]
        if args.adam_steps is not None:
            cmd += ["--adam-steps", str(args.adam_steps)]
        if args.lbfgs_iter is not None:
            cmd += ["--lbfgs-iter", str(args.lbfgs_iter)]
        t = time.perf_counter()
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            rec = {
                "task": task,
                "seed": seed,
                "failed": True,
                "error": proc.stderr.strip()[-2000:],
            }
        else:
            rec = json.loads(proc.stdout.strip().splitlines()[-1])
            if not rec.get("finite", True):
                rec["failed"] = True
                rec["error"] = "non-finite loss or parameters"
        rec["subprocess_wall_seconds"] = time.perf_counter() - t
        _dump(out / f"pinn_{task}_seed{seed}.json", rec)
        msg = (
            "FAILED"
            if rec.get("failed")
            else f"train {rec['train_seconds']:.0f} s, peak RSS {rec['peak_rss_mib']:.0f} MiB"
        )
        print(f"PINN {task} seed {seed}: {msg}", flush=True)
        return rec

    tasks = [("B", s) for s in args.seeds] + [("A", s) for s in args.seeds]
    with ThreadPoolExecutor(max_workers=jobs) as ex:
        recs = list(ex.map(run_one, tasks))
    pinn_b = [r for r in recs if r["task"] == "B"]
    pinn_a = [r for r in recs if r["task"] == "A"]

    # ---- load normalisation factors (1.0 when the run was sequential)
    factors = {}
    for r in recs:
        if r.get("failed") or jobs == 1:
            f = 1.0
        else:
            f = cal[f"{r['task']}_adam_step_s"] / r["adam_step_seconds"]
        factors[(r["task"], r["seed"])] = f
    fa = {s: factors[("A", s)] for s in args.seeds}
    fb = {s: factors[("B", s)] for s in args.seeds}

    r1 = assemble_r1(fv, pinn_a, fa)
    r2 = assemble_r2(fv, pinn_b, fb)
    r1["load_factors"], r2["load_factors"] = fa, fb
    r2["symmetric_degeneracy"] = fv["symmetric_degeneracy"]
    _dump(out / "r1.json", r1)
    _dump(out / "r2.json", r2)

    # ---- criteria
    m1a = fv["M1_A"]
    m1_ok = (
        m1a["psi_max_rel"] <= 1e-6
        and m1a["psi_l2_rel"] <= 1e-6
        and m1a["J_change_rel_to_int_abs_psi"] <= 1e-6
        and all(
            m["psi_max_rel"] <= 1e-6
            and m["psi_l2_rel"] <= 1e-6
            and m["J_change_rel_to_absJ"] <= 1e-6
            for m in fv["M1_B"]
        )
    )
    m3_vals = [fv["M3_A"]] + fv["M3_B"]
    m3_ok = fv["eval_points"] >= 4001 and all(
        max(m["interp_error_bound"], m["interp_error_empirical_est"]) <= 1e-7 for m in m3_vals
    )
    full = list(args.seeds) == SEEDS
    ok_recs = [r for r in recs if not r.get("failed")]
    m4_ok = (
        all(
            r.get("train_seconds") is not None
            and r.get("peak_rss_mib") is not None
            and r.get("threads", {}).get("torch_threads") is not None
            for r in ok_recs
        )
        and "1024" in fv["subprocess_rss"]
        and fv["timing_ref_mesh_forward_s"] > 0
    )
    crit = {
        "M1": m1_ok,
        "M2": full and len(pinn_a) == len(SEEDS) and len(pinn_b) == len(SEEDS),
        "M3": m3_ok,
        "M4": m4_ok,
        "M5": fv["M5"]["max_rel_err_vs_R2"] <= 1e-8,
    }
    crit = {k: bool(v) for k, v in crit.items()}
    failures = [
        {"task": r["task"], "seed": r["seed"], "error": r.get("error")}
        for r in recs
        if r.get("failed")
    ]
    concl = conclusion(fv, r1, r2)
    total = time.perf_counter() - t_start
    summary = {
        "criteria": crit,
        "gate_passed": all(crit.values()),
        "pre_registered_configuration": not overridden,
        "failures": failures,
        "R3_conclusion": concl,
        "total_runtime_seconds": total,
        "pinn_jobs": jobs,
        "pinn_threads_per_job": threads,
    }
    _dump(out / "summary.json", summary)

    try:
        figures(
            out,
            fv,
            r1,
            r2,
            pinn_a,
            pinn_b,
            f"stage 11 | {manifest['git_commit']} | seeds {args.seeds}",
            save_repo=(out == (REPO / "results" / "stage11").resolve()),
        )
    except Exception as e:  # figures must not hide the numerical results
        print("figure generation failed:", repr(e), flush=True)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
