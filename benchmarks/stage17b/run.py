"""Stage-17b benchmark (docs/benchmarks/stage17b_classification_constraint.md).

Stage 17 (benchmarks/stage17/run.py, imported unchanged) plus the classification constraints
c_k(theta) = (N_A,k - N_D,k(theta)) / N_A,k - m >= 0, m = 0.02, for every gated surface column that
is p-type at theta = 0 on BOTH the h/2 (on-state) and h/4 (re-check) meshes, with analytic gradients
(widegap.adjoint_channel.ClassificationConstraint).  Criteria: J1, J2, J3, J5, J6, J7 as in stage 17
(same thresholds, same job functions), J4b (sign guard at every SLSQP iterate / final / Pareto
designs on h/2 and h/4, clamp margins as J4; crossing line-search trial points are reported),
J8 (constraint-gradient FD <= 1e-8), R4 (active constraints, BV cost vs the stage-17 final,
re-classified R_on,ch of the stage-17 final).  No field-line / W-measure reports (not in 17b).

``--stage17-results`` points to the stage-17 results.json holding the stage-17 final theta.
``--dev`` marks the run DEV (refused when --out lies under results/).
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
from types import SimpleNamespace

import numpy as np
from scipy.optimize import minimize

from widegap.adjoint_channel import ClassificationConstraint, ron_reclassified
from widegap.manifest import build_manifest

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[2]
DEFAULT_OUT = REPO / "results" / "stage17b"
DEFAULT_S17 = REPO / "results" / "stage17" / "results.json"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


S17 = _load("stage17_run", REPO / "benchmarks/stage17/run.py")  # fixed name: fork pickling
M13D, M = S17.M13D, S17.M
DESIGN, TH0, Bounds, MEASURES = S17.DESIGN, S17.TH0, S17.Bounds, S17.MEASURES
tracked, job_diag = S17.tracked, S17.job_diag
MARGIN = 0.02
J8_TOL = 1e-8


def _lst(a):
    return np.asarray(a, dtype=float).tolist()


def make_cons():
    evs = [M13D.make_ev13d(S17.LEVEL_ON, with_cuts=False),
           M13D.make_ev13d(S17.LEVEL_RECHECK, with_cuts=False)]  # fmt: skip
    return ClassificationConstraint(evs, MARGIN)


def sign_status(cons, th):
    """(min value + margin, number of crossing columns) of the sign guard at theta."""
    c = cons.values(th) + MARGIN
    return float(c.min()), int(np.sum(c <= 0.0))


# ============================================================ SLSQP with the extra constraints
def run_slsqp_b(comp, cons, th_start, bv0, r0, e0, ratio, bounds, maxiter, ftol=1e-10):
    pr = M13D.Problem(comp, bv0, r0, e0, ratio, bounds)
    iterates: list = []
    cons_list = [
        {"type": "ineq", "fun": pr.c_ron, "jac": pr.c_ron_jac},
        {"type": "ineq", "fun": pr.c_ox, "jac": pr.c_ox_jac},
        {"type": "ineq", "fun": lambda t: bounds.A @ t + bounds.b, "jac": lambda t: bounds.A},
        {"type": "ineq", "fun": cons.values, "jac": cons.jac},
    ]

    def cb(xk):
        xk = np.asarray(xk, float)
        m, n = sign_status(cons, xk)
        iterates.append({"theta": xk.tolist(), "min_margin_plus_m": m, "n_crossing": n})

    try:
        res = minimize(
            pr.fun, np.asarray(th_start, float), jac=pr.jac, method="SLSQP", bounds=bounds.box,
            constraints=cons_list, options={"maxiter": maxiter, "ftol": ftol}, callback=cb,
        )  # fmt: skip
        return np.asarray(res.x, float), res, pr, iterates
    except M.EvalAbort as exc:
        best, bv_best = np.asarray(th_start, float), -np.inf
        for h in pr.history:
            t = np.asarray(h["theta"])
            if (h["R"] <= ratio * r0 * (1.0 + M.FEAS_TOL) and h["E_ox"] <= e0 * (1.0 + M.FEAS_TOL)
                    and bounds.violation(t) <= 1e-9 and sign_status(cons, t)[1] == 0
                    and h["bv_obj"] > bv_best):  # fmt: skip
                best, bv_best = t, h["bv_obj"]
        res = SimpleNamespace(status=-1, message=f"ABORTED (forward solve failed): {exc}",
                              nit=len(iterates), nfev=pr.n_eval, njev=pr.n_eval, x=best)  # fmt: skip
        return best, res, pr, iterates


# ============================================================ jobs
@tracked
def job_pareto17b(args):
    th_start, ratio, bv0, r0, e0, maxiter = args
    t0 = time.perf_counter()
    ev_off = S17.make_ev_off()
    cr = S17.make_cr(S17.LEVEL_ON)
    cons = make_cons()
    comp = S17.Eval17(ev_off, cr)
    bounds = Bounds(DESIGN, ev_off.sc)
    th, res, pr, its = run_slsqp_b(comp, cons, th_start, bv0, r0, e0, ratio, bounds,
                                   min(40, maxiter))  # fmt: skip
    b = ev_off.objective(th, grad=False)
    f = ev_off.fields(th, grad=False)
    mn, nc = sign_status(cons, th)
    return {
        "ratio": ratio, "theta": th.tolist(), "bv_obj": b["BV"], "bv_k": b["bv_k"].tolist(),
        "R": cr.ron(th, grad=False)["R"], "E_ox": f["E_ox"]["S"], "status": int(res.status),
        "message": res.message, "nit": int(res.nit), "bound_violation": bounds.violation(th),
        "iterates": its, "final_sign": {"min_margin_plus_m": mn, "n_crossing": nc},
        "n_trial_evals_crossing": int(sum(g["sign_changes"] > 0 for g in cr.guards)),
        "failure": pr.failure, "n_evals": pr.n_eval, "diag": job_diag([cr], [ev_off]),
        "seconds": time.perf_counter() - t0,
    }  # fmt: skip


@tracked
def job_reclass(args):
    th, level, tag = args
    t0 = time.perf_counter()
    cr = S17.make_cr(level)
    out = {"tag": tag, "level": level, "R_frozen": cr.ron(np.asarray(th, float), grad=False)["R"],
           **{"reclassified": ron_reclassified(cr.ev, cr.P, np.asarray(th, float))}}  # fmt: skip
    out["diag"] = job_diag([cr])
    out["seconds"] = time.perf_counter() - t0
    return out


def j8_check(cons, th):
    """Analytic constraint Jacobian vs Richardson central FD (eps = 1e-4, 5e-5)."""
    th = np.asarray(th, float)
    jac = cons.jac(th)
    fd = np.zeros_like(jac)
    for k in range(5):
        u = np.zeros(5)
        u[k] = 1.0

        def f1(e, u=u):
            return (cons.values(th + e * u) - cons.values(th - e * u)) / (2.0 * e)

        fd[:, k] = (4.0 * f1(5e-5) - f1(1e-4)) / 3.0
    return {"rel_err_fro": float(np.linalg.norm(jac - fd) / np.linalg.norm(fd)),
            "max_abs_err": float(np.max(np.abs(jac - fd))), "n_constraints": int(jac.shape[0])}  # fmt: skip


# ============================================================ figures
def plot_all(figdir, res, label):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figdir = Path(figdir)
    figdir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8.2, 5.0), constrained_layout=True)
    for key, c, name in (("initial_L0", "#0072B2", "initial"), ("final_L0", "#D55E00", "17b"),
                         ("s17_L0", "#009E73", "stage-17 final")):  # fmt: skip
        d = res["dense"].get(key)
        if d:
            ax.plot(d["x_um"], d["bv"], color=c, label=f"{name}: min {d['bv_min']:.0f} V")
    ax.set_xlabel("x [um]")
    ax.set_ylabel("per-cut BV on h [V]")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.text(0.5, 0.005, label, ha="center", fontsize=7, color="0.4")
    fig.savefig(figdir / "stage17b_per_cut_bv.png", dpi=130)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(7.2, 4.0), constrained_layout=True)
    x = np.arange(5)
    ax.bar(x - 0.2, res["R4"]["theta_17_final"], 0.4, label="stage 17")
    ax.bar(x + 0.2, res["final"]["theta"], 0.4, label="stage 17b")
    ax.set_xticks(x, [f"theta_{i}" for i in range(4)] + ["theta_J"])
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")
    fig.text(0.5, 0.005, label, ha="center", fontsize=7, color="0.4")
    fig.savefig(figdir / "stage17b_theta.png", dpi=130)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(7.2, 4.0), constrained_layout=True)
    par = res["R4"]["pareto_17b"]
    ax.plot([q["R_over_R0"] for q in par], [q["bv_true_min"] for q in par], "o-", label="17b")
    ax.set_xlabel("R_on,ch / R_on,ch,0 (h/2)")
    ax.set_ylabel("dense BV on h [V]")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.text(0.5, 0.005, label, ha="center", fontsize=7, color="0.4")
    fig.savefig(figdir / "stage17b_pareto.png", dpi=130)
    plt.close(fig)


# ============================================================ main
def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--stage17-results", type=Path, default=DEFAULT_S17)
    ap.add_argument("--maxiter", type=int, default=M.MAX_ITER)
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
    man = build_manifest(REPO, sys.argv, list(M.SEEDS))
    man["dev"] = bool(args.dev)
    man["stage17_results"] = str(args.stage17_results)
    (out / "manifest.json").write_text(json.dumps(man, indent=2) + "\n")
    t_start = time.perf_counter()
    th17 = np.asarray(json.loads(args.stage17_results.read_text())["final"]["theta"], float)
    ev_off = S17.make_ev_off()
    cr = S17.make_cr(S17.LEVEL_ON)
    cons = make_cons()
    comp = S17.Eval17(ev_off, cr)
    bounds = Bounds(DESIGN, ev_off.sc)
    res: dict = {"dev": bool(args.dev)}
    diags: dict = {}
    T: dict = {}

    def dump():
        (out / "results.json").write_text(json.dumps(res, indent=2, default=float) + "\n")

    res["config"] = {
        "margin_m": MARGIN, "n_constraints": cons.n, "constraint_columns": cons.labels(),
        "max_iter": args.maxiter, "workers": args.workers, "stage17_final_theta": th17.tolist(),
        "tolerances": {"J1": S17.J1_TOL, "J2": S17.J2_TOL, "J5_ratio": S17.J5_RATIO,
                       "J6": S17.J6_TOL, "J8": J8_TOL},
        "registry": S17.registry_check(),
    }  # fmt: skip
    o0 = ev_off.objective(TH0, grad=False)
    f0 = ev_off.fields(TH0, grad=False)
    r0o = cr.ron(TH0, grad=False)
    bv0, e0, r0 = o0["BV"], f0["E_ox"]["S"], r0o["R"]
    c0 = cons.values(TH0)
    res["initial"] = {"bv_obj": bv0, "R_channel_h2": r0, "E_ox": e0,
                      "min_constraint_value": float(c0.min())}  # fmt: skip
    print(f"initial: BV_obj {bv0:.3f}, R_ch {r0:.6e}, E_ox {e0:.4f}, min c {c0.min():.4f}",
          flush=True)  # fmt: skip

    pool = ProcessPoolExecutor(args.workers, mp_context=mp.get_context("fork"))
    fut: dict = {}

    def submit(key, fn, a):
        fut[pool.submit(fn, a)] = key

    for k in range(5):
        submit(("j1_fd", k), S17.job_fd_ron, (TH0, k, S17.LEVEL_ON, "initial"))
    submit(("gc0", None), S17.job_grad_off, (TH0, e0, "initial"))
    submit(("ron_pt", "initial_L1"), S17.job_ron_point, (TH0, S17.LEVEL_ON, "initial", True))
    submit(("ron_pt", "initial_L2"), S17.job_ron_point, (TH0, S17.LEVEL_RECHECK, "initial", False))
    for lev in (0, 1):
        submit(("dense", f"initial_L{lev}"), S17.job_dense17, (lev, TH0, "initial", False))

    t0 = time.perf_counter()
    th1, sl, pr, its = run_slsqp_b(comp, cons, TH0, bv0, r0, e0, 1.0, bounds, args.maxiter)
    T["slsqp_main"] = time.perf_counter() - t0
    o1 = ev_off.objective(th1, grad=False)
    ro1 = cr.ron(th1, grad=False)
    f1 = ev_off.fields(th1, grad=False)
    bv1, r1, e1 = o1["BV"], ro1["R"], f1["E_ox"]["S"]
    c1 = cons.values(th1)
    print(f"main: BV_obj {bv0:.3f} -> {bv1:.3f}, R_ch {r0:.5e} -> {r1:.5e}, "
          f"E_ox {e0:.4f} -> {e1:.4f}, min c {c1.min():.4g}, {sl.message}", flush=True)  # fmt: skip
    n_trial_cross = int(sum(g["sign_changes"] > 0 for g in cr.guards))
    j8 = {"initial": j8_check(cons, TH0), "final": j8_check(cons, th1)}
    res["final"] = {
        "theta": th1.tolist(), "bv_obj": bv1, "R_channel_h2": r1, "E_ox": e1,
        "R_over_R0_minus_1": r1 / r0 - 1.0, "E_ox_over_E0_minus_1": e1 / e0 - 1.0,
        "bound_violation": bounds.violation(th1), "constraint_values": _lst(c1),
    }  # fmt: skip
    res["slsqp_main"] = {"status": int(sl.status), "message": sl.message, "nit": int(sl.nit),
                         "iterates": its, "evaluations": pr.history,
                         "n_trial_evals_crossing": n_trial_cross}  # fmt: skip
    diags["main"] = job_diag([cr], [ev_off])
    dump()

    for k in range(5):
        submit(("j2_fd", k), S17.job_fd_ron, (th1, k, S17.LEVEL_ON, "final"))
    submit(("gc1", None), S17.job_grad_off, (th1, e0, "final"))
    submit(("ron_pt", "final_L1"), S17.job_ron_point, (th1, S17.LEVEL_ON, "final", True))
    submit(("ron_pt", "final_L2"), S17.job_ron_point, (th1, S17.LEVEL_RECHECK, "final", False))
    for rt in S17.PARETO_RATIOS:
        submit(("pareto", rt), job_pareto17b, (th1, rt, bv0, r0, e0, args.maxiter))
    for lev in (0, 1):
        submit(("dense", f"final_L{lev}"), S17.job_dense17, (lev, th1, "final", False))
        submit(("dense", f"s17_L{lev}"), S17.job_dense17, (lev, th17, "s17", False))
    submit(("reclass", 1), job_reclass, (th17, S17.LEVEL_ON, "s17_final"))
    submit(("reclass", 2), job_reclass, (th17, S17.LEVEL_RECHECK, "s17_final"))
    jobs: dict = {k: {} for k in ("j1_fd", "j2_fd", "pareto", "dense", "ron_pt", "reclass")}
    jobs["gc0"] = jobs["gc1"] = None
    pending = set(fut)
    t0 = time.perf_counter()
    while pending:
        for f in as_completed(list(pending)):
            pending.discard(f)
            kind, key = fut[f]
            try:
                r = f.result()
            except Exception as exc:
                r = {"error": f"{type(exc).__name__}: {exc}", "seconds": float("nan")}
            if kind in ("gc0", "gc1"):
                jobs[kind] = r
            else:
                jobs[kind][key] = r
            if "diag" in r:
                diags[f"{kind}|{key}"] = r["diag"]
            print(f"done {kind} {key} {r.get('seconds', float('nan')):.0f}s"
                  f"{' ERROR ' + r['error'] if 'error' in r else ''}", flush=True)  # fmt: skip
            if kind == "pareto" and "error" not in r:
                nf = pool.submit(
                    S17.job_dense17, (0, np.asarray(r["theta"]), f"ratio_{key:g}", False)
                )
                fut[nf] = ("dense", f"ratio_{key:g}_L0")
                pending.add(nf)
                nf2 = pool.submit(S17.job_ron_point, (np.asarray(r["theta"]), S17.LEVEL_RECHECK,
                                                      f"ratio_{key:g}", False))  # fmt: skip
                fut[nf2] = ("ron_pt", f"ratio_{key:g}_L2")
                pending.add(nf2)
            break
    pool.shutdown()
    T["phase2_parallel"] = time.perf_counter() - t0
    dense = jobs["dense"]
    res["dense"], res["jobs"], res["pareto_runs"] = dense, jobs, jobs["pareto"]
    dump()

    # ================================================================ criteria
    crit, vals = {}, {}
    g0 = jobs["ron_pt"]["initial_L1"]
    fdv = np.array([jobs["j1_fd"][k]["richardson"] for k in range(5)])
    j1_err = float(np.linalg.norm(np.array(g0["grad"]) - fdv) / np.linalg.norm(fdv))
    gc0, gf0 = jobs["gc0"]["grad_check"], jobs["gc0"]["grad_check_fields"]
    d2_bv = max(r["bv"]["rel_err"] for r in gc0["directions"])
    d2_e = max(r["E_ox"]["rel_err"] for r in gf0["directions"])
    crit["J1"] = bool(j1_err <= S17.J1_TOL and d2_bv <= S17.D2_TOL_REL and d2_e <= S17.D2_TOL_REL)
    vals["J1"] = {"rel_err_norm_R_on_ch": j1_err, "worst_rel_BV_obj": d2_bv, "worst_rel_E_ox": d2_e}
    gc1, gf1 = jobs["gc1"]["grad_check"], jobs["gc1"]["grad_check_fields"]
    g1 = np.array(jobs["ron_pt"]["final_L1"]["grad"])
    fd1 = np.array([jobs["j2_fd"][k]["richardson"] for k in range(5)])
    a_r = float(np.max(np.abs(g1 - fd1)) / r0)
    n_r, n_bv, n_e = (float(np.linalg.norm(g1) / r0), gc1["grad_bv_norm"] / bv0,
                      gf1["grad_E_over_E0_norm"])  # fmt: skip
    a_bv = max(r["bv"]["abs_diff"] for r in gc1["directions"]) / bv0
    a_e = max(r["E_ox"]["abs_diff"] for r in gf1["directions"]) / e0
    tol = {k: S17.J2_TOL * max(1.0, n) for k, n in (("bv", n_bv), ("r", n_r), ("e", n_e))}
    crit["J2"] = bool(a_bv <= tol["bv"] and a_r <= tol["r"] and a_e <= tol["e"])
    vals["J2"] = {"abs_BV_over_BV0": a_bv, "abs_R_over_R0": a_r, "abs_E_over_E0": a_e, "tol": tol,
                  "grad_norms_normalized": [n_bv, n_r, n_e]}  # fmt: skip
    viol = bounds.violation(th1)
    r1_h2 = jobs["ron_pt"]["final_L1"]["R"]
    crit["J3"] = bool(r1 <= r0 * (1 + 1e-8) and r1_h2 <= r0 * (1 + 1e-8)
                      and e1 <= e0 * (1 + 1e-8) and viol <= 1e-9)  # fmt: skip
    vals["J3"] = {"R_over_R0_minus_1": r1 / r0 - 1.0, "E_over_E0_minus_1": e1 / e0 - 1.0,
                  "bound_violation": viol}  # fmt: skip
    # J4b
    pr_runs = jobs["pareto"]
    iters_main = its
    n_it = len(iters_main) + sum(len(q["iterates"]) for q in pr_runs.values())
    cross_it = sum(i["n_crossing"] > 0 for i in iters_main) + sum(
        i["n_crossing"] > 0 for q in pr_runs.values() for i in q["iterates"])  # fmt: skip
    cross_fin = int(sign_status(cons, th1)[1] > 0) + sum(
        q["final_sign"]["n_crossing"] > 0 for q in pr_runs.values())  # fmt: skip
    g_all = [d["guards"] for d in diags.values() if "guards" in d and d["guards"]["n"]]
    fin_guards = [jobs["ron_pt"][k]["guard"] for k in jobs["ron_pt"] if k.endswith("_L1")
                  or k.endswith("_L2")]  # fmt: skip
    fin_clamp = sum(g["gated_near_ns_clamp"] + g["gated_near_eeff_clamp"] for g in fin_guards)
    n_trial_total = n_trial_cross + sum(q["n_trial_evals_crossing"] for q in pr_runs.values())
    crit["J4b"] = bool(
        cross_it == 0
        and cross_fin == 0
        and fin_clamp == 0
        and bool(g_all)
        and min(g["min_gated_ns_cm2"] for g in g_all) > 2e11
        and min(g["min_gated_eeff_raw_MV_cm"] for g in g_all) > 2e-2
    )
    vals["J4b"] = {
        "n_iterates_checked": n_it, "n_iterates_crossing": cross_it,
        "n_final_designs_crossing": cross_fin, "clamp_margin_violations_at_designs": fin_clamp,
        "n_trial_evaluations_crossing_reported": n_trial_total,
        "min_gated_ns_cm2": min(g["min_gated_ns_cm2"] for g in g_all),
        "min_gated_eeff_raw_MV_cm": min(g["min_gated_eeff_raw_MV_cm"] for g in g_all),
        "min_constraint_value_final": float(c1.min()),
        "note": "sign guard checked on h/2 and h/4 columns at every iterate (callback), final and "
        "Pareto designs; evaluations at line-search trial points that cross are counted, not gated",
    }  # fmt: skip
    crit["J8"] = bool(max(j8["initial"]["rel_err_fro"], j8["final"]["rel_err_fro"]) <= J8_TOL)
    vals["J8"] = j8
    ini = {lev: dense[f"initial_L{lev}"] for lev in (0, 1)}
    fin = {lev: dense[f"final_L{lev}"] for lev in (0, 1)}
    imp = {lev: fin[lev]["bv_min"] - ini[lev]["bv_min"] for lev in (0, 1)}
    crit["J5"] = bool(imp[0] > 0 and imp[1] > 0 and imp[1] >= S17.J5_RATIO * imp[0]
                      and len(ini[0]["x_um"]) >= 40)  # fmt: skip
    vals["J5"] = {"improvement_h": imp[0], "improvement_h2": imp[1],
                  "min_bv_initial_h_h2": [ini[0]["bv_min"], ini[1]["bv_min"]],
                  "min_bv_final_h_h2": [fin[0]["bv_min"], fin[1]["bv_min"]]}  # fmt: skip
    r_i4, r_f4 = jobs["ron_pt"]["initial_L2"]["R"], jobs["ron_pt"]["final_L2"]["R"]
    crit["J6"] = bool(r_f4 <= r_i4 * (1.0 + S17.J6_TOL))
    vals["J6"] = {"R_h4_initial": r_i4, "R_h4_final": r_f4, "rel": r_f4 / r_i4 - 1.0}
    # J7
    nf_tot = sum(d["nonfinite_calls"] for d in diags.values())
    nf_failed = sum(d["NF"]["failed_newton"] for d in diags.values())
    nf_adj = sum(d["NF"]["adjoint"] for d in diags.values())
    nf_conv = sum(d["NF"]["converged_newton"] for d in diags.values())

    def agg(a, b):
        return sum(d[a].get(b, 0) for d in diags.values() if a in d)

    stderr = log_path.read_text(errors="replace") if log_path.exists() else ""
    entry = [ln for ln in stderr.splitlines() if "On entry to" in ln]
    errs = [f"{kd}|{k}" for kd in ("j1_fd", "j2_fd", "pareto", "dense", "ron_pt", "reclass")
            for k, v in jobs[kd].items() if "error" in v]  # fmt: skip
    errs += [k for k in ("gc0", "gc1") if jobs[k] is None or "error" in jobs[k]]
    crit["J7"] = bool(
        sl.status == 0 and sl.nit <= args.maxiter and nf_tot == nf_failed and nf_adj == 0
        and nf_conv == 0 and agg("on_state", "failed_picard") == 0
        and agg("on_state", "failed_evals") == 0 and agg("off_state", "failed_evals") == 0
        and not entry and not errs and all(d["resolved_all"] for d in dense.values())
    )  # fmt: skip
    vals["J7"] = {"slsqp_status": int(sl.status), "slsqp_message": sl.message, "nit": int(sl.nit),
                  "pareto_status": {str(k): v["status"] for k, v in pr_runs.items()},
                  "nonfinite_calls_total": nf_tot, "nonfinite_in_failed_newton": nf_failed,
                  "nonfinite_in_adjoint": nf_adj, "nonfinite_in_converged_newton": nf_conv,
                  "failed_picard": agg("on_state", "failed_picard"),
                  "failed_evals_on": agg("on_state", "failed_evals"),
                  "failed_evals_off": agg("off_state", "failed_evals"),
                  "on_entry_to_lines": len(entry), "job_errors": errs}  # fmt: skip
    # ================================================================ R4
    lab = cons.labels()
    order = np.argsort(c1)
    active = [{**lab[i], "value": float(c1[i])} for i in order if c1[i] < 1e-6]
    near = [{**lab[i], "value": float(c1[i])} for i in order[:6]]
    par = [{"ratio": 1.0, "theta": th1.tolist(), "bv_obj": bv1, "R": r1, "E_ox": e1,
            "bv_true_min": fin[0]["bv_min"]}]  # fmt: skip
    for rt in S17.PARETO_RATIOS:
        q = dict(pr_runs[rt])
        q["bv_true_min"] = dense[f"ratio_{rt:g}_L0"]["bv_min"]
        par.append(q)
    for q in par:
        q["R_over_R0"] = q["R"] / r0
    rc = jobs["reclass"]
    res["R4"] = {
        "active_constraints": active, "nearest_constraints": near,
        "active_R": bool(abs(r1 / r0 - 1) < 1e-6), "active_E_ox": bool(abs(e1 / e0 - 1) < 1e-6),
        "theta_17b": th1.tolist(), "theta_17_final": th17.tolist(),
        "bv_true_min_h": {"17b": fin[0]["bv_min"], "17": dense["s17_L0"]["bv_min"],
                          "initial": ini[0]["bv_min"]},
        "bv_true_min_h2": {"17b": fin[1]["bv_min"], "17": dense["s17_L1"]["bv_min"],
                           "initial": ini[1]["bv_min"]},
        "bv_cost_of_constraint_h": fin[0]["bv_min"] - dense["s17_L0"]["bv_min"],
        "bv_obj_17b": bv1,
        "reclassified_R_s17_final": {str(k): v.get("reclassified") for k, v in rc.items()},
        "frozen_R_s17_final": {str(k): v.get("R_frozen") for k, v in rc.items()},
        "reclass_rel_change_h2": (rc[1]["reclassified"]["R"] / rc[1]["R_frozen"] - 1.0)
        if 1 in rc and "reclassified" in rc[1] else None,
        "reclass_vs_R0_h2": (rc[1]["reclassified"]["R"] / r0 - 1.0)
        if 1 in rc and "reclassified" in rc[1] else None,
        "pareto_17b": par,
    }  # fmt: skip
    res["diagnostics"], res["seconds"] = diags, T
    res["runtime_seconds"] = time.perf_counter() - t_start
    figdir = (REPO / "figures") if out == DEFAULT_OUT.resolve() else out / "figures"
    label = f"stage 17b | commit {str(man.get('git_commit'))[:8]} | {man['timestamp_utc'][:19]}Z"
    try:
        plot_all(figdir, res, label)
    except Exception as exc:
        res["figures_error"] = f"{type(exc).__name__}: {exc}"
    summary = {"criteria": crit, "gate_passed": bool(all(crit.values())), "dev": bool(args.dev),
               "values": {"bv_obj_initial": bv0, "bv_obj_final": bv1, "R_ch_initial": r0,
                          "R_ch_final": r1, **vals, "R4": res["R4"]},
               "seconds": T, "runtime_seconds": res["runtime_seconds"]}  # fmt: skip
    dump()
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    print(json.dumps(crit, indent=2))
    print(f"runtime {res['runtime_seconds']:.1f} s")


if __name__ == "__main__":
    main()
