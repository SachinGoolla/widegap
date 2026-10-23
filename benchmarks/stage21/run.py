"""Stage-21 benchmark (docs/benchmarks/stage21_safe_geometry_opt.md): S1-S6 (gating), R1-R2.

The stage-17b problem (channel-model R_on on h/2, E_ox at 400 V on h, 17b classification
constraints on h/2 and h/4, soft-min BV on h) on the stage-20 punch-through-safe geometry
d_body = 1.4 um (``widegap.bodyredesign.body_design`` and its mesh rule), with the punch-through
guard V_pt,es(theta) >= 1.2 BV checked (NOT a constraint) at every SLSQP iterate, the final design
and every Pareto design on h, and at the final design on h/2.

Geometry mapping (recorded in results.json["config"]): the device only gets deeper (the drift starts
below the body), x is unchanged, so the 13d cut x positions, the 49 dense lines and the x extent are
used as is.  The 13d diagnostic windows W1/W2 (y boxes around the 0.8 um junction) are shifted
down by 0.6 um (W3, the JFET surface, is unchanged); they are reported nowhere in stage 21.

All stage-13d / 13b / 15b / 17 / 17b code is reused by import; the geometry is installed by
replacing the ``DESIGN`` globals of those modules (and ``make_ev13d``) before any worker forks.
Guard BV: BV_obj at iterates (the problem text), max(BV_obj, dense true-min BV) at the final and
Pareto designs.  V_pt,es: smallest V_D with the minimum n+ -> body barrier < 10 V_T (stage 20).

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

from widegap.adjoint_channel import D_CH_UM
from widegap.adjoint_sic_oxide import OxideAdjoint
from widegap.bodyredesign import OffMarch, body_design, column_index_matrix, source_columns
from widegap.channel_mobility import ChannelMap, acc_mobility_from_registry, picard_solve
from widegap.manifest import build_manifest

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[2]
DEFAULT_OUT = REPO / "results" / "stage21"
DEFAULT_S17B = REPO / "results" / "stage17b" / "results.json"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


# stage 20 first (it re-registers the 13d chain); stage 17b/17 last so the final registration is
# the one whose functions are pickled to the fork workers
S20 = _load("stage20_run", REPO / "benchmarks/stage20/run.py")
S17B = _load("stage17b_run", REPO / "benchmarks/stage17b/run.py")  # fixed names: fork pickling
S17 = S17B.S17
M13D, M, S15B = S17.M13D, S17.M, S17.S15B
tracked, job_diag = S17.tracked, S17.job_diag

D_BODY_UM = 1.4
GUARD_RATIO = 1.2
BARRIER_LEVEL_VT = 10.0
V_CEIL = 1200.0
S5_TOL = 0.03
Y_SHIFT_UM = D_BODY_UM - 0.8
DESIGN0 = M13D.DESIGN
DESIGN = body_design(DESIGN0, D_BODY_UM)
WINDOWS_UM = {
    "W1": ((0.5, 1.1), (-1.1 - Y_SHIFT_UM, -0.5 - Y_SHIFT_UM)),
    "W2": ((1.3, 2.0), (-1.2 - Y_SHIFT_UM, -0.6 - Y_SHIFT_UM)),
    "W3": ((0.0, 0.4), (-0.2, 0.0)),
}


def make_ev13d_deep(level, with_cuts=True):
    sc = M13D.SiC4H(M13D.TracedRegistry.from_json())
    return OxideAdjoint(sc, DESIGN, level, r_um=M.R_UM, cuts_um=M.CUTS_UM if with_cuts else (),
                        tau_v=M.TAU_V, windows_um=WINDOWS_UM)  # fmt: skip


def install_geometry():
    """Replace the DESIGN globals of every reused module and the 13d evaluator factory."""
    for mod in (M, M13D, S17, S17B, S15B, S15B.M):
        mod.DESIGN = DESIGN
    M13D.make_ev13d = make_ev13d_deep


install_geometry()
TH0, Bounds, MEASURES = S17.TH0, S17.Bounds, S17.MEASURES


def _lst(a):
    return np.asarray(a, dtype=float).tolist()


# ============================================================ punch-through guard (a check)
class PtGuard:
    """V_pt,es(theta) on one mesh level; columns frozen at theta = 0 (stage 20)."""

    def __init__(self, level):
        self.ev = S20.make_ev(D_BODY_UM, level)
        self.ev.set_theta(TH0)
        self.cols = source_columns(self.ev)
        self.idx = column_index_matrix(self.ev, self.cols)
        self.level = level

    def vpt(self, th):
        try:
            r = OffMarch(self.ev, np.asarray(th, float)).vpt_es(
                self.idx, BARRIER_LEVEL_VT, v_max_v=V_CEIL
            )
        except Exception as exc:  # recorded as a guard failure
            return {"V_pt_es": None, "ceiling": False, "error": f"{type(exc).__name__}: {exc}"}
        v = r["V_pt_es"]
        return {"V_pt_es": float(v) if v is not None else V_CEIL, "ceiling": v is None,
                "error": None}  # fmt: skip


def guard_ok(vpt, bv):
    return bool(vpt.get("error") is None and vpt["V_pt_es"] >= GUARD_RATIO * bv)


# ============================================================ SLSQP with constraints + guard
def run_slsqp21(comp, cons, pt, th_start, bv0, r0, e0, ratio, bounds, maxiter, ftol=1e-10):
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
        m, n = S17B.sign_status(cons, xk)
        bv = comp.ev_off.objective(xk, grad=False)["BV"]
        v = pt.vpt(xk)
        iterates.append({"theta": xk.tolist(), "min_margin_plus_m": m, "n_crossing": n,
                         "bv_obj": float(bv), "V_pt_es": v["V_pt_es"], "ceiling": v["ceiling"],
                         "pt_error": v["error"], "ratio": None if v["error"] else v["V_pt_es"] / bv,
                         "guard_ok": guard_ok(v, bv)})  # fmt: skip

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
                    and bounds.violation(t) <= 1e-9 and S17B.sign_status(cons, t)[1] == 0
                    and h["bv_obj"] > bv_best):  # fmt: skip
                best, bv_best = t, h["bv_obj"]
        res = SimpleNamespace(status=-1, message=f"ABORTED (forward solve failed): {exc}",
                              nit=len(iterates), nfev=pr.n_eval, njev=pr.n_eval, x=best)  # fmt: skip
        return best, res, pr, iterates


# ============================================================ jobs
@tracked
def job_pareto21(args):
    th_start, ratio, bv0, r0, e0, maxiter = args
    t0 = time.perf_counter()
    ev_off = S17.make_ev_off()
    cr = S17.make_cr(S17.LEVEL_ON)
    cons = S17B.make_cons()
    pt = PtGuard(0)
    comp = S17.Eval17(ev_off, cr)
    bounds = Bounds(DESIGN, ev_off.sc)
    th, res, pr, its = run_slsqp21(comp, cons, pt, th_start, bv0, r0, e0, ratio, bounds,
                                   min(40, maxiter))  # fmt: skip
    b = ev_off.objective(th, grad=False)
    f = ev_off.fields(th, grad=False)
    mn, nc = S17B.sign_status(cons, th)
    v = pt.vpt(th)
    return {
        "ratio": ratio, "theta": th.tolist(), "bv_obj": b["BV"], "R": cr.ron(th, grad=False)["R"],
        "E_ox": f["E_ox"]["S"], "status": int(res.status), "message": res.message,
        "nit": int(res.nit), "bound_violation": bounds.violation(th), "iterates": its,
        "final_sign": {"min_margin_plus_m": mn, "n_crossing": nc},
        "final_vpt": v, "n_trial_evals_crossing": int(sum(g["sign_changes"] > 0 for g in cr.guards)),
        "failure": pr.failure, "n_evals": pr.n_eval, "diag": job_diag([cr], [ev_off]),
        "seconds": time.perf_counter() - t0,
    }  # fmt: skip


@tracked
def job_vpt21(args):
    th, level, tag = args
    t0 = time.perf_counter()
    pt = PtGuard(level)
    v = pt.vpt(th)
    if v["error"] is None:
        v.pop("error")  # an "error" key marks a failed job
    out = {"tag": tag, "level": level, **v, "diag": job_diag()}
    out["seconds"] = time.perf_counter() - t0
    return out


@tracked
def job_vt21(args):
    th, tag = args
    t0 = time.perf_counter()
    S20.TH0 = np.asarray(th, float)  # vt_es reads the module-level theta
    try:
        vt = S20.vt_es(D_BODY_UM, 1)
    finally:
        S20.TH0 = TH0
    return {"tag": tag, "V_T_es_V": float(vt), "diag": job_diag(),
            "seconds": time.perf_counter() - t0}  # fmt: skip


@tracked
def job_acc21(args):
    """18b L' (excess accumulation mobility, Zhao) and U (stage-16 channel) R_on at theta."""
    th, level, tag = args
    t0 = time.perf_counter()
    cr = S17.make_cr(level)
    ev, sc = cr.ev, cr.ev.sc
    th = np.asarray(th, float)
    ru = cr.ron(th, grad=False)["R"]
    mu_acc = acc_mobility_from_registry(sc.reg)
    ev.set_theta(th)
    sol_b, bc = ev.on_state()
    cm = ChannelMap(sc, ev.mod, ev.nd, ev.na, D_CH_UM,
                    gate_x_max=ev.prm.gate_x1, mu_acc_cm2=mu_acc, acc_mode="excess")  # fmt: skip
    ev.mod.mu_edge = cm.mu_bulk.copy()
    r = picard_solve(cm, cr.P, sol_b.X, bc)
    out = {"tag": tag, "level": level, "R_U": ru, "converged": bool(r["converged"]),
           "n_acc_edges": int(len(cm.acc_edges))}  # fmt: skip
    if r["converged"]:
        i_sc = ev.mod.contact_current(r["sol"], ev.base.dev.nodes["drain"])
        out["R_Lp"] = float(ev.V_DS * ev.x_cm / ev.base.current_a_per_cm(abs(i_sc)))
    ev.mod.mu_edge = cm.mu_bulk.copy()
    out["diag"] = job_diag([cr])
    out["seconds"] = time.perf_counter() - t0
    return out


@tracked
def job_ox21(args):
    th, tag = args
    t0 = time.perf_counter()
    out = {"tag": tag}
    for lv in (0, 1):
        ev = make_ev13d_deep(lv, with_cuts=False)
        out[f"E_ox_L{lv}"] = float(ev.fields(th, grad=False, which=("E_ox",))["E_ox"]["S"])
    out["seconds"] = time.perf_counter() - t0
    return out


# ============================================================ figures
def plot_all(figdir, res, label):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figdir = Path(figdir)
    figdir.mkdir(parents=True, exist_ok=True)
    foot = lambda fig: fig.text(0.5, 0.005, label, ha="center", fontsize=7, color="0.4")  # noqa: E731
    fig, ax = plt.subplots(figsize=(8.0, 4.6), constrained_layout=True)
    for key, c, nm in (("initial_L0", "#0072B2", "initial"), ("final_L0", "#D55E00", "stage 21")):
        d = res["dense"].get(key)
        if d:
            ax.plot(d["x_um"], d["bv"], color=c, label=f"{nm}: min {d['bv_min']:.0f} V")
    ax.set_xlabel("x [um]")
    ax.set_ylabel("per-cut BV on h [V]")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    foot(fig)
    fig.savefig(figdir / "stage21_per_cut_bv.png", dpi=130)
    plt.close(fig)
    its = res["slsqp_main"]["iterates"]
    fig, ax = plt.subplots(figsize=(7.2, 4.2), constrained_layout=True)
    ax.plot([i["bv_obj"] for i in its], [i["V_pt_es"] for i in its], "o-", label="SLSQP iterates")
    xs = np.array([i["bv_obj"] for i in its] + [res["initial"]["bv_obj"]])
    ax.plot(xs, GUARD_RATIO * xs, "k--", label="V_pt = 1.2 BV")
    ax.set_xlabel("BV_obj [V]")
    ax.set_ylabel("V_pt,es [V]")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    foot(fig)
    fig.savefig(figdir / "stage21_pt_guard.png", dpi=130)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(7.2, 4.0), constrained_layout=True)
    par = res["R1"]["pareto"]
    ax.plot([q["R_over_R0"] for q in par], [q["bv_true_min"] for q in par], "o-")
    ax.set_xlabel("R_on,ch / R_on,ch,0 (h/2)")
    ax.set_ylabel("dense BV on h [V]")
    ax.grid(alpha=0.3)
    foot(fig)
    fig.savefig(figdir / "stage21_pareto.png", dpi=130)
    plt.close(fig)


# ============================================================ main
def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--stage17b-results", type=Path, default=DEFAULT_S17B)
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
    (out / "manifest.json").write_text(json.dumps(man, indent=2) + "\n")
    t_start = time.perf_counter()
    ev_off = S17.make_ev_off()
    sc = ev_off.sc
    cr = S17.make_cr(S17.LEVEL_ON)
    cons = S17B.make_cons()
    pt = PtGuard(0)
    comp = S17.Eval17(ev_off, cr)
    bounds = Bounds(DESIGN, sc)
    res: dict = {"dev": bool(args.dev)}
    diags: dict = {}
    T: dict = {}

    def dump():
        (out / "results.json").write_text(json.dumps(res, indent=2, default=float) + "\n")

    res["config"] = {
        "d_body_um": D_BODY_UM, "ny_body": DESIGN.mesh_kw["ny_body"],
        "mesh_nodes_off_L0": int(ev_off.mesh.n_nodes), "mesh_nodes_on_L1": int(cr.ev.mesh.n_nodes),
        "cut_x_um": ev_off.cut_x_um, "windows_um_shifted": {k: [list(a), list(b)] for k, (a, b)
                                                              in WINDOWS_UM.items()},
        "n_classification_constraints": cons.n, "constraint_columns": cons.labels(),
        "guard": {"ratio": GUARD_RATIO, "level_VT": BARRIER_LEVEL_VT, "ceiling_V": V_CEIL,
                  "columns_L0_um": [float(sc.length_cm(ev_off.mesh.x[i])) * 1e4
                                    for i in pt.cols]},
        "max_iter": args.maxiter, "workers": args.workers, "registry": S17.registry_check(),
    }  # fmt: skip
    o0 = ev_off.objective(TH0, grad=False)
    f0 = ev_off.fields(TH0, grad=False)
    r0o = cr.ron(TH0, grad=False)
    bv0, e0, r0 = o0["BV"], f0["E_ox"]["S"], r0o["R"]
    v0 = pt.vpt(TH0)
    res["initial"] = {"bv_obj": bv0, "R_channel_h2": r0, "E_ox": e0, "V_pt_es_h": v0,
                      "guard_ok": guard_ok(v0, bv0), "min_constraint_value": float(cons.values(TH0).min()),
                      "cut_limiting_x_um": ev_off.cut_x_um[o0["limiting_cut"]]}  # fmt: skip
    print(f"initial: BV_obj {bv0:.3f}, R_ch {r0:.6e}, E_ox {e0:.4f}, V_pt {v0['V_pt_es']:.1f}",
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
        submit(("fl", f"initial_L{lev}"), S17.job_fl17,
               (lev, TH0, "initial", S17.N_S_FL, S17.DS_FL_UM))  # fmt: skip
        submit(("vpt", f"initial_L{lev}"), job_vpt21, (TH0, lev, "initial"))
    submit(("vt", "initial"), job_vt21, (TH0, "initial"))
    submit(("ox", "initial"), job_ox21, (TH0, "initial"))
    for lv in (1, 2):
        submit(("acc", f"initial_L{lv}"), job_acc21, (TH0, lv, "initial"))

    t0 = time.perf_counter()
    th1, sl, pr, its = run_slsqp21(comp, cons, pt, TH0, bv0, r0, e0, 1.0, bounds, args.maxiter)
    T["slsqp_main"] = time.perf_counter() - t0
    o1 = ev_off.objective(th1, grad=False)
    ro1 = cr.ron(th1, grad=False)
    f1 = ev_off.fields(th1, grad=False)
    bv1, r1, e1 = o1["BV"], ro1["R"], f1["E_ox"]["S"]
    c1 = cons.values(th1)
    v1 = pt.vpt(th1)
    print(f"main: BV_obj {bv0:.3f} -> {bv1:.3f}, R {r0:.5e} -> {r1:.5e}, E_ox {e0:.4f} -> {e1:.4f},"
          f" V_pt {v0['V_pt_es']:.1f} -> {v1['V_pt_es']:.1f}, min c {c1.min():.3g}, {sl.message}",
          flush=True)  # fmt: skip
    j8 = {"initial": S17B.j8_check(cons, TH0), "final": S17B.j8_check(cons, th1)}
    res["final"] = {
        "theta": th1.tolist(), "bv_obj": bv1, "R_channel_h2": r1, "E_ox": e1,
        "R_over_R0_minus_1": r1 / r0 - 1.0, "E_ox_over_E0_minus_1": e1 / e0 - 1.0,
        "bound_violation": bounds.violation(th1), "constraint_values": _lst(c1),
        "V_pt_es_h": v1, "cut_limiting_x_um": ev_off.cut_x_um[o1["limiting_cut"]],
    }  # fmt: skip
    res["slsqp_main"] = {"status": int(sl.status), "message": sl.message, "nit": int(sl.nit),
                         "iterates": its, "evaluations": pr.history,
                         "n_trial_evals_crossing": int(sum(g["sign_changes"] > 0
                                                           for g in cr.guards))}  # fmt: skip
    diags["main"] = job_diag([cr], [ev_off])
    dump()

    for k in range(5):
        submit(("j2_fd", k), S17.job_fd_ron, (th1, k, S17.LEVEL_ON, "final"))
    submit(("gc1", None), S17.job_grad_off, (th1, e0, "final"))
    submit(("ron_pt", "final_L1"), S17.job_ron_point, (th1, S17.LEVEL_ON, "final", True))
    submit(("ron_pt", "final_L2"), S17.job_ron_point, (th1, S17.LEVEL_RECHECK, "final", False))
    for rt in S17.PARETO_RATIOS:
        submit(("pareto", rt), job_pareto21, (th1, rt, bv0, r0, e0, args.maxiter))
    for lev in (0, 1):
        submit(("dense", f"final_L{lev}"), S17.job_dense17, (lev, th1, "final", False))
        submit(("fl", f"final_L{lev}"), S17.job_fl17,
               (lev, th1, "final", S17.N_S_FL, S17.DS_FL_UM))  # fmt: skip
        submit(("vpt", f"final_L{lev}"), job_vpt21, (th1, lev, "final"))
    submit(("vt", "final"), job_vt21, (th1, "final"))
    submit(("ox", "final"), job_ox21, (th1, "final"))
    for lv in (1, 2):
        submit(("acc", f"final_L{lv}"), job_acc21, (th1, lv, "final"))
    kinds = ("j1_fd", "j2_fd", "pareto", "dense", "ron_pt", "fl", "vpt", "vt", "ox", "acc")
    jobs: dict = {k: {} for k in kinds}
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
                tp = np.asarray(r["theta"])
                for a in ((S17.job_dense17, "dense", f"ratio_{key:g}_L0", (0, tp, f"ratio_{key:g}", False)),
                          (S17.job_ron_point, "ron_pt", f"ratio_{key:g}_L2", (tp, S17.LEVEL_RECHECK, f"ratio_{key:g}", False)),
                          (job_vpt21, "vpt", f"ratio_{key:g}_L0", (tp, 0, f"ratio_{key:g}"))):  # fmt: skip
                    nf = pool.submit(a[0], a[3])
                    fut[nf] = (a[1], a[2])
                    pending.add(nf)
            break
    pool.shutdown()
    T["phase2_parallel"] = time.perf_counter() - t0
    dense, vpt, fl = jobs["dense"], jobs["vpt"], jobs["fl"]
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
    gc1, gf1 = jobs["gc1"]["grad_check"], jobs["gc1"]["grad_check_fields"]
    g1 = np.array(jobs["ron_pt"]["final_L1"]["grad"])
    fd1 = np.array([jobs["j2_fd"][k]["richardson"] for k in range(5)])
    a_r = float(np.max(np.abs(g1 - fd1)) / r0)
    n_r, n_bv, n_e = (float(np.linalg.norm(g1) / r0), gc1["grad_bv_norm"] / bv0,
                      gf1["grad_E_over_E0_norm"])  # fmt: skip
    a_bv = max(r["bv"]["abs_diff"] for r in gc1["directions"]) / bv0
    a_e = max(r["E_ox"]["abs_diff"] for r in gf1["directions"]) / e0
    tol = {k: S17.J2_TOL * max(1.0, n) for k, n in (("bv", n_bv), ("r", n_r), ("e", n_e))}
    crit["S1"] = bool(j1_err <= S17.J1_TOL and d2_bv <= S17.D2_TOL_REL and d2_e <= S17.D2_TOL_REL
                      and a_bv <= tol["bv"] and a_r <= tol["r"] and a_e <= tol["e"])  # fmt: skip
    vals["S1"] = {"J1_rel_err_R_on_ch": j1_err, "J1_rel_BV_obj": d2_bv, "J1_rel_E_ox": d2_e,
                  "J2_abs_BV_over_BV0": a_bv, "J2_abs_R_over_R0": a_r, "J2_abs_E_over_E0": a_e,
                  "J2_tol": tol, "J8_constraint_gradient_reported": j8}  # fmt: skip
    # S2 = J3 + J4b
    viol = bounds.violation(th1)
    r1_h2 = jobs["ron_pt"]["final_L1"]["R"]
    j3 = bool(r1 <= r0 * (1 + 1e-8) and r1_h2 <= r0 * (1 + 1e-8) and e1 <= e0 * (1 + 1e-8)
              and viol <= 1e-9)  # fmt: skip
    pr_runs = jobs["pareto"]
    cross_it = sum(i["n_crossing"] > 0 for i in its) + sum(
        i["n_crossing"] > 0 for q in pr_runs.values() for i in q["iterates"])  # fmt: skip
    cross_fin = int(S17B.sign_status(cons, th1)[1] > 0) + sum(
        q["final_sign"]["n_crossing"] > 0 for q in pr_runs.values())  # fmt: skip
    g_all = [d["guards"] for d in diags.values() if "guards" in d and d["guards"]["n"]]
    fin_g = [jobs["ron_pt"][k]["guard"] for k in jobs["ron_pt"]]
    fin_clamp = sum(g["gated_near_ns_clamp"] + g["gated_near_eeff_clamp"] for g in fin_g)
    j4b = bool(cross_it == 0 and cross_fin == 0 and fin_clamp == 0 and g_all
               and min(g["min_gated_ns_cm2"] for g in g_all) > 2e11
               and min(g["min_gated_eeff_raw_MV_cm"] for g in g_all) > 2e-2)  # fmt: skip
    crit["S2"] = bool(j3 and j4b)
    vals["S2"] = {"J3": j3, "J4b": j4b, "R_over_R0_minus_1": r1 / r0 - 1.0,
                  "E_over_E0_minus_1": e1 / e0 - 1.0, "bound_violation": viol,
                  "n_iterates_crossing": cross_it, "n_final_designs_crossing": cross_fin,
                  "clamp_violations_at_designs": fin_clamp,
                  "n_trial_evaluations_crossing_reported": int(
                      res["slsqp_main"]["n_trial_evals_crossing"]
                      + sum(q["n_trial_evals_crossing"] for q in pr_runs.values())),
                  "min_constraint_value_final": float(c1.min())}  # fmt: skip
    # S3
    ini = {lev: dense[f"initial_L{lev}"] for lev in (0, 1)}
    fin = {lev: dense[f"final_L{lev}"] for lev in (0, 1)}
    imp = {lev: fin[lev]["bv_min"] - ini[lev]["bv_min"] for lev in (0, 1)}
    crit["S3"] = bool(imp[0] > 0 and imp[1] > 0 and imp[1] >= S17.J5_RATIO * imp[0]
                      and len(ini[0]["x_um"]) >= 40)  # fmt: skip
    vals["S3"] = {"improvement_h": imp[0], "improvement_h2": imp[1],
                  "ratio": imp[1] / imp[0] if imp[0] else None,
                  "min_bv_initial_h_h2": [ini[0]["bv_min"], ini[1]["bv_min"]],
                  "min_bv_final_h_h2": [fin[0]["bv_min"], fin[1]["bv_min"]]}  # fmt: skip
    # S4: iterates (BV_obj), final + Pareto (max of BV_obj and dense true-min BV), h/2 final
    it_all = list(its) + [i for q in pr_runs.values() for i in q["iterates"]]
    it_bad = [i for i in it_all if not i["guard_ok"]]
    bv_ref_final = max(bv1, fin[0]["bv_min"])
    fin_ok = guard_ok(v1, bv_ref_final)
    par_rows = []
    for rt, q in pr_runs.items():
        bvd = dense.get(f"ratio_{rt:g}_L0", {}).get("bv_min", float("nan"))
        ref = max(q["bv_obj"], bvd)
        par_rows.append({"ratio": rt, "V_pt_es": q["final_vpt"]["V_pt_es"], "bv_obj": q["bv_obj"],
                         "bv_dense": bvd, "V_pt_over_BV": q["final_vpt"]["V_pt_es"] / ref,
                         "ok": guard_ok(q["final_vpt"], ref)})  # fmt: skip
    v_h2 = vpt["final_L1"]
    h2_ok = guard_ok(v_h2, max(bv1, fin[1]["bv_min"]))
    crit["S4"] = bool(not it_bad and fin_ok and all(p["ok"] for p in par_rows) and h2_ok)
    vals["S4"] = {
        "n_iterates_checked": len(it_all), "n_iterates_violating": len(it_bad),
        "min_iterate_ratio": min((i["ratio"] for i in it_all if i["ratio"] is not None),
                                 default=None),
        "first_violating_iterates": [{k: i[k] for k in ("theta", "bv_obj", "V_pt_es", "ratio")}
                                     for i in it_bad[:3]],
        "final_h": {"V_pt_es": v1["V_pt_es"], "BV_ref": bv_ref_final,
                    "ratio": v1["V_pt_es"] / bv_ref_final, "ok": fin_ok},
        "final_h2": {"V_pt_es": v_h2["V_pt_es"], "ratio": v_h2["V_pt_es"] / max(bv1, fin[1]["bv_min"]),
                     "ok": h2_ok},
        "pareto": par_rows, "ceiling_hits": int(sum(i["ceiling"] for i in it_all)),
        "iterate_vpt_errors": int(sum(i["pt_error"] is not None for i in it_all)),
    }  # fmt: skip
    # S5
    s5 = {}
    for lv in (0, 1):
        f_ = fl[f"final_L{lv}"]
        d_ = fin[lv]["bv_min"]
        s5[f"L{lv}"] = {"fieldline": f_.get("bv"), "dense_cut_min": d_,
                        "rel_diff": abs(f_.get("bv", np.nan) / d_ - 1.0)}  # fmt: skip
    crit["S5"] = bool(all(v["rel_diff"] <= S5_TOL for v in s5.values()))
    vals["S5"] = s5
    # S6
    nf_tot = sum(d["nonfinite_calls"] for d in diags.values())
    nf_failed = sum(d["NF"]["failed_newton"] for d in diags.values())
    nf_adj = sum(d["NF"]["adjoint"] for d in diags.values())
    nf_conv = sum(d["NF"]["converged_newton"] for d in diags.values())

    def agg(a, b):
        return sum(d[a].get(b, 0) for d in diags.values() if a in d)

    stderr = log_path.read_text(errors="replace") if log_path.exists() else ""
    entry = [ln for ln in stderr.splitlines() if "On entry to" in ln]
    errs = [f"{kd}|{k}" for kd in kinds for k, v in jobs[kd].items() if "error" in v]
    errs += [k for k in ("gc0", "gc1") if jobs[k] is None or "error" in jobs[k]]
    acc_bad = [k for k, v in jobs["acc"].items() if "error" not in v and not v["converged"]]
    crit["S6"] = bool(
        sl.status == 0 and sl.nit <= args.maxiter and nf_tot == nf_failed and nf_adj == 0
        and nf_conv == 0 and agg("on_state", "failed_picard") == 0
        and agg("on_state", "failed_evals") == 0 and agg("off_state", "failed_evals") == 0
        and not entry and not errs and not acc_bad
        and all(d["resolved_all"] for d in dense.values())
        and all(v.get("converged", True) for v in fl.values())
    )  # fmt: skip
    vals["S6"] = {"slsqp_status": int(sl.status), "slsqp_message": sl.message, "nit": int(sl.nit),
                  "pareto_status": {str(k): v["status"] for k, v in pr_runs.items()},
                  "nonfinite_calls_total": nf_tot, "nonfinite_in_failed_newton": nf_failed,
                  "nonfinite_in_adjoint": nf_adj, "nonfinite_in_converged_newton": nf_conv,
                  "failed_picard": agg("on_state", "failed_picard"),
                  "failed_evals_on": agg("on_state", "failed_evals"),
                  "failed_evals_off": agg("off_state", "failed_evals"),
                  "on_entry_to_lines": len(entry), "job_errors": errs, "acc_unconverged": acc_bad}  # fmt: skip

    # ================================================================ R1 / R2
    def jr(kind, key, field):
        r = jobs[kind].get(key, {})
        return r.get(field)

    par = [{"ratio": None, "note": "initial", "R": r0, "bv_true_min": ini[0]["bv_min"]},
           {"ratio": 1.0, "theta": th1.tolist(), "R": r1, "bv_obj": bv1,
            "bv_true_min": fin[0]["bv_min"]}]  # fmt: skip
    for rt in S17.PARETO_RATIOS:
        q = dict(pr_runs[rt])
        q["bv_true_min"] = dense[f"ratio_{rt:g}_L0"]["bv_min"]
        q["V_pt_es"] = q["final_vpt"]["V_pt_es"]
        par.append(q)
    for q in par:
        q["R_over_R0"] = q["R"] / r0
    res["R1"] = {
        "bv_dense_h": [ini[0]["bv_min"], fin[0]["bv_min"]],
        "bv_dense_h2": [ini[1]["bv_min"], fin[1]["bv_min"]],
        "bv_fieldline_h": [jr("fl", "initial_L0", "bv"), jr("fl", "final_L0", "bv")],
        "bv_fieldline_h2": [jr("fl", "initial_L1", "bv"), jr("fl", "final_L1", "bv")],
        "bv_obj": [bv0, bv1],
        "R_on_ch_h2": [r0, r1_h2],
        "R_on_ch_h4": [jr("ron_pt", "initial_L2", "R"), jr("ron_pt", "final_L2", "R")],
        "accumulation_bracket_Lp_U": {k: {"R_U": v.get("R_U"), "R_Lp": v.get("R_Lp")}
                                      for k, v in jobs["acc"].items()},
        "V_T_es_V": [jr("vt", "initial", "V_T_es_V"), jr("vt", "final", "V_T_es_V")],
        "E_ox_400V": {"initial": jobs["ox"].get("initial"), "final": jobs["ox"].get("final")},
        "V_pt_es": {"initial_h": jr("vpt", "initial_L0", "V_pt_es"),
                    "final_h": jr("vpt", "final_L0", "V_pt_es"),
                    "initial_h2": jr("vpt", "initial_L1", "V_pt_es"),
                    "final_h2": jr("vpt", "final_L1", "V_pt_es")},
        "V_pt_over_BV_dense": {"initial_h": v0["V_pt_es"] / ini[0]["bv_min"],
                               "final_h": v1["V_pt_es"] / fin[0]["bv_min"]},
        "theta_initial": TH0.tolist(), "theta_final": th1.tolist(), "pareto": par,
        "active_constraints": [
            {**cons.labels()[i], "value": float(c1[i])} for i in np.argsort(c1) if c1[i] < 1e-6],
        "active_R": bool(abs(r1 / r0 - 1) < 1e-6), "active_E_ox": bool(abs(e1 / e0 - 1) < 1e-6),
    }  # fmt: skip
    r17b = None
    if args.stage17b_results.exists():
        d17 = json.loads(args.stage17b_results.read_text())
        r17b = {"source": str(args.stage17b_results), "dev": d17.get("dev"),
                "theta": d17.get("final", {}).get("theta"),
                "bv_obj": [d17.get("initial", {}).get("bv_obj"), d17.get("final", {}).get("bv_obj")],
                "bv_dense_h": d17.get("R4", {}).get("bv_true_min_h"),
                "note": "DIFFERENT GEOMETRY: 0.8 um body (punches through at ~90-130 V)"}  # fmt: skip
    res["R2"] = {
        "stage17b_0p8um": r17b,
        "stage21_1p4um": {
            "theta": th1.tolist(),
            "bv_obj": [bv0, bv1],
            "bv_dense_h": [ini[0]["bv_min"], fin[0]["bv_min"]],
        },
    }
    res["diagnostics"], res["seconds"] = diags, T
    res["runtime_seconds"] = time.perf_counter() - t_start
    figdir = (REPO / "figures") if out == DEFAULT_OUT.resolve() else out / "figures"
    label = f"stage 21 | commit {str(man.get('git_commit'))[:8]} | {man['timestamp_utc'][:19]}Z"
    try:
        plot_all(figdir, res, label)
    except Exception as exc:
        res["figures_error"] = f"{type(exc).__name__}: {exc}"
    summary = {"criteria": crit, "gate_passed": bool(all(crit.values())), "dev": bool(args.dev),
               "values": {"bv_obj_initial": bv0, "bv_obj_final": bv1, "R_ch_initial": r0,
                          "R_ch_final": r1, **vals, "R1": res["R1"], "R2": res["R2"]},
               "seconds": T, "runtime_seconds": res["runtime_seconds"]}  # fmt: skip
    dump()
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    print(json.dumps(crit, indent=2))
    print(f"runtime {res['runtime_seconds']:.1f} s")


if __name__ == "__main__":
    main()
