"""Stage-24 benchmark (docs/benchmarks/stage24_guard_ablation.md): X1-X4 (gating), R1-R2, P24a-c.

Experiment A: the stage-21 problem on the stage-20 d_body = 1.2 um geometry, solved twice.
Arm U has no punch-through information.  Arm G adds the SLSQP constraint
B_tau(theta; V* = 650 V) >= 10 V_T (``widegap.adjoint_barrier``: nested log-sum-exp soft-minima,
tau = 0.5 V_T, over y and over the stage-20 barrier columns frozen at theta = 0; 12c off-state at
fixed V_D on h; off-state adjoint).  If B_tau(0) < 10 V_T the pre-registered restoration runs
first: SLSQP minimizing -B_tau/10 under the common constraints (13d settings, maxiter 60: the spec
says "short" but gives no iteration cap, so it runs to SLSQP convergence); its final iterate is the
Arm-G start (recorded).  The predictions P24a-c are evaluated on h (and reported on h/2); they are
reported whatever the gate outcome.

Geometry plumbing is stage 21's (``benchmarks/stage21/run.py``), re-installed for d_body = 1.2 um
(13d machinery on the deeper device; W1/W2 shifted down by 0.4 um, unused here).

Writes <out>/{results,summary,manifest}.json, stderr.log (default results/stage24) and figures
stage24_*.png.  ``--dev`` marks a DEV run (refused under results/).
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

from widegap.adjoint_barrier import LEVEL_VT, TAU_VT, V_STAR_V, BarrierAdjoint
from widegap.adjoint_sic import richardson_fd
from widegap.bodyredesign import body_design, column_index_matrix, source_columns
from widegap.manifest import build_manifest

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[2]
DEFAULT_OUT = REPO / "results" / "stage24"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


S21 = _load("stage21_run", REPO / "benchmarks/stage21/run.py")  # fixed name: fork pickling
S17B, S17, S20 = S21.S17B, S21.S17, S21.S20
M13D, M = S17.M13D, S17.M
tracked, job_diag = S17.tracked, S17.job_diag

D_BODY_UM = 1.2
Y_SHIFT_UM = D_BODY_UM - 0.8
S21.D_BODY_UM = D_BODY_UM
S21.Y_SHIFT_UM = Y_SHIFT_UM
S21.DESIGN = body_design(S21.DESIGN0, D_BODY_UM)
S21.WINDOWS_UM = {
    "W1": ((0.5, 1.1), (-1.1 - Y_SHIFT_UM, -0.5 - Y_SHIFT_UM)),
    "W2": ((1.3, 2.0), (-1.2 - Y_SHIFT_UM, -0.6 - Y_SHIFT_UM)),
    "W3": ((0.0, 0.4), (-0.2, 0.0)),
}
S21.install_geometry()
DESIGN, TH0, Bounds = S21.DESIGN, S17.TH0, S17.Bounds

X1_TOL, X4_BV_TOL, X4_VPT_TOL = 1e-4, 0.03, 0.05
P24C_GAIN_MARGIN_V = 5.0
P24B_FACTOR = 1.05
FEAS = 1e-8


def _lst(a):
    return np.asarray(a, dtype=float).tolist()


def make_barrier():
    ev = M13D.make_ev13d(0, with_cuts=False)
    ev.set_theta(TH0)
    idx = column_index_matrix(ev, source_columns(ev))
    return BarrierAdjoint(ev, idx, V_STAR_V, TAU_VT), ev


# ============================================================ SLSQP (arms and restoration)
class BarrierCons:
    """c(theta) = B_tau / 10 V_T - 1 >= 0 with a one-point cache."""

    def __init__(self, ba):
        self.ba, self.key, self.val = ba, None, None

    def _get(self, th):
        key = np.asarray(th, float).tobytes()
        if key != self.key:
            self.val, self.key = self.ba.evaluate(th), key
        return self.val

    def fun(self, th):
        return self._get(th)["B"] / LEVEL_VT - 1.0

    def jac(self, th):
        return self._get(th)["grad"] / LEVEL_VT


def run_arm_slsqp(comp, cons, ba, th_start, bv0, r0, e0, bounds, maxiter, guarded, restore=False,
                  ftol=1e-10):  # fmt: skip
    pr = M13D.Problem(comp, bv0, r0, e0, 1.0, bounds)
    bc = BarrierCons(ba)
    cl = [
        {"type": "ineq", "fun": pr.c_ron, "jac": pr.c_ron_jac},
        {"type": "ineq", "fun": pr.c_ox, "jac": pr.c_ox_jac},
        {"type": "ineq", "fun": lambda t: bounds.A @ t + bounds.b, "jac": lambda t: bounds.A},
        {"type": "ineq", "fun": cons.values, "jac": cons.jac},
    ]
    if guarded and not restore:
        cl.append({"type": "ineq", "fun": bc.fun, "jac": bc.jac})
    fun, jac = (pr.fun, pr.jac) if not restore else (
        lambda t: -bc.fun(t) - 1.0, lambda t: -bc.jac(t))  # fmt: skip
    iterates: list = []

    def cb(xk):
        xk = np.asarray(xk, float)
        m, n = S17B.sign_status(cons, xk)
        iterates.append({"theta": xk.tolist(), "min_margin_plus_m": m, "n_crossing": n,
                         "B_tau": float(ba.evaluate(xk, grad=False)["B"])})  # fmt: skip

    try:
        res = minimize(fun, np.asarray(th_start, float), jac=jac, method="SLSQP",
                       bounds=bounds.box, constraints=cl,
                       options={"maxiter": maxiter, "ftol": ftol}, callback=cb)  # fmt: skip
        return np.asarray(res.x, float), res, pr, iterates
    except M.EvalAbort as exc:
        res = SimpleNamespace(status=-1, message=f"ABORTED (forward solve failed): {exc}",
                              nit=len(iterates), nfev=pr.n_eval, njev=pr.n_eval,
                              x=np.asarray(th_start, float))  # fmt: skip
        return np.asarray(th_start, float), res, pr, iterates


@tracked
def job_arm(args):
    arm, bv0, r0, e0, maxiter = args
    t0 = time.perf_counter()
    ev_off = S17.make_ev_off()
    cr = S17.make_cr(S17.LEVEL_ON)
    cons = S17B.make_cons()
    ba, ev_bar = make_barrier()
    comp = S17.Eval17(ev_off, cr)
    bounds = Bounds(DESIGN, ev_off.sc)
    out = {"arm": arm}
    start = TH0
    if arm == "G":
        b0 = ba.evaluate(TH0, grad=False)
        out["B_tau_at_theta0"] = b0["B"]
        out["hard_barrier_at_theta0"] = b0["hard"]
        if b0["B"] < LEVEL_VT * (1.0 - FEAS):
            th, rs, _, its = run_arm_slsqp(comp, cons, ba, TH0, bv0, r0, e0, bounds, 60, True,
                                           restore=True)  # fmt: skip
            bb = ba.evaluate(th, grad=False)
            out["restoration"] = {"status": int(rs.status), "message": rs.message,
                                  "nit": int(rs.nit), "theta": th.tolist(), "B_tau": bb["B"],
                                  "hard": bb["hard"], "iterates": its,
                                  "reached_feasible_B": bool(bb["B"] >= LEVEL_VT * (1 - FEAS))}  # fmt: skip
            start = th
        else:
            out["restoration"] = None
    th, res, pr, its = run_arm_slsqp(comp, cons, ba, start, bv0, r0, e0, bounds, maxiter,
                                     arm == "G")  # fmt: skip
    b = ev_off.objective(th, grad=False)
    f = ev_off.fields(th, grad=False)
    bb = ba.evaluate(th, grad=False)
    mn, nc = S17B.sign_status(cons, th)
    out.update({
        "start_theta": np.asarray(start).tolist(), "theta": th.tolist(), "bv_obj": b["BV"],
        "R": cr.ron(th, grad=False)["R"], "E_ox": f["E_ox"]["S"], "B_tau": bb["B"],
        "hard_barrier_650": bb["hard"], "status": int(res.status), "message": res.message,
        "nit": int(res.nit), "bound_violation": bounds.violation(th),
        "constraint_values": _lst(cons.values(th)),
        "final_sign": {"min_margin_plus_m": mn, "n_crossing": nc}, "iterates": its,
        "n_trial_evals_crossing": int(sum(g["sign_changes"] > 0 for g in cr.guards)),
        "failure": pr.failure, "n_evals": pr.n_eval,
        "diag": job_diag([cr], [ev_off, ev_bar]), "seconds": time.perf_counter() - t0,
    })  # fmt: skip
    return out


@tracked
def job_barrier_fd(args):
    """B_tau adjoint gradient and Richardson FD along one unit direction (kind 'fd') at theta."""
    th, k, tag = args
    t0 = time.perf_counter()
    ba, ev = make_barrier()
    th = np.asarray(th, float)
    if k is None:
        o = ba.evaluate(th)
        out = {"tag": tag, "B": o["B"], "hard": o["hard"], "grad": _lst(o["grad"])}
    else:
        u = np.zeros(5)
        u[k] = 1.0
        r, f1, f2 = richardson_fd(lambda t: ba.evaluate(t, grad=False)["B"], th, u, 1e-3)
        out = {"tag": tag, "k": k, "richardson": float(r), "fd_eps": float(f1),
               "fd_eps_half": float(f2)}  # fmt: skip
    out["diag"] = job_diag([], [ev])
    out["seconds"] = time.perf_counter() - t0
    return out


@tracked
def job_bar650(args):
    th, tag = args
    t0 = time.perf_counter()
    ba, ev = make_barrier()
    o = ba.evaluate(np.asarray(th, float), grad=False)
    return {"tag": tag, "B_tau": o["B"], "hard_barrier_650": o["hard"],
            "diag": job_diag([], [ev]), "seconds": time.perf_counter() - t0}  # fmt: skip


# ============================================================ figures
def plot_all(figdir, res, label):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figdir = Path(figdir)
    figdir.mkdir(parents=True, exist_ok=True)
    foot = lambda fig: fig.text(0.5, 0.005, label, ha="center", fontsize=7, color="0.4")  # noqa: E731
    fig, ax = plt.subplots(figsize=(7.2, 4.4), constrained_layout=True)
    pts = res["R1"]["blocking"]
    labs, bv, vp = [], [], []
    for k in ("initial", "U", "G"):
        if k in pts:
            labs.append(k)
            bv.append(pts[k]["bv_dense_h"])
            vp.append(pts[k]["V_pt_es_h"])
    x = np.arange(len(labs))
    ax.bar(x - 0.2, bv, 0.4, label="dense avalanche BV (h)")
    ax.bar(x + 0.2, vp, 0.4, label="V_pt,es (h)")
    ax.set_xticks(x, labs)
    ax.set_ylabel("[V]")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")
    foot(fig)
    fig.savefig(figdir / "stage24_blocking.png", dpi=130)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(7.2, 4.0), constrained_layout=True)
    for arm, c in (("U", "#D55E00"), ("G", "#0072B2")):
        its = res["arms"].get(arm, {}).get("iterates", [])
        ax.plot([i["B_tau"] for i in its], "o-", color=c, label=f"arm {arm}")
    ax.axhline(LEVEL_VT, color="k", ls="--", label="10 V_T")
    ax.set_xlabel("SLSQP iterate")
    ax.set_ylabel("B_tau(650 V) [V_T]")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    foot(fig)
    fig.savefig(figdir / "stage24_barrier_iterates.png", dpi=130)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(7.2, 4.0), constrained_layout=True)
    th = {k: res["arms"][k]["theta"] for k in ("U", "G") if k in res["arms"]}
    x = np.arange(5)
    for j, (k, v) in enumerate(th.items()):
        ax.bar(x + 0.4 * (j - 0.5), v, 0.4, label=f"arm {k}")
    ax.set_xticks(x, [f"theta_{i}" for i in range(4)] + ["theta_J"])
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")
    foot(fig)
    fig.savefig(figdir / "stage24_theta.png", dpi=130)
    plt.close(fig)


# ============================================================ main
def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
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
    ba, ev_bar = make_barrier()
    bounds = Bounds(DESIGN, sc)
    res: dict = {"dev": bool(args.dev)}
    diags: dict = {}

    def dump():
        (out / "results.json").write_text(json.dumps(res, indent=2, default=float) + "\n")

    o0 = ev_off.objective(TH0, grad=False)
    f0 = ev_off.fields(TH0, grad=False)
    r0 = cr.ron(TH0, grad=False)["R"]
    bv0, e0 = o0["BV"], f0["E_ox"]["S"]
    b0 = ba.evaluate(TH0, grad=False)
    res["config"] = {
        "d_body_um": D_BODY_UM, "V_star_V": V_STAR_V, "tau_VT": TAU_VT, "level_VT": LEVEL_VT,
        "n_barrier_columns": int(ba.measure.idx.shape[0]), "n_window_nodes": int(ba.measure.idx.shape[1]),
        "n_classification_constraints": cons.n, "max_iter": args.maxiter, "workers": args.workers,
        "mesh_nodes_off_L0": int(ev_off.mesh.n_nodes), "registry": S17.registry_check(),
        "restoration_maxiter": 60,
    }  # fmt: skip
    res["initial"] = {"bv_obj": bv0, "R_channel_h2": r0, "E_ox": e0, "B_tau_650": b0["B"],
                      "hard_barrier_650": b0["hard"]}  # fmt: skip
    print(f"initial: BV_obj {bv0:.3f}, R {r0:.6e}, E_ox {e0:.4f}, B_tau(650) {b0['B']:.3f} VT "
          f"(hard {b0['hard']:.3f})", flush=True)  # fmt: skip

    pool = ProcessPoolExecutor(args.workers, mp_context=mp.get_context("fork"))
    fut: dict = {}

    def submit(key, fn, a):
        fut[pool.submit(fn, a)] = key

    for arm in ("U", "G"):
        submit(("arm", arm), job_arm, (arm, bv0, r0, e0, args.maxiter))
    submit(("bfd0", "adj"), job_barrier_fd, (TH0, None, "initial"))
    for k in range(5):
        submit(("bfd0", k), job_barrier_fd, (TH0, k, "initial"))
        submit(("j1_fd", k), S17.job_fd_ron, (TH0, k, S17.LEVEL_ON, "initial"))
    submit(("gc0", None), S17.job_grad_off, (TH0, e0, "initial"))
    submit(("ron_pt", "initial_L1"), S17.job_ron_point, (TH0, S17.LEVEL_ON, "initial", True))
    for lev in (0, 1):
        submit(("dense", f"initial_L{lev}"), S17.job_dense17, (lev, TH0, "initial", False))
        submit(("vpt", f"initial_L{lev}"), S21.job_vpt21, (TH0, lev, "initial"))
        submit(("fl", f"initial_L{lev}"), S17.job_fl17,
               (lev, TH0, "initial", S17.N_S_FL, S17.DS_FL_UM))  # fmt: skip
    submit(("bar650", "initial"), job_bar650, (TH0, "initial"))
    kinds = ("bfd0", "j1_fd", "ron_pt", "dense", "vpt", "fl", "bar650", "arm", "ox")
    jobs: dict = {k: {} for k in kinds}
    for arm in ("U", "G"):
        for k in ("j2_fd", "bfd", "gc"):
            jobs[f"{k}_{arm}"] = {}
    jobs["gc0"] = None
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
            if kind == "gc0":
                jobs["gc0"] = r
            else:
                jobs.setdefault(kind, {})[key] = r
            if "diag" in r:
                diags[f"{kind}|{key}"] = r["diag"]
            print(f"done {kind} {key} {r.get('seconds', float('nan')):.0f}s"
                  f"{' ERROR ' + r['error'] if 'error' in r else ''}", flush=True)  # fmt: skip
            if kind == "arm" and "error" not in r:
                arm, th = key, np.asarray(r["theta"])
                new = []
                for k in range(5):
                    new.append(
                        (f"j2_fd_{arm}", k, S17.job_fd_ron, (th, k, S17.LEVEL_ON, f"final{arm}"))
                    )
                    new.append((f"bfd_{arm}", k, job_barrier_fd, (th, k, f"final{arm}")))
                new.append((f"bfd_{arm}", "adj", job_barrier_fd, (th, None, f"final{arm}")))
                new.append((f"gc_{arm}", None, S17.job_grad_off, (th, e0, f"final{arm}")))
                new.append(("ron_pt", f"final{arm}_L1", S17.job_ron_point,
                            (th, S17.LEVEL_ON, f"final{arm}", True)))  # fmt: skip
                new.append(("bar650", f"final{arm}", job_bar650, (th, f"final{arm}")))
                new.append(("ox", f"final{arm}", S21.job_ox21, (th, f"final{arm}")))
                for lev in (0, 1):
                    new.append(("dense", f"final{arm}_L{lev}", S17.job_dense17,
                                (lev, th, f"final{arm}", False)))  # fmt: skip
                    new.append(("vpt", f"final{arm}_L{lev}", S21.job_vpt21,
                                (th, lev, f"final{arm}")))  # fmt: skip
                    new.append(("fl", f"final{arm}_L{lev}", S17.job_fl17,
                                (lev, th, f"final{arm}", S17.N_S_FL, S17.DS_FL_UM)))  # fmt: skip
                for kd, key2, fn, a in new:
                    nf = pool.submit(fn, a)
                    fut[nf] = (kd, key2)
                    pending.add(nf)
            break
    pool.shutdown()
    arms = {a: jobs["arm"][a] for a in ("U", "G")}
    res["arms"], res["jobs"] = arms, jobs
    res["seconds"] = {"parallel": time.perf_counter() - t0}
    dump()

    # ================================================================ criteria
    dense, vpt = jobs["dense"], jobs["vpt"]
    crit, vals = {}, {}

    def fdnorm(adj, fd):
        return float(np.linalg.norm(np.array(adj) - np.array(fd)) / np.linalg.norm(fd))

    x1 = {}
    for tag, bj, adjk in (("theta0", jobs["bfd0"], "adj"), ("G_final", jobs["bfd_G"], "adj"),
                          ("U_final", jobs["bfd_U"], "adj")):  # fmt: skip
        fdv = [bj[k]["richardson"] for k in range(5)]
        x1[tag] = {"rel_err_norm": fdnorm(bj[adjk]["grad"], fdv), "B": bj[adjk]["B"],
                   "grad_norm": float(np.linalg.norm(bj[adjk]["grad"]))}  # fmt: skip
    # stage-21 S1-type checks for the common gradients: theta = 0 (J1) and each final (J2)
    g0 = jobs["ron_pt"]["initial_L1"]
    fd0 = np.array([jobs["j1_fd"][k]["richardson"] for k in range(5)])
    j1_err = float(np.linalg.norm(np.array(g0["grad"]) - fd0) / np.linalg.norm(fd0))
    gc0, gf0 = jobs["gc0"]["grad_check"], jobs["gc0"]["grad_check_fields"]
    d2_bv = max(r["bv"]["rel_err"] for r in gc0["directions"])
    d2_e = max(r["E_ox"]["rel_err"] for r in gf0["directions"])
    s1 = bool(j1_err <= S17.J1_TOL and d2_bv <= S17.D2_TOL_REL and d2_e <= S17.D2_TOL_REL)
    s1_arm = {}
    for arm in ("U", "G"):
        gcf, gff = (
            jobs[f"gc_{arm}"][None]["grad_check"],
            jobs[f"gc_{arm}"][None]["grad_check_fields"],
        )
        g1 = np.array(jobs["ron_pt"][f"final{arm}_L1"]["grad"])
        fd1 = np.array([jobs[f"j2_fd_{arm}"][k]["richardson"] for k in range(5)])
        a_r = float(np.max(np.abs(g1 - fd1)) / r0)
        n_r, n_bv, n_e = (float(np.linalg.norm(g1) / r0), gcf["grad_bv_norm"] / bv0,
                          gff["grad_E_over_E0_norm"])  # fmt: skip
        a_bv = max(r["bv"]["abs_diff"] for r in gcf["directions"]) / bv0
        a_e = max(r["E_ox"]["abs_diff"] for r in gff["directions"]) / e0
        tol = {k: S17.J2_TOL * max(1.0, n) for k, n in (("bv", n_bv), ("r", n_r), ("e", n_e))}
        s1_arm[arm] = {"ok": bool(a_bv <= tol["bv"] and a_r <= tol["r"] and a_e <= tol["e"]),
                       "abs_BV": a_bv, "abs_R": a_r, "abs_E": a_e, "tol": tol}  # fmt: skip
        s1 = s1 and s1_arm[arm]["ok"]
    crit["X1"] = bool(x1["theta0"]["rel_err_norm"] <= X1_TOL
                      and x1["G_final"]["rel_err_norm"] <= X1_TOL and s1)  # fmt: skip
    vals["X1"] = {"B_tau": x1, "tol": X1_TOL, "S1_theta0": {"R": j1_err, "BV": d2_bv, "E_ox": d2_e},
                  "S1_finals": s1_arm}  # fmt: skip
    # X2
    x2 = {}
    for arm, a in arms.items():
        viol = a["bound_violation"]
        r_h2 = jobs["ron_pt"][f"final{arm}_L1"]["R"]
        ok = bool(a["R"] <= r0 * (1 + FEAS) and r_h2 <= r0 * (1 + FEAS) and a["E_ox"] <= e0 * (1 + FEAS)
                  and viol <= 1e-9 and a["final_sign"]["n_crossing"] == 0
                  and all(i["n_crossing"] == 0 for i in a["iterates"]))  # fmt: skip
        if arm == "G":
            ok = bool(ok and a["B_tau"] >= LEVEL_VT * (1 - FEAS))
        x2[arm] = {"ok": ok, "R_over_R0_minus_1": a["R"] / r0 - 1.0,
                   "E_over_E0_minus_1": a["E_ox"] / e0 - 1.0, "bound_violation": viol,
                   "B_tau_650_VT": a["B_tau"], "min_constraint_value": min(a["constraint_values"]),
                   "n_iterates_crossing": int(sum(i["n_crossing"] > 0 for i in a["iterates"])),
                   "n_trial_evals_crossing_reported": a["n_trial_evals_crossing"]}  # fmt: skip
    clamp = sum(g["gated_near_ns_clamp"] + g["gated_near_eeff_clamp"]
                for k, v in jobs["ron_pt"].items() for g in [v["guard"]])  # fmt: skip
    crit["X2"] = bool(all(v["ok"] for v in x2.values()) and clamp == 0)
    vals["X2"] = {**x2, "clamp_margin_violations_at_designs": clamp}
    # X3
    nf_tot = sum(d["nonfinite_calls"] for d in diags.values())
    nf_failed = sum(d["NF"]["failed_newton"] for d in diags.values())
    nf_adj = sum(d["NF"]["adjoint"] for d in diags.values())
    nf_conv = sum(d["NF"]["converged_newton"] for d in diags.values())

    def agg(a_, b_):
        return sum(d[a_].get(b_, 0) for d in diags.values() if a_ in d)

    stderr = log_path.read_text(errors="replace") if log_path.exists() else ""
    entry = [ln for ln in stderr.splitlines() if "On entry to" in ln]
    errs = [f"{kd}|{k}" for kd, d in jobs.items() if isinstance(d, dict)
            for k, v in d.items() if isinstance(v, dict) and "error" in v]  # fmt: skip
    errs += [f"{kd}|{k}" for kd, d in jobs.items() if isinstance(d, dict) for k, v in d.items()
             if isinstance(v, dict) and v.get("error") is not None and "seconds" not in v]  # fmt: skip
    errs += ["gc0"] if jobs["gc0"] is None or "error" in jobs["gc0"] else []
    rest = arms["G"].get("restoration")
    st_ok = all(a["status"] == 0 for a in arms.values()) and (rest is None or rest["status"] == 0)
    crit["X3"] = bool(st_ok and nf_tot == nf_failed and nf_adj == 0 and nf_conv == 0
                      and agg("on_state", "failed_picard") == 0 and agg("on_state", "failed_evals") == 0
                      and agg("off_state", "failed_evals") == 0 and not entry and not errs
                      and all(d["resolved_all"] for d in dense.values())
                      and all(v.get("converged", True) for v in jobs["fl"].values()))  # fmt: skip
    vals["X3"] = {"status": {a: arms[a]["status"] for a in arms}, "messages": {a: arms[a]["message"]
                                                                              for a in arms},
                  "nit": {a: arms[a]["nit"] for a in arms},
                  "restoration": None if rest is None else {k: rest[k] for k in
                                                           ("status", "message", "nit", "B_tau", "hard",
                                                            "reached_feasible_B")},
                  "nonfinite_calls_total": nf_tot, "nonfinite_in_failed_newton": nf_failed,
                  "nonfinite_in_adjoint": nf_adj, "nonfinite_in_converged_newton": nf_conv,
                  "failed_picard": agg("on_state", "failed_picard"),
                  "failed_evals_on": agg("on_state", "failed_evals"),
                  "failed_evals_off": agg("off_state", "failed_evals"),
                  "on_entry_to_lines": len(entry), "job_errors": errs}  # fmt: skip
    # X4 + per-arm table
    table: dict = {}
    names = {"initial": "initial", "U": "finalU", "G": "finalG"}
    for nm, tag in names.items():
        d0, d1 = dense[f"{tag}_L0"], dense[f"{tag}_L1"]
        v0, v1 = vpt[f"{tag}_L0"], vpt[f"{tag}_L1"]
        fl0, fl1 = jobs["fl"][f"{tag}_L0"], jobs["fl"][f"{tag}_L1"]
        row = {"bv_dense_h": d0["bv_min"], "bv_dense_h2": d1["bv_min"],
               "bv_fieldline_h": fl0.get("bv"), "bv_fieldline_h2": fl1.get("bv"),
               "V_pt_es_h": v0["V_pt_es"], "V_pt_es_h2": v1["V_pt_es"],
               "V_pt_ceiling_h": v0["ceiling"], "blocking_h": min(d0["bv_min"], v0["V_pt_es"]),
               "blocking_h2": min(d1["bv_min"], v1["V_pt_es"]),
               "hard_barrier_650_VT": jobs["bar650"][tag if nm != "initial" else "initial"]["hard_barrier_650"],
               "B_tau_650_VT": jobs["bar650"][tag if nm != "initial" else "initial"]["B_tau"],
               "x_argmin_dense_h_um": d0["x_argmin_um"]}  # fmt: skip
        if nm != "initial":
            a = arms[nm]
            row.update({"theta": a["theta"], "bv_obj": a["bv_obj"], "R_on_ch_h2": a["R"],
                        "E_ox": a["E_ox"], "start_theta": a["start_theta"],
                        "active_constraints": {
                            "R": bool(abs(a["R"] / r0 - 1) < 1e-6),
                            "E_ox": bool(abs(a["E_ox"] / e0 - 1) < 1e-6),
                            "classification": [i for i, c in enumerate(a["constraint_values"]) if c < 1e-6],
                            "B_tau": bool(nm == "G" and abs(a["B_tau"] / LEVEL_VT - 1) < 1e-6),
                            "bounds": [i for i, (t, (lo, hi)) in enumerate(zip(a["theta"], bounds.box, strict=True))
                                       if min(t - lo, hi - t) < 1e-9]}})  # fmt: skip
        table[nm] = row
    x4 = {}
    for nm in ("U", "G"):
        r_ = table[nm]
        x4[nm] = {
            "dV_pt_rel": abs(r_["V_pt_es_h2"] / r_["V_pt_es_h"] - 1.0),
            "dBV_rel": abs(r_["bv_dense_h2"] / r_["bv_dense_h"] - 1.0),
        }
        x4[nm]["ok"] = bool(x4[nm]["dV_pt_rel"] <= X4_VPT_TOL and x4[nm]["dBV_rel"] <= X4_BV_TOL)
    crit["X4"] = bool(all(v["ok"] for v in x4.values()))
    vals["X4"] = x4

    def preds(suffix):
        i, u, g = table["initial"], table["U"], table["G"]
        bv, vp, bl = f"bv_dense_{suffix}", f"V_pt_es_{suffix}", f"blocking_{suffix}"
        gain_u, gain_g = u[bv] - i[bv], g[bv] - i[bv]
        return {
            "P24a": {"V_pt_U": u[vp], "BV_U": u[bv], "confirmed": bool(u[vp] < u[bv])},
            "P24b": {
                "blocking_U": u[bl],
                "blocking_initial": i[bl],
                "ratio": u[bl] / i[bl],
                "confirmed": bool(u[bl] < P24B_FACTOR * i[bl]),
            },
            "P24c": {
                "V_pt_G": g[vp],
                "BV_G": g[bv],
                "gain_U": gain_u,
                "gain_G": gain_g,
                "gain_U_minus_gain_G": gain_u - gain_g,
                "blocking_G": g[bl],
                "blocking_U": u[bl],
                "confirmed": bool(
                    g[vp] >= g[bv] and gain_u - gain_g >= P24C_GAIN_MARGIN_V and g[bl] > u[bl]
                ),
                "parts": {
                    "G_not_punchthrough_limited": bool(g[vp] >= g[bv]),
                    "gain_gap_ge_5V": bool(gain_u - gain_g >= P24C_GAIN_MARGIN_V),
                    "blocking_G_gt_U": bool(g[bl] > u[bl]),
                },
            },
        }

    res["R1"] = {
        "blocking": table,
        "theta": {a: arms[a]["theta"] for a in arms},
        "G_restoration": rest,
        "G_start_theta": arms["G"]["start_theta"],
        "G_B_tau_at_theta0": arms["G"].get("B_tau_at_theta0"),
    }
    res["R2"] = {"h": preds("h"), "h2": preds("h2")}
    res["diagnostics"] = diags
    res["runtime_seconds"] = time.perf_counter() - t_start
    figdir = (REPO / "figures") if out == DEFAULT_OUT.resolve() else out / "figures"
    label = f"stage 24 | commit {str(man.get('git_commit'))[:8]} | {man['timestamp_utc'][:19]}Z"
    try:
        plot_all(figdir, res, label)
    except Exception as exc:
        res["figures_error"] = f"{type(exc).__name__}: {exc}"
    summary = {"criteria": crit, "gate_passed": bool(all(crit.values())), "dev": bool(args.dev),
               "values": {**vals, "R1": res["R1"], "R2": res["R2"]},
               "runtime_seconds": res["runtime_seconds"]}  # fmt: skip
    dump()
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    print(json.dumps(crit, indent=2))
    print(json.dumps({k: {p: v[p]["confirmed"] for p in v} for k, v in res["R2"].items()}))
    print(f"runtime {res['runtime_seconds']:.1f} s")


if __name__ == "__main__":
    main()
