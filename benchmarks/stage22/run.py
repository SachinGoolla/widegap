"""Stage-22 benchmark (docs/benchmarks/stage22_adjoint_simplifications.md): Q1-Q3 (gating), R1-R2.

How much do common adjoint simplifications cost?  On the stage-20 device with d_body = 1.4 um
(stage-20 mesh rule) the exact adjoint E is compared with the variants A (drop F_mu mu_X), B (drop
d mu_bulk / dN in dF/dtheta), C (complete-ionization derivative), D = A + B + C, built by
``widegap.adjoint_variants`` from the SAME converged state as E, at theta = 0 and two seeded random
interior points.  P22b: SLSQP of the 17b-type problem (R_on,ch(h/2) <= R0, E_ox(400 V, h) <= E0,
17b classification constraints, bounds) started from theta = 0 with D gradients (exact function
values), compared with the E-gradient run (stage 21's main run if recorded, else run here).

Writes <out>/{results,summary,manifest}.json, <out>/stderr.log (default results/stage22) and
figures stage22_*.png (repo figures/ for the default output, else <out>/figures).  ``--dev`` marks
the run DEV (refused when --out lies under results/).

Interpretations fixed before any run (recorded in results.json["interpretations"]):
* Random points: ``np.random.default_rng(seed).uniform(-0.1, 0.1, 5)`` for seed 0 and 1, clipped
  to the 13d box; the linear drift-bound violation and the channel guard are recorded at each.
* Q1: relative Euclidean norm of (E - Richardson FD) over the five unit directions
  (eps = 1e-3 and 5e-4), for R_on,ch (h/2), the soft-min BV objective (h) and E_ox at 400 V (h).
* Q2: every variant gradient finite; two runs with fresh evaluators and the same call order
  identical bit for bit (gradients, values, hashes of the converged states); every assembly
  re-builds E through the variant machinery and compares it bitwise with the existing exact
  gradient (a mismatch raises).  Variants are assembled from the stored converged states of the
  single forward evaluation that produced E (no re-solve).
* Q3: gated on evaluations, as stage 17 J7: no failed evaluation, no non-converged Picard loop, no
  aborted SLSQP forward solve, every dense BV line bracketed and resolved, non-finite SuperLU calls
  only inside rejected (failed) Newton trial steps, no "On entry to" line in stderr.log.  Rejected
  Newton trial steps of the cold marches (step halving) are reported, not gated.
* A, B apply to R_on only (spec).  The off-state figures of merit have a Dirichlet electron row, so
  the mobility part of dF/dtheta vanishes identically (B is evaluated on them to record this); D
  equals C there.  N_A does not depend on theta, so dN_A^-/dN_A := 1 changes nothing.
* P22a: max over the three figures of merit of ||g_D - g_E|| / ||g_E|| at theta = 0 >= 1 %.
* P22b (i)-(iv) evaluated at the D-run final design: (i) |dense true-min BV(h) of D - that of the
  E run| >= 1 V; (ii) any of R_on,ch/R0 - 1, E_ox/E0 - 1, -(classification c_k), bound violation
  > 1e-6; (iii) KKT residual with E gradients (NNLS, widegap.adjoint_variants.kkt_residual; active
  = constraint value <= 1e-6, all constraint families + box) / ||grad BV_obj|| >= 0.1; (iv)
  SLSQP exit status != 0.  The E-run values are reported as baseline.
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

from widegap.adjoint_channel import ChannelRon, ClassificationConstraint
from widegap.adjoint_sic import richardson_fd
from widegap.adjoint_sic_corner import DenseCuts, dense_bv
from widegap.adjoint_sic_oxide import OxideAdjoint
from widegap.adjoint_variants import (
    VARIANTS,
    VariantComposite,
    VariantStats,
    fields_gradients,
    kkt_residual,
    objective_gradients,
    ron_gradients,
)
from widegap.bodyredesign import body_design
from widegap.channel_mobility import NoguchiParams
from widegap.manifest import build_manifest
from widegap.sic import SiC4H, TracedRegistry

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[2]
DEFAULT_OUT = REPO / "results" / "stage22"
DEFAULT_S21 = REPO / "results" / "stage21" / "results.json"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


S17B = _load("stage17b_run", REPO / "benchmarks/stage17b/run.py")  # fixed name: fork pickling
S17 = S17B.S17
M13D, M = S17.M13D, S17.M
DESIGN, TH0, Bounds = M13D.DESIGN, M13D.TH0, M13D.Bounds
tracked, job_diag = S17.tracked, S17.job_diag

D_BODY_UM = 1.4
DESIGN14 = body_design(DESIGN, D_BODY_UM)
LEVEL_ON, LEVEL_OFF, LEVEL_RECHECK = S17.LEVEL_ON, S17.LEVEL_OFF, S17.LEVEL_RECHECK
MARGIN = S17B.MARGIN
EPS_FD = M.EPS_FD
Q1_TOL = 1e-5
P22A_THRESH = 0.01
KKT_TOL, VIOL_TOL, KKT_RATIO, BV_DIFF_V = 1e-6, 1e-6, 0.1, 1.0
FOMS = ("R", "BV", "E_ox")
FOM_NAMES = {"R": "R_on,ch (h/2)", "BV": "soft-min BV (h)", "E_ox": "E_ox 400 V (h)"}
SEEDS_PT = (0, 1)
N_DENSE = M.N_DENSE
LOOSE_TOL_V = 1.0


def _lst(a):
    return np.asarray(a, dtype=float).tolist()


def design_points(bounds):
    pts = [("theta0", np.zeros(5), None)]
    lo = np.array([b[0] for b in bounds.box])
    hi = np.array([b[1] for b in bounds.box])
    for s in SEEDS_PT:
        th = np.clip(np.random.default_rng(s).uniform(-0.1, 0.1, 5), lo, hi)
        pts.append((f"seed{s}", th, s))
    return pts


# ============================================================ evaluators (1.4 um geometry)
def make_sc():
    return SiC4H(TracedRegistry.from_json())


def make_ev_off(with_cuts=True):
    return OxideAdjoint(make_sc(), DESIGN14, LEVEL_OFF, r_um=M.R_UM,
                        cuts_um=M.CUTS_UM if with_cuts else (), tau_v=M.TAU_V)  # fmt: skip


def make_ev_level(level, with_cuts=False):
    return OxideAdjoint(make_sc(), DESIGN14, level, r_um=M.R_UM,
                        cuts_um=M.CUTS_UM if with_cuts else (), tau_v=M.TAU_V)  # fmt: skip


def make_cr(level=LEVEL_ON):
    ev = make_ev_level(level)
    return ChannelRon(ev, NoguchiParams.from_registry(ev.sc.reg))


def make_cons():
    return ClassificationConstraint([make_ev_level(LEVEL_ON), make_ev_level(LEVEL_RECHECK)], MARGIN)


# ============================================================ jobs
@tracked
def job_grads(args):
    """All variants of all three figures of merit at ONE theta (fresh evaluators, fixed order)."""
    th, tag = args
    t0 = time.perf_counter()
    th = np.asarray(th, float)
    ev_off, cr = make_ev_off(), make_cr()
    st = VariantStats()
    o = objective_gradients(ev_off, th, ("E", "B", "C", "D"), st)
    f = fields_gradients(ev_off, th, ("E_ox",), ("E", "B", "C", "D"), st)
    r = ron_gradients(cr, th, VARIANTS, st)
    out = {
        "tag": tag, "theta": _lst(th),
        "values": {"R": r["R"], "BV": o["BV"], "E_ox": f["E_ox"]["S"]},
        "grad": {"R": {v: _lst(g) for v, g in r["grad"].items()},
                 "BV": {v: _lst(g) for v, g in o["grad"].items()},
                 "E_ox": {v: _lst(g) for v, g in f["E_ox"]["grad"].items()}},
        "state_sha": {"R": r["state_sha"], "BV": [_sha_state(c["state"].X) for c in o["res"]["outs"]],
                      "E_ox": _sha_state(f["X"])},
        "bv_k": _lst(o["bv_k"]), "softmin_w": _lst(o["w"]), "guard": r["guard"],
        "picard_iterations": r["out"]["iterations"], "variant_stats": st.as_dict(),
        "n_nodes_on": int(cr.ev.mesh.n_nodes), "n_nodes_off": int(ev_off.mesh.n_nodes),
    }  # fmt: skip
    out["diag"] = job_diag([cr], [ev_off])
    out["seconds"] = time.perf_counter() - t0
    return out


def _sha_state(x):
    import hashlib

    return hashlib.sha256(np.ascontiguousarray(x, dtype=float).tobytes()).hexdigest()[:16]


@tracked
def job_fd(args):
    """Richardson FD (eps = 1e-3, 5e-4) of one figure of merit along one unit direction."""
    fom, th, k = args
    t0 = time.perf_counter()
    th = np.asarray(th, float)
    u = np.zeros(5)
    u[k] = 1.0
    if fom == "R":
        cr = make_cr()
        fun = lambda t: cr.ron(t, grad=False)["R"]  # noqa: E731
        diag = lambda: job_diag([cr])  # noqa: E731
    else:
        ev = make_ev_off()
        if fom == "BV":
            fun = lambda t: ev.objective(t, grad=False)["BV"]  # noqa: E731
        else:
            fun = lambda t: ev.fields(t, grad=False, which=("E_ox",))["E_ox"]["S"]  # noqa: E731
        diag = lambda: job_diag([], [ev])  # noqa: E731
    rich, f1, f2 = richardson_fd(fun, th, u, EPS_FD)
    return {"fom": fom, "k": k, "richardson": float(rich), "fd_eps": float(f1),
            "fd_eps_half": float(f2), "spread": float(abs(f1 - f2)), "diag": diag(),
            "seconds": time.perf_counter() - t0}  # fmt: skip


@tracked
def job_slsqp(args):
    """17b-type SLSQP from theta = 0 on the 1.4 um geometry with E or D gradients."""
    variant, maxiter = args
    t0 = time.perf_counter()
    ev_off, cr = make_ev_off(), make_cr()
    cons = make_cons()
    bounds = Bounds(DESIGN14, ev_off.sc)
    comp = S17.Eval17(ev_off, cr) if variant == "E" else VariantComposite(ev_off, cr, variant)
    o0 = ev_off.objective(TH0, grad=False)
    f0 = ev_off.fields(TH0, grad=False)
    r0 = cr.ron(TH0, grad=False)["R"]
    bv0, e0 = o0["BV"], f0["E_ox"]["S"]
    th, res, pr, its = S17B.run_slsqp_b(comp, cons, TH0, bv0, r0, e0, 1.0, bounds, maxiter)
    out = {
        "variant": variant, "initial": {"bv_obj": bv0, "R": r0, "E_ox": e0},
        "theta": _lst(th), "status": int(res.status), "message": str(res.message),
        "nit": int(res.nit), "iterates": its, "evaluations": pr.history,
        "n_evals": pr.n_eval, "failure": pr.failure, "max_bound_violation": float(pr.max_violation),
        "n_trial_evals_crossing": int(sum(g["sign_changes"] > 0 for g in cr.guards)),
        "bitwise_checks": comp.stats.as_dict() if hasattr(comp, "stats") else None,
        "diag": job_diag([cr], [ev_off]), "seconds": time.perf_counter() - t0,
    }  # fmt: skip
    return out


@tracked
def job_final(args):
    """Everything evaluated at one design: values, constraints, E (and D) gradients, KKT with E
    gradients, dense true-minimum BV on h."""
    th, tag, r0, e0, with_dense = args
    t0 = time.perf_counter()
    th = np.asarray(th, float)
    ev_off, cr = make_ev_off(), make_cr()
    cons = make_cons()
    bounds = Bounds(DESIGN14, ev_off.sc)
    st = VariantStats()
    o = objective_gradients(ev_off, th, ("E", "D"), st)
    f = fields_gradients(ev_off, th, ("E_ox",), ("E", "D"), st)
    r = ron_gradients(cr, th, ("E", "D"), st)
    cv, cj = cons.values(th), cons.jac(th)
    lin_c = bounds.A @ th + bounds.b
    lo = np.array([b[0] for b in bounds.box])
    hi = np.array([b[1] for b in bounds.box])
    eye = np.eye(5)
    box_c = np.r_[th - lo, hi - th]
    box_j = np.vstack([eye, -eye])

    def kkt(grads, tol):
        g_bv, g_r, g_e = grads
        fam = [("R_on", np.array([1.0 - r["R"] / r0]), -g_r[None, :] / r0),
               ("E_ox", np.array([1.0 - f["E_ox"]["S"] / e0]), -g_e[None, :] / e0),
               ("linear_drift", lin_c, bounds.A), ("classification", cv, cj),
               ("box", box_c, box_j)]  # fmt: skip
        return kkt_residual(g_bv, fam, tol)

    gE = (o["grad"]["E"], r["grad"]["E"], f["E_ox"]["grad"]["E"])
    gD = (o["grad"]["D"], r["grad"]["D"], f["E_ox"]["grad"]["D"])
    out = {
        "tag": tag, "theta": _lst(th),
        "values": {"bv_obj": o["BV"], "R": r["R"], "E_ox": f["E_ox"]["S"]},
        "violations": {"R_rel": r["R"] / r0 - 1.0, "E_ox_rel": f["E_ox"]["S"] / e0 - 1.0,
                       "classification_min_c": float(cv.min()),
                       "classification_violation": float(max(0.0, -cv.min())),
                       "linear_and_box": float(bounds.violation(th))},
        "constraint_values": {"classification": _lst(cv), "linear_drift_min": float(lin_c.min())},
        "grads_E": {"bv_obj": _lst(gE[0]), "R": _lst(gE[1]), "E_ox": _lst(gE[2])},
        "grads_D": {"bv_obj": _lst(gD[0]), "R": _lst(gD[1]), "E_ox": _lst(gD[2])},
        "kkt_E": {f"{t:g}": kkt(gE, t) for t in (KKT_TOL, 1e-4, 1e-8)},
        "kkt_D_own_gradients": kkt(gD, KKT_TOL), "variant_stats": st.as_dict(),
        "guard": r["guard"],
    }  # fmt: skip
    out["diag"] = job_diag([cr], [ev_off])
    if with_dense:
        ev_d = make_ev_level(LEVEL_OFF)
        xs = np.linspace(0.0, float(ev_d.sc.length_cm(ev_d.prm.x_max)) * 1e4, N_DENSE)
        d = dense_bv(ev_d, th, DenseCuts(ev_d, xs), tol_v=0.02, loose_tol_v=LOOSE_TOL_V)
        bv = np.asarray(d["bv"], float)
        out["dense"] = {"bv_min": float(np.nanmin(bv)), "x_argmin_um": float(xs[np.nanargmin(bv)]),
                        "bv": bv.tolist(), "x_um": xs.tolist(), "n_unbracketed": int(np.isnan(bv).sum()),
                        "resolved_all": bool(d["resolved_all"])}  # fmt: skip
        out["diag"]["off_dense"] = dict(ev_d.stats)
    out["seconds"] = time.perf_counter() - t0
    return out


# ============================================================ stage-21 reuse
def load_stage21(path: Path):
    """The recorded stage-21 main (E-gradient) run, or None if absent / incomplete / DEV."""
    try:
        d = json.loads(Path(path).read_text())
        if "runtime_seconds" not in d or d.get("dev", True):
            return None
        sl = d["slsqp_main"]
        return {"variant": "E", "source": str(path), "theta": d["final"]["theta"],
                "status": int(sl["status"]), "message": str(sl["message"]), "nit": int(sl["nit"]),
                "iterates": sl.get("iterates"), "evaluations": sl.get("evaluations"),
                "initial": d.get("initial"), "failure": None}  # fmt: skip
    except (OSError, KeyError, ValueError, TypeError):
        return None


# ============================================================ analysis helpers
def rel_err(g, ref):
    g, ref = np.asarray(g, float), np.asarray(ref, float)
    return float(np.linalg.norm(g - ref) / np.linalg.norm(ref))


def r1_entry(gv, ge):
    gv, ge = np.asarray(gv, float), np.asarray(ge, float)
    d = gv - ge
    return {
        "rel_norm_error": float(np.linalg.norm(d) / np.linalg.norm(ge)),
        "cosine": float(gv @ ge / (np.linalg.norm(gv) * np.linalg.norm(ge))),
        "component_error_over_norm_E": (d / np.linalg.norm(ge)).tolist(),
        "component_rel_error": (d / ge).tolist(),
        "grad": gv.tolist(), "norm": float(np.linalg.norm(gv)),
    }  # fmt: skip


# ============================================================ figures
def plot_all(figdir, res, label):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figdir = Path(figdir)
    figdir.mkdir(parents=True, exist_ok=True)
    pts = list(res["r1"].keys())
    fig, axs = plt.subplots(1, 3, figsize=(12.0, 4.0), constrained_layout=True, sharey=True)
    cols = {"A": "#0072B2", "B": "#009E73", "C": "#E69F00", "D": "#D55E00"}
    for ax, fom in zip(axs, FOMS, strict=True):
        vs = [v for v in ("A", "B", "C", "D") if v in res["r1"][pts[0]][fom]]
        w = 0.8 / len(pts)
        for j, p in enumerate(pts):
            ax.bar(np.arange(len(vs)) + j * w, [max(res["r1"][p][fom][v]["rel_norm_error"], 1e-12)
                                                for v in vs], w, label=p,
                   color=[cols[v] for v in vs], alpha=0.45 + 0.25 * j)  # fmt: skip
        ax.set_xticks(np.arange(len(vs)) + 0.3, vs)
        ax.set_yscale("log")
        ax.axhline(P22A_THRESH, color="k", ls="--", lw=0.8)
        ax.set_title(FOM_NAMES[fom], fontsize=9)
        ax.grid(alpha=0.3, axis="y")
    axs[0].set_ylabel("||g_variant - g_E|| / ||g_E||")
    fig.text(0.5, 0.005, label + " | bars: theta0, seed0, seed1 (light to dark)", ha="center",
             fontsize=7, color="0.4")  # fmt: skip
    fig.savefig(figdir / "stage22_gradient_error.png", dpi=130)
    plt.close(fig)
    runs = res["p22b"]["runs"]
    fig, axs = plt.subplots(1, 3, figsize=(12.0, 3.8), constrained_layout=True)
    for tag, c in (("E", "#0072B2"), ("D", "#D55E00")):
        run = runs.get(tag)
        if not run or not run.get("evaluations"):
            continue
        ev = run["evaluations"]
        x = np.arange(len(ev))
        axs[0].plot(x, [e["bv_obj"] for e in ev], "o-", color=c, ms=3, label=f"{tag} gradients")
        axs[1].plot(x, [e["R"] / ev[0]["R"] for e in ev], "o-", color=c, ms=3)
        axs[2].plot(x, [e["E_ox"] / ev[0]["E_ox"] for e in ev], "o-", color=c, ms=3)
    for ax, t in zip(axs, ("BV_obj [V]", "R_on,ch / R0", "E_ox / E0"), strict=True):
        ax.set_xlabel("function evaluation")
        ax.set_ylabel(t)
        ax.grid(alpha=0.3)
    axs[0].legend(fontsize=8)
    fig.text(0.5, 0.005, label, ha="center", fontsize=7, color="0.4")
    fig.savefig(figdir / "stage22_slsqp_trajectories.png", dpi=130)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(7.2, 4.0), constrained_layout=True)
    x = np.arange(5)
    for j, (tag, c) in enumerate((("E", "#0072B2"), ("D", "#D55E00"))):
        run = runs.get(tag)
        if run:
            ax.bar(x + (j - 0.5) * 0.4, run["theta"], 0.4, color=c, label=f"{tag} gradients")
    ax.set_xticks(x, [f"theta_{i}" for i in range(4)] + ["theta_J"])
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")
    fig.text(0.5, 0.005, label, ha="center", fontsize=7, color="0.4")
    fig.savefig(figdir / "stage22_final_theta.png", dpi=130)
    plt.close(fig)


# ============================================================ main
def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--stage21-results", type=Path, default=DEFAULT_S21)
    ap.add_argument("--no-reuse", action="store_true", help="always run the E-gradient SLSQP here")
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
    man = build_manifest(REPO, sys.argv, [0, 1])
    man["dev"] = bool(args.dev)
    (out / "manifest.json").write_text(json.dumps(man, indent=2) + "\n")
    t_start = time.perf_counter()
    sc = make_sc()
    bounds = Bounds(DESIGN14, sc)
    pts = design_points(bounds)
    s21 = None if args.no_reuse else load_stage21(args.stage21_results)
    res: dict = {"dev": bool(args.dev)}
    res["interpretations"] = [ln.strip() for ln in __doc__.split("Interpretations")[1].splitlines()[1:]
                              if ln.strip()]  # fmt: skip
    res["config"] = {
        "d_body_um": D_BODY_UM, "mesh_kw": dict(DESIGN14.mesh_kw), "levels": {"on": LEVEL_ON,
        "off": LEVEL_OFF, "recheck": LEVEL_RECHECK}, "eps_fd": EPS_FD, "workers": args.workers,
        "max_iter": args.maxiter, "margin_m": MARGIN, "q1_tol": Q1_TOL, "p22a_threshold": P22A_THRESH,
        "p22b": {"bv_diff_V": BV_DIFF_V, "violation_tol": VIOL_TOL, "kkt_ratio": KKT_RATIO,
                 "kkt_active_tol": KKT_TOL},
        "design_points": [{"tag": t, "theta": th.tolist(), "seed": s,
                           "bound_violation": bounds.violation(th)} for t, th, s in pts],
        "stage21_reuse": None if s21 is None else s21["source"],
        "bounds_identical_to_0.8um": bool(
            np.array_equal(bounds.A, Bounds(DESIGN, sc).A)
            and np.array_equal(bounds.b, Bounds(DESIGN, sc).b)
        ),
        "registry": S17.registry_check(),
    }  # fmt: skip
    (out / "results.json").write_text(json.dumps(res, indent=2, default=float) + "\n")

    pool = ProcessPoolExecutor(args.workers, mp_context=mp.get_context("fork"))
    fut: dict = {}

    def submit(key, fn, a):
        fut[pool.submit(fn, a)] = key

    submit(("slsqp", "D"), job_slsqp, ("D", args.maxiter))
    if s21 is None:
        submit(("slsqp", "E"), job_slsqp, ("E", args.maxiter))
    for tag, th, _ in pts:
        submit(("grads", tag, 1), job_grads, (th, tag))
        submit(("grads", tag, 2), job_grads, (th, tag))
        for fom in FOMS:
            for k in range(5):
                submit(("fd", tag, fom, k), job_fd, (fom, th, k))
    jobs: dict = {"slsqp": {}, "grads": {}, "fd": {}, "final": {}}
    diags: dict = {}
    errors: list = []
    pending = set(fut)
    sub_done: set = set()
    if s21 is not None:
        jobs["slsqp"]["E"] = s21
    t0 = time.perf_counter()

    def maybe_submit_finals():
        i0 = (jobs["slsqp"].get("D") or {}).get("initial") or (s21 or {}).get("initial")
        if i0 is None:
            return
        r0, e0 = i0.get("R", i0.get("R_channel_h2")), i0.get("E_ox")
        if r0 is None or e0 is None:
            return
        for tag in ("D", "E"):
            run_ = jobs["slsqp"].get(tag)
            if run_ is not None and "theta" in run_ and tag not in sub_done:
                sub_done.add(tag)
                nf = pool.submit(job_final, (run_["theta"], f"final_{tag}", r0, e0, True))
                fut[nf] = ("final", f"final_{tag}")
                pending.add(nf)
        if "init" not in sub_done:
            sub_done.add("init")
            nf = pool.submit(job_final, (TH0, "initial", r0, e0, True))
            fut[nf] = ("final", "initial")
            pending.add(nf)

    maybe_submit_finals()
    while pending:
        for f in as_completed(list(pending)):
            pending.discard(f)
            key = fut[f]
            try:
                r = f.result()
            except Exception as exc:
                r = {"error": f"{type(exc).__name__}: {exc}", "seconds": float("nan")}
                errors.append({"job": list(map(str, key)), "error": r["error"]})
            kind = key[0]
            if kind == "slsqp":
                jobs["slsqp"][key[1]] = r
            elif kind == "grads":
                jobs["grads"].setdefault(key[1], {})[key[2]] = r
            elif kind == "fd":
                jobs["fd"].setdefault(key[1], {}).setdefault(key[2], {})[key[3]] = r
            else:
                jobs["final"][key[1]] = r
            if "diag" in r:
                diags["|".join(map(str, key))] = r["diag"]
            print(f"done {key} {r.get('seconds', float('nan')):.0f}s"
                  f"{' ERROR ' + r['error'] if 'error' in r else ''}", flush=True)  # fmt: skip
            if kind == "slsqp" and "error" not in r:
                print(f"  slsqp {key[1]}: status {r['status']} nit {r['nit']} {r['message']}",
                      flush=True)  # fmt: skip
            maybe_submit_finals()
            break
    pool.shutdown()
    T = {"parallel": time.perf_counter() - t0}
    res["jobs"] = {"slsqp": jobs["slsqp"], "grads": jobs["grads"], "final": jobs["final"]}

    # ================================================================ Q1, R1, Q2
    crit, vals = {}, {}
    q1, r1, q2_rows = {}, {}, []
    finite_all, det_all, sha_all = True, True, True
    for tag, _th, _ in pts:
        g1, g2 = jobs["grads"][tag][1], jobs["grads"][tag][2]
        if "error" in g1 or "error" in g2:
            continue
        q1[tag], r1[tag] = {}, {}
        for fom in FOMS:
            fd = np.array([jobs["fd"][tag][fom][k]["richardson"] for k in range(5)])
            spread = max(jobs["fd"][tag][fom][k]["spread"] for k in range(5))
            ge = np.array(g1["grad"][fom]["E"])
            q1[tag][fom] = {"rel_norm_error": rel_err(ge, fd), "fd_grad": fd.tolist(),
                            "max_fd_spread": float(spread), "E_grad": ge.tolist()}  # fmt: skip
            r1[tag][fom] = {}
            for v in VARIANTS[1:]:
                if v in g1["grad"][fom]:
                    r1[tag][fom][v] = r1_entry(g1["grad"][fom][v], ge)
            finite_all &= all(np.all(np.isfinite(g)) for g in g1["grad"][fom].values())
        det = all(
            np.array_equal(np.asarray(g1["grad"][f_][v]), np.asarray(g2["grad"][f_][v]))
            for f_ in FOMS for v in g1["grad"][f_]
        ) and g1["values"] == g2["values"]  # fmt: skip
        sha = g1["state_sha"] == g2["state_sha"]
        det_all &= det
        sha_all &= sha
        q2_rows.append({"tag": tag, "bitwise_identical_rerun": det, "same_states": sha,
                        "variant_stats": g1["variant_stats"]})  # fmt: skip
    q1_max = max((v["rel_norm_error"] for p in q1.values() for v in p.values()), default=np.inf)
    crit["Q1"] = bool(q1 and len(q1) == len(pts) and q1_max <= Q1_TOL)
    vals["Q1"] = {"max_rel_norm_error": q1_max, "tolerance": Q1_TOL, "by_point": q1}
    res["r1"], res["q1"] = r1, q1
    runs = jobs["slsqp"]
    d_stats = (runs.get("D") or {}).get("bitwise_checks") or {}
    n_chk = sum(r["variant_stats"]["n_bitwise_checks"] for r in q2_rows) + d_stats.get(
        "n_bitwise_checks", 0)  # fmt: skip
    crit["Q2"] = bool(len(q2_rows) == len(pts) and finite_all and det_all and sha_all and n_chk > 0
                      and "error" not in runs.get("D", {"error": 1}))  # fmt: skip
    vals["Q2"] = {"finite": bool(finite_all), "bitwise_identical_reruns": bool(det_all),
                  "same_converged_states": bool(sha_all), "bitwise_E_checks_total": n_chk,
                  "D_run_bitwise_E_checks": d_stats, "rows": q2_rows}  # fmt: skip

    # ================================================================ P22a
    p22a = {}
    if "theta0" in r1:
        p22a = {f: r1["theta0"][f]["D"]["rel_norm_error"] for f in FOMS}
    p22a_holds = bool(p22a and max(p22a.values()) >= P22A_THRESH)
    res["p22a"] = {"D_vs_E_rel_norm_theta0": p22a, "threshold": P22A_THRESH,
                   "prediction_holds": p22a_holds}  # fmt: skip

    # ================================================================ P22b
    fin = jobs["final"]

    def kk(tag, tol=KKT_TOL):
        return fin[tag]["kkt_E"][f"{tol:g}"] if tag in fin and "kkt_E" in fin[tag] else None

    p22b = {"runs": {}}
    for tag in ("D", "E"):
        rn = runs.get(tag)
        if rn and "error" not in rn:
            p22b["runs"][tag] = {k: rn.get(k) for k in ("theta", "status", "message", "nit",
                                 "iterates", "evaluations", "n_evals", "failure", "source",
                                 "max_bound_violation", "initial", "seconds")}  # fmt: skip
            f_ = fin.get(f"final_{tag}")
            if f_ and "error" not in f_:
                p22b["runs"][tag].update(
                    {"final": {k: f_[k] for k in ("values", "violations", "constraint_values",
                                                  "kkt_E", "kkt_D_own_gradients", "grads_E")},
                     "dense_bv_min": f_["dense"]["bv_min"], "dense_x_argmin_um": f_["dense"]["x_argmin_um"]})  # fmt: skip
    if "D" in p22b["runs"] and "E" in p22b["runs"] and "final" in p22b["runs"]["D"] \
            and "final" in p22b["runs"]["E"]:  # fmt: skip
        d, e = p22b["runs"]["D"], p22b["runs"]["E"]
        viol = d["final"]["violations"]
        worst = max(viol["R_rel"], viol["E_ox_rel"], viol["classification_violation"],
                    viol["linear_and_box"])  # fmt: skip
        cond = {
            "i_dense_bv_diff_ge_1V": abs(d["dense_bv_min"] - e["dense_bv_min"]) >= BV_DIFF_V,
            "ii_constraint_violated_gt_1e-6": bool(worst > VIOL_TOL),
            "iii_kkt_E_ratio_ge_0.1": bool(d["final"]["kkt_E"][f"{KKT_TOL:g}"]["ratio"] >= KKT_RATIO),
            "iv_exit_status_nonzero": bool(d["status"] != 0),
        }  # fmt: skip
        p22b.update({
            "conditions": cond, "prediction_holds": bool(any(cond.values())),
            "dense_bv_min_D": d["dense_bv_min"], "dense_bv_min_E": e["dense_bv_min"],
            "dense_bv_diff_D_minus_E": d["dense_bv_min"] - e["dense_bv_min"],
            "worst_violation_D": worst, "kkt_ratio_D": d["final"]["kkt_E"][f"{KKT_TOL:g}"]["ratio"],
            "kkt_ratio_E_baseline": e["final"]["kkt_E"][f"{KKT_TOL:g}"]["ratio"],
            "theta_diff_norm": float(np.linalg.norm(np.array(d["theta"]) - np.array(e["theta"]))),
            "exit_status_D": d["status"], "exit_status_E": e["status"],
            "E_run_source": e.get("source") or "run here (identical settings)",
        })  # fmt: skip
    res["p22b"] = p22b
    if "initial" in fin and "dense" in fin["initial"]:
        res["initial_dense_bv_min"] = fin["initial"]["dense"]["bv_min"]

    # ================================================================ Q3
    def agg(a, b):
        return sum(d[a].get(b, 0) for d in diags.values() if a in d)

    stderr = log_path.read_text(errors="replace") if log_path.exists() else ""
    entry = [ln for ln in stderr.splitlines() if "On entry to" in ln]
    nf_tot = sum(d["nonfinite_calls"] for d in diags.values())
    nf_failed = sum(d["NF"]["failed_newton"] for d in diags.values())
    nf_adj = sum(d["NF"]["adjoint"] for d in diags.values())
    nf_conv = sum(d["NF"]["converged_newton"] for d in diags.values())
    dense_ok = all(f_["dense"]["resolved_all"] and f_["dense"]["n_unbracketed"] == 0
                   for f_ in fin.values() if "dense" in f_)  # fmt: skip
    aborted = [t for t, r in runs.items() if r.get("failure")]
    failed_dense = sum(f_["diag"].get("off_dense", {}).get("failed_evals", 0) for f_ in fin.values()
                       if "diag" in f_)  # fmt: skip
    crit["Q3"] = bool(
        not errors and not entry and not aborted and dense_ok and nf_tot == nf_failed
        and nf_adj == 0 and nf_conv == 0 and agg("on_state", "failed_picard") == 0
        and agg("on_state", "failed_evals") == 0 and agg("off_state", "failed_evals") == 0
        and failed_dense == 0 and len(fin) >= 2 and all(r.get("n_evals", 1) > 0 for r in runs.values())
    )  # fmt: skip
    vals["Q3"] = {"job_errors": errors, "on_entry_to_lines": len(entry), "aborted_slsqp": aborted,
                  "dense_all_resolved": bool(dense_ok), "nonfinite_calls_total": nf_tot,
                  "nonfinite_in_failed_newton": nf_failed, "nonfinite_in_adjoint": nf_adj,
                  "nonfinite_in_converged_newton": nf_conv,
                  "failed_picard": agg("on_state", "failed_picard"),
                  "failed_evals_on": agg("on_state", "failed_evals"),
                  "failed_evals_off": agg("off_state", "failed_evals"),
                  "failed_newton_attempts_reported": agg("off_state", "failed_attempts")
                  + agg("on_state", "failed_attempts_on"),
                  "guards_ok": bool(all(d["guards"]["all_ok"] for d in diags.values()
                                        if "guards" in d))}  # fmt: skip
    res["diagnostics"], res["seconds"] = diags, T
    res["runtime_seconds"] = time.perf_counter() - t_start
    figdir = (REPO / "figures") if out == DEFAULT_OUT.resolve() else out / "figures"
    label = f"stage 22 | commit {str(man.get('git_commit'))[:8]} | {man['timestamp_utc'][:19]}Z"
    try:
        if r1:
            plot_all(figdir, res, label)
    except Exception as exc:
        res["figures_error"] = f"{type(exc).__name__}: {exc}"
    r1_table = {t: {f: {v: {"rel_norm_error": e["rel_norm_error"], "cosine": e["cosine"],
                            "component_error_over_norm_E": e["component_error_over_norm_E"]}
                        for v, e in fv.items()} for f, fv in pt.items()} for t, pt in r1.items()}  # fmt: skip
    summary = {"criteria": crit, "gate_passed": bool(all(crit.values())), "dev": bool(args.dev),
               "values": vals, "R1": r1_table, "P22a": res["p22a"],
               "P22b": {k: v for k, v in p22b.items() if k != "runs"}
               | {"runs": {t: {k: v for k, v in r.items() if k not in ("iterates", "evaluations")}
                            for t, r in p22b["runs"].items()}},
               "seconds": T, "runtime_seconds": res["runtime_seconds"]}  # fmt: skip
    (out / "results.json").write_text(json.dumps(res, indent=2, default=float) + "\n")
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    print(json.dumps(crit, indent=2))
    print(f"runtime {res['runtime_seconds']:.1f} s")


if __name__ == "__main__":
    main()
