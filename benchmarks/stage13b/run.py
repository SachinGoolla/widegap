"""Stage-13b benchmark (docs/benchmarks/stage13b_corner_aware.md): U1-U7 (gating), R1-R2.

Corner-aware adjoint optimization of the stage-12 4H-SiC DMOSFET with rounded (erfc-product)
p-body and n+ source: maximize the soft-minimum BV_obj = -tau log sum_k exp(-BV_k/tau) over K
fixed vertical cuts (each BV_k the Niwa electron-initiated ionization-integral breakdown) subject
to R_on,sp <= R_on,sp,0 (V_GS = 15 V, V_DS = 0.1 V), SLSQP on adjoint gradients
(widegap.adjoint_sic_corner on top of widegap.adjoint_sic).

Writes <out>/{results,summary,manifest}.json (default results/stage13b) and the figures
stage13b_per_cut_bv.png, stage13b_optimized_profile.png, stage13b_sensitivity_maps.png
(repo figures/ for the default output, else <out>/figures).

Study choices, fixed BEFORE the official/dev run and recorded in results.json:

* device: the stage-12 12b design (benchmarks/stage12/run.py), mesh level 0 (h) for the
  optimization, level 1 (h/2) and level 2 (h/4) for verification.
* rounding length r = 0.1 um (standard deviation of the erfc edge; the spec's example value).
  History, disclosed: a first dev run used r = 0.2 um (chosen from a resolvability scan of the
  initial design only); the on-state solution then proved path dependent (R_on differed 1.7e-4
  between cold starts: the Gaussian source tail leaves a thin, compensated p layer under the
  n+ source so the hole system of the channel region floats) which made the R_on Richardson check
  meaningless (U2 20% error).  r = 0.05, 0.1, 0.15 are deterministic; r = 0.1 um was therefore
  used.  The cut positions below were selected on the r = 0.2 initial BV(x) profile and were
  NOT re-selected for r = 0.1 (so they straddle the sharper r = 0.1 dip at x ~ 0.85 um).
* K = 10 objective cuts at level-0 mesh columns (so the same physical lines exist on h/2, h/4):
  x = 0, 0.4, 0.667, 0.8, 0.933, 1.067, 1.2, 1.6, 2.249, 3.0 um, i.e. from the JFET symmetry line
  through the JFET/body corner region (the initial minimum is at x ~ 0.93 um) to the planar region
  under the body contact (x_bc = 3.0 um, the stage-13 cut).  tau = 5 V.
* U5 dense set: 49 equispaced vertical lines on [0, x_max] (bilinear interpolation of the mesh
  potential, same physical positions on h and h/2), independent of the objective cuts.
* theta bounds and Pareto levels (R_on <= {1, 1.05, 1.10, 1.20} R_on,0) as stage 13; the three
  looser levels start from the main solution (independent runs, executed in parallel).
* U1/U6 reference voltage: U1 uses the initial design's true (dense-set) BV on mesh h as a fixed
  physical drain voltage for h, h/2, h/4; U6 uses the initial design's own true BV on each mesh.
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
from types import SimpleNamespace

import numpy as np
from scipy.optimize import minimize

from widegap.adjoint_sic import EvalError, richardson_fd
from widegap.adjoint_sic_corner import (
    CornerAdjoint,
    DenseCuts,
    dense_bv,
    solve_at_voltage,
)
from widegap.manifest import build_manifest
from widegap.sic import SiC4H, TracedRegistry

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[2]
DEFAULT_OUT = REPO / "results" / "stage13b"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


S13 = _load("stage13_run", REPO / "benchmarks/stage13/run.py")
DESIGN, FLAGS12 = S13.DESIGN, S13.FLAGS12
Bounds, directions, unipolar_limit, ideal_r_at = (
    S13.Bounds, S13.directions, S13.unipolar_limit, S13.ideal_r_at,
)  # fmt: skip

# ---- study parameters (docstring) -----------------------------------------------------------
R_UM = 0.1
CUTS_UM = (0.0, 0.4, 0.667, 0.8, 0.933, 1.067, 1.2, 1.6, 2.249, 3.0)
TAU_V = 5.0
N_DENSE = 49
PARETO_RATIOS = (1.05, 1.10, 1.20)
SEEDS = (0, 1, 2)
EPS_FD = 1e-3
FEAS_TOL = 1e-8
MAX_ITER = 60
TH0 = np.zeros(5)
U1_LEVELS = (0, 1, 2)
U1_TOL = 0.02
U2_TOL_BV, U2_TOL_R = 1e-4, 1e-5
LOOSE_TOL_V = {0: 1.0, 1: 3.0}


def rss_mib():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def make_ev(level, with_cuts=True):
    sc = SiC4H(TracedRegistry.from_json())
    return CornerAdjoint(sc, DESIGN, level, r_um=R_UM, cuts_um=CUTS_UM if with_cuts else (),
                         tau_v=TAU_V)  # fmt: skip


def dense_positions(sc, ev):
    x_max_um = float(sc.length_cm(ev.prm.x_max)) * 1e4
    return np.linspace(0.0, x_max_um, N_DENSE)


# ============================================================ optimization problem
class EvalAbort(RuntimeError):
    """A forward solve failed inside SLSQP; the run is stopped and reported."""


class Problem:
    failure = None

    def __init__(self, ev: CornerAdjoint, bv0, r0, ratio, bounds: Bounds):
        self.ev, self.bv0, self.r0, self.ratio, self.bounds = ev, bv0, r0, ratio, bounds
        self.last = (None, None)
        self.max_violation = -np.inf
        self.n_eval = 0
        self.history: list = []

    def _get(self, th):
        th = np.asarray(th, float)
        key = th.tobytes()
        if self.last[0] != key:
            self.max_violation = max(self.max_violation, self.bounds.violation(th))
            try:
                b = self.ev.objective(th)
                r = self.ev.ron(th)
            except EvalError as exc:  # a line-search probe whose forward solve does not converge
                self.failure = {"theta": th.tolist(), "error": str(exc), "n_eval": self.n_eval}
                raise EvalAbort(str(exc)) from exc
            self.n_eval += 1
            self.history.append({"theta": th.tolist(), "bv_obj": b["BV"], "bv_k": b["bv_k"].tolist(),
                                 "R": r["R"]})  # fmt: skip
            self.last = (key, (b, r))
        return self.last[1]

    def fun(self, th):
        return -self._get(th)[0]["BV"] / self.bv0

    def jac(self, th):
        return -self._get(th)[0]["grad"] / self.bv0

    def c_fun(self, th):
        return self.ratio - self._get(th)[1]["R"] / self.r0

    def c_jac(self, th):
        return -self._get(th)[1]["grad"] / self.r0


def run_slsqp(ev, th_start, bv0, r0, ratio, bounds, maxiter=MAX_ITER, ftol=1e-10):
    pr = Problem(ev, bv0, r0, ratio, bounds)
    cons = [
        {"type": "ineq", "fun": pr.c_fun, "jac": pr.c_jac},
        {"type": "ineq", "fun": lambda t: bounds.A @ t + bounds.b, "jac": lambda t: bounds.A},
    ]
    its: list = []
    try:
        res = minimize(
            pr.fun, np.asarray(th_start, float), jac=pr.jac, method="SLSQP", bounds=bounds.box,
            constraints=cons, options={"maxiter": maxiter, "ftol": ftol},
            callback=lambda xk: its.append(np.asarray(xk).tolist()),
        )  # fmt: skip
        return np.asarray(res.x, float), res, pr, its
    except EvalAbort as exc:
        # keep the best feasible evaluated iterate (strictly feasible, bound violation <= 1e-9)
        best, bv_best = np.asarray(th_start, float), -np.inf
        for h in pr.history:
            t = np.asarray(h["theta"])
            if (h["R"] <= ratio * r0 * (1.0 + FEAS_TOL) and bounds.violation(t) <= 1e-9
                    and h["bv_obj"] > bv_best):  # fmt: skip
                best, bv_best = t, h["bv_obj"]
        res = SimpleNamespace(status=-1, message=f"ABORTED (forward solve failed): {exc}",
                              nit=len(its), nfev=pr.n_eval, njev=pr.n_eval, x=best)  # fmt: skip
        return best, res, pr, its


# ============================================================ gradient checks (U2, U3)
def grad_check(ev: CornerAdjoint, th):
    """Adjoint vs Richardson central FD (eps = 1e-3, 5e-4) of BV_obj [V] and R_on [ohm cm^2] along
    3 seeded directions."""
    th = np.asarray(th, float)
    b = ev.objective(th)
    r = ev.ron(th)
    gb, gr, r0 = b["grad"], r["grad"], r["R"]
    rows = []
    for s, u in zip(SEEDS, directions(), strict=True):
        rec = {"seed": s}
        for name, fun, g, f0 in (
            ("bv", lambda t: ev.objective(t, grad=False)["BV"], gb, b["BV"]),
            ("ron", lambda t: ev.ron(t, grad=False)["R"], gr, r["R"]),
        ):
            rich, f1, f2 = richardson_fd(fun, th, u, EPS_FD)
            a = float(g @ u)
            rec[name] = {
                "adjoint": a, "richardson": rich, "fd_eps": f1, "fd_eps_half": f2,
                "abs_diff": abs(a - rich), "rel_err": abs(a - rich) / abs(rich),
                "fd_spread": abs(f1 - f2),
                "roundoff_floor_est": abs(f0) * 1e-13 / (EPS_FD / 2.0),
            }  # fmt: skip
        rec["ron_norm"] = {
            "adjoint": rec["ron"]["adjoint"] / r0, "richardson": rec["ron"]["richardson"] / r0,
            "abs_diff": rec["ron"]["abs_diff"] / r0,
        }  # fmt: skip
        rows.append(rec)
    ev.objective(th)
    ev.ron(th)
    return {
        "bv_obj_V": b["BV"], "bv_k": b["bv_k"].tolist(), "weights": b["w"].tolist(),
        "R_ohm_cm2": r0, "grad_bv_norm": float(np.linalg.norm(gb)),
        "grad_ron_norm": float(np.linalg.norm(gr)),
        "grad_ron_over_R0_norm": float(np.linalg.norm(gr) / r0),
        "grad_bv": gb.tolist(), "grad_ron": gr.tolist(), "directions": rows,
    }  # fmt: skip


# ============================================================ parallel jobs (module level)
def job_dense(args):
    """Dense BV profile (+ objective-cut columns) of one design on one mesh, fresh solves."""
    level, theta, tag, with_ron = args
    t0 = time.perf_counter()
    ev = make_ev(level, with_cuts=False)
    sc = ev.sc
    xs = dense_positions(sc, ev)
    xs_all = np.r_[xs, np.asarray(CUTS_UM)]
    dc = DenseCuts(ev, xs_all)
    th = np.asarray(theta, float)
    d = dense_bv(ev, th, dc, tol_v=0.02, loose_tol_v=LOOSE_TOL_V.get(level, 3.0))
    bv = d["bv"]
    bd = bv[:N_DENSE]
    jm = int(np.nanargmin(bd))
    out = {
        "tag": tag, "level": level, "x_um": xs.tolist(), "bv": bd.tolist(),
        "cut_x_um": list(CUTS_UM), "bv_at_objective_columns": bv[N_DENSE:].tolist(),
        "bv_min": float(bd[jm]), "argmin_index": jm, "x_argmin_um": float(xs[jm]),
        "n_unbracketed": int(np.isnan(bd).sum()), "resolved_all": d["resolved_all"],
        "n_refine_solves": d["n_refine_solves"], "stats": dict(ev.stats),
    }  # fmt: skip
    if with_ron:
        out["R"] = ev.ron(th, grad=False)["R"]
    out["seconds"] = time.perf_counter() - t0
    return out


def job_peak(args):
    """Global peak field (blocking edges), and the source-junction-box peak, at drain voltage V."""
    level, theta, tag, v_volt = args
    t0 = time.perf_counter()
    ev = make_ev(level, with_cuts=False)
    sc = ev.sc
    th = np.asarray(theta, float)
    X = solve_at_voltage(ev, th, v_volt)
    psi = X[0::3]
    pk = ev.base.peak_field(psi)
    mod, mesh = ev.mod, ev.mesh
    e = np.abs(psi[mod.b] - psi[mod.a]) / mesh.h
    pr_, m3 = (
        ev.prm,
        3.0 * ev.base.r_s,
    )  # n+ source junction region (3 r margin), semiconductor edges
    mx = 0.5 * (mesh.node_xy[mod.a, 0] + mesh.node_xy[mod.b, 0])
    my = 0.5 * (mesh.node_xy[mod.a, 1] + mesh.node_xy[mod.b, 1])
    box = ((ev.dev.reg.w_s > 0) & (mx >= pr_.x_src0 - m3) & (mx <= pr_.x_src1 + m3)
           & (my >= -pr_.d_src - m3) & (my <= 0.0))  # fmt: skip
    ebox = float(np.max(np.where(box, e, 0.0)))
    k = int(np.argmax(np.where(ev.dev.blocking_edge_mask(), e, 0.0)))
    mid = 0.5 * (mesh.node_xy[mod.a[k]] + mesh.node_xy[mod.b[k]])
    f = lambda v: float(sc.field_phys(v)) * 1e-6  # noqa: E731
    um = lambda v: float(sc.length_cm(v)) * 1e4  # noqa: E731
    return {
        "tag": tag, "level": level, "v_volt": v_volt, "E_peak_MV_cm": f(pk[0]),
        "x_um": um(pk[1]), "y_um": um(pk[2]), "dir": pk[3],
        "E_source_region_peak_MV_cm": f(ebox), "n_nodes": int(mesh.n_nodes),
        "seconds": time.perf_counter() - t0, "mid_check": [um(mid[0]), um(mid[1])],
    }  # fmt: skip


def job_pareto(args):
    """Independent SLSQP run (level 0) from ``th_start`` at R_on <= ratio R_on,0."""
    th_start, ratio, bv0, r0, maxiter = args
    t0 = time.perf_counter()
    ev = make_ev(0)
    sc = ev.sc
    bounds = Bounds(DESIGN, sc)
    th, res, pr, its = run_slsqp(ev, th_start, bv0, r0, ratio, bounds, maxiter=min(40, maxiter))
    b = ev.objective(th, grad=False)
    return {
        "ratio": ratio, "theta": th.tolist(), "bv_obj": b["BV"], "bv_k": b["bv_k"].tolist(),
        "R": ev.ron(th, grad=False)["R"], "status": int(res.status), "message": res.message,
        "nit": int(res.nit), "bound_violation": bounds.violation(th),
        "failed_evals": ev.stats["failed_evals"], "seconds": time.perf_counter() - t0,
        "failure": pr.failure, "n_evals": pr.n_eval,
    }  # fmt: skip


# ============================================================ derived analyses
def u1_richardson(vals):
    e1, e2 = vals[1] - vals[0], vals[2] - vals[1]
    if e1 == 0.0 or e2 == 0.0 or e1 * e2 < 0.0:
        return {"order": None, "richardson": None, "rel_dist": None, "monotone": False,
                "pass": False}  # fmt: skip
    p = float(np.log2(abs(e1 / e2)))
    rich = vals[2] + e2 / (2.0**p - 1.0)
    rel = abs(vals[2] - rich) / abs(rich)
    return {"order": p, "richardson": float(rich), "rel_dist": float(rel), "monotone": True,
            "pass": bool(rel <= U1_TOL)}  # fmt: skip


def sensitivity_maps(ev: CornerAdjoint, th):
    """gamma(x, y) = dJ/dN_D,i / vol_s,i, sign D = N_A - N_D (gamma_D = -gamma_ND)."""
    b = ev.objective(th)
    r = ev.ron(th)
    vol = ev.dev.reg.vol_s
    ok = vol > 0
    out = {}
    for name, g in (("BV", b["g_nd"]), ("R_on", r["g_nd"] / r["R"])):
        gd = np.zeros_like(g)
        gd[ok] = -g[ok] / vol[ok]
        out[name] = gd
    return out


def plot_per_cut(dense, cuts_um, s13_dense, out_png, man):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=(8.2, 5.0), constrained_layout=True)
    sty = {("initial", 0): ("#0072B2", "-"), ("initial", 1): ("#0072B2", "--"),
           ("final", 0): ("#D55E00", "-"), ("final", 1): ("#D55E00", "--")}  # fmt: skip
    for (tag, lev), (col, ls) in sty.items():
        d = dense.get((tag, lev))
        if d is None:
            continue
        ax.plot(d["x_um"], d["bv"], ls, color=col, lw=1.2,
                label=f"{tag} design, {'h' if lev == 0 else 'h/2'}: min {d['bv_min']:.0f} V "
                      f"@ {d['x_argmin_um']:.2f} um")  # fmt: skip
        if lev == 0:
            ax.plot(d["cut_x_um"], d["bv_at_objective_columns"], "o", color=col, ms=4)
    if s13_dense is not None:
        ax.plot(s13_dense["x_um"], s13_dense["bv"], ":", color="#009E73", lw=1.2,
                label=f"stage-13 optimum on this device, h: min {s13_dense['bv_min']:.0f} V")  # fmt: skip
    for c in cuts_um:
        ax.axvline(c, color="0.85", lw=0.6, zorder=0)
    ax.set_xlabel("x  [um]   (JFET centre 0, JFET/body edge 0.8, body contact 3.0)")
    ax.set_ylabel("avalanche BV of the vertical cut  [V]")
    ax.set_title("per-cut BV (49 dense lines; dots = 10 objective columns; grey = objective cuts)")
    ax.set_ylim(300.0, 900.0)
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)
    fig.text(0.5, 0.005, man, ha="center", fontsize=7, color="0.4")
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


def plot_profile(ev: CornerAdjoint, sc, designs, out_png, man):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    dev, prm = ev.dev, ev.prm
    ys = np.linspace(-(prm.d_body + prm.w_drift), -prm.d_body, 200)
    um = sc.l_d * 1e4
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6), constrained_layout=True)
    cols = ["#0072B2", "#D55E00", "#E69F00", "#CC79A7", "#009E73", "0.4"]
    for (name, th), col in zip(designs.items(), cols, strict=False):
        nd = dev.drift_profile(th, ys) * sc.n0
        ax[0].plot(nd, -ys * um, color=col, label=f"{name} (JFET {dev.jfet_level(th) * sc.n0:.2e})")
        ev.reset_caches()
        o = ev.objective(th, grad=False)
        k = o["limiting_cut"]
        X = ev.caches[k][1]
        nodes = ev.cuts[k]
        depth = -ev.mesh.node_xy[nodes, 1]
        e = np.abs(np.diff(X[0::3][nodes]) / np.diff(depth)) * sc.e_unit * 1e-6
        ax[1].plot(e, 0.5 * (depth[1:] + depth[:-1]) * um, color=col,
                   label=f"{name}: limiting cut x={ev.cut_x_um[k]:.2f} um, BV_k={o['bv_k'][k]:.0f} V")  # fmt: skip
    ax[0].set_xscale("log")
    ax[0].set_ylim(ax[0].get_ylim()[::-1])
    ax[0].set_xlabel("drift donor density [cm$^{-3}$]")
    ax[0].set_ylabel("depth below the surface [um]")
    ax[0].set_title("drift profile (body junction at top)")
    ax[0].legend(fontsize=7)
    ax[1].set_ylim(ax[1].get_ylim()[::-1])
    ax[1].set_xlabel("|E| on the limiting cut at its BV  [MV/cm]")
    ax[1].set_ylabel("depth [um]")
    ax[1].set_title("field on each design's limiting objective cut at BV_k")
    ax[1].legend(fontsize=7)
    fig.text(0.5, 0.005, man, ha="center", fontsize=7, color="0.4")
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


# ============================================================ main
def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--maxiter", type=int, default=MAX_ITER)
    ap.add_argument("--workers", type=int, default=min(10, os.cpu_count() or 1))
    args = ap.parse_args(argv)
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    man = build_manifest(REPO, sys.argv, list(SEEDS))
    (out / "manifest.json").write_text(json.dumps(man, indent=2) + "\n")
    t_start = time.perf_counter()
    reg = TracedRegistry.from_json()
    sc = SiC4H(reg)
    ev = CornerAdjoint(sc, DESIGN, 0, r_um=R_UM, cuts_um=CUTS_UM, tau_v=TAU_V)
    bounds = Bounds(DESIGN, sc)
    T: dict = {}
    res: dict = {}
    ctx = mp.get_context("fork")

    def dump(name="results.json"):
        (out / name).write_text(json.dumps(res, indent=2, default=float) + "\n")

    def timed(name, fn, *a, **k):
        t0 = time.perf_counter()
        r = fn(*a, **k)
        T[name] = time.perf_counter() - t0
        print(f"{name}: {T[name]:.1f}s", flush=True)
        return r

    res["config"] = {
        "design": {k: getattr(DESIGN, k) for k in DESIGN.__dataclass_fields__},
        "rounding_length_um": R_UM, "rounding_definition": "F(s) = erfc(-s/(sqrt2 r))/2 (sigma = r)",
        "objective_cuts_um_requested": list(CUTS_UM), "objective_cuts_um_snapped": ev.cut_x_um,
        "K": len(CUTS_UM), "tau_V": TAU_V, "n_dense_cuts": N_DENSE,
        "u1_levels": list(U1_LEVELS), "u1_tol": U1_TOL,
        "mesh_level0": [int(ev.mesh.nx), int(ev.mesh.ny)], "n_nodes_level0": int(ev.mesh.n_nodes),
        "drift_bounds_cm3": S13.ND_DRIFT_BOUNDS, "jfet_bounds_cm3": S13.ND_JFET_BOUNDS,
        "pareto_ratios": [1.0, *PARETO_RATIOS], "eps_fd": EPS_FD, "seeds": list(SEEDS),
        "max_iter": args.maxiter, "V_GS_V": ev.V_GS, "V_DS_V": ev.V_DS,
        "coefficients": "Niwa 2015 (Okuto-Crowell, electron-initiated I_n)",
        "doping": "compensated nd/na node fields (CV averages of erfc-product cell values); "
                  "JFET donors = n_jfet e^theta_J F(x_jfet - x) (complement of the body edge)",
        "offstate_node_classification": "theta = 0 device; n-type nodes in the source box enlarged by 2.5 r tied to the source",
        "x_half_cell_um": ev.x_cm * 1e4, "workers": args.workers,
    }  # fmt: skip

    # ---- initial design (phase 1, serial) --------------------------------------------------
    o0 = timed("initial_objective", lambda: ev.objective(TH0, grad=False))
    r0 = timed("initial_ron", lambda: ev.ron(TH0, grad=False)["R"])
    bv0 = o0["BV"]
    res["initial"] = {"theta": TH0.tolist(), "bv_obj": bv0, "bv_k": o0["bv_k"].tolist(),
                      "weights": o0["w"].tolist(), "R": r0, "limiting_cut": o0["limiting_cut"],
                      "limiting_cut_x_um": ev.cut_x_um[o0["limiting_cut"]]}  # fmt: skip
    gc0 = timed("grad_check_initial", grad_check, ev, TH0)
    res["grad_check_initial"] = gc0
    gam = timed("sensitivity", sensitivity_maps, ev, TH0)
    dump()

    th1, sl, pr, its = timed("slsqp_main", run_slsqp, ev, TH0, bv0, r0, 1.0, bounds, args.maxiter)
    o1, rr1 = ev.objective(th1, grad=False), ev.ron(th1, grad=False)
    bv1, r1 = o1["BV"], rr1["R"]
    print(
        f"main: BV_obj {bv0:.3f} -> {bv1:.3f} V, R {r0:.5e} -> {r1:.5e}, {sl.message}", flush=True
    )
    res["final"] = {
        "theta": th1.tolist(), "bv_obj": bv1, "bv_k": o1["bv_k"].tolist(),
        "weights": o1["w"].tolist(), "R": r1, "R_over_R0_minus_1": r1 / r0 - 1.0,
        "limiting_cut": o1["limiting_cut"], "limiting_cut_x_um": ev.cut_x_um[o1["limiting_cut"]],
        "bound_violation": bounds.violation(th1), "max_violation_during_run": pr.max_violation,
        "flips": o1["flips"],
    }  # fmt: skip
    res["slsqp_main"] = {"status": int(sl.status), "message": sl.message, "nit": int(sl.nit),
                         "nfev": int(sl.nfev), "njev": int(sl.njev), "iterates": its,
                         "evaluations": pr.history}  # fmt: skip
    gc1 = timed("grad_check_final", grad_check, ev, th1)
    res["grad_check_final"] = gc1
    dump()

    # stage-13 optimum (for the comparison), if its results exist
    s13_theta, s13_summary = None, None
    p13 = REPO / "results" / "stage13" / "results.json"
    if p13.exists():
        r13 = json.loads(p13.read_text())
        s13_theta = np.asarray(r13["final"]["theta"], float)
        s13_summary = {"bv_planar_initial": r13["initial"]["bv"], "bv_planar_final": r13["final"]["bv"],
                       "R_initial": r13["initial"]["R"], "R_final": r13["final"]["R"],
                       "corner_diag": r13.get("corner_diag")}  # fmt: skip

    # ---- phase 2 (parallel): Pareto SLSQP runs + dense verification of initial/final ------------
    jobs_dense = [(0, TH0, "initial", True), (0, th1, "final", True),
                  (1, TH0, "initial", True), (1, th1, "final", True)]  # fmt: skip
    if s13_theta is not None:
        jobs_dense.append((0, s13_theta, "stage13_final", True))
    t_ph2 = time.perf_counter()
    pareto_res: dict = {}
    dense: dict = {}
    with ProcessPoolExecutor(args.workers, mp_context=ctx) as pool:
        futs = {}
        for rt in PARETO_RATIOS:
            futs[pool.submit(job_pareto, (th1, rt, bv0, r0, args.maxiter))] = ("pareto", rt)
        for a in jobs_dense:
            futs[pool.submit(job_dense, a)] = ("dense", (a[2], a[0]))
        for f in as_completed(futs):
            kind, key = futs[f]
            r = f.result()
            if kind == "pareto":
                pareto_res[key] = r
            else:
                dense[key] = r
            print(f"phase2 done: {kind} {key} {r['seconds']:.0f}s", flush=True)
    T["phase2_parallel"] = time.perf_counter() - t_ph2
    res["dense"] = {f"{k[0]}_L{k[1]}": v for k, v in dense.items()}
    dump()

    # ---- phase 3 (parallel): U1/U6 peak-field solves, Pareto designs' dense BV --------------
    v_ref = {lev: dense[("initial", lev)]["bv_min"] for lev in (0, 1)}
    jobs_peak = [(lev, TH0, "initial_u1", v_ref[0]) for lev in U1_LEVELS]
    for lev in (0, 1):
        jobs_peak += [(lev, TH0, "initial_u6", v_ref[lev]), (lev, th1, "final_u6", v_ref[lev])]
    jobs_d3 = [
        (0, np.asarray(pareto_res[rt]["theta"]), f"ratio_{rt:g}", True) for rt in PARETO_RATIOS
    ]
    t_ph3 = time.perf_counter()
    peaks: dict = {}
    with ProcessPoolExecutor(args.workers, mp_context=ctx) as pool:
        futs = {pool.submit(job_peak, a): ("peak", (a[2], a[0])) for a in jobs_peak}
        futs.update({pool.submit(job_dense, a): ("dense", (a[2], a[0])) for a in jobs_d3})
        for f in as_completed(futs):
            kind, key = futs[f]
            r = f.result()
            (peaks if kind == "peak" else dense)[key] = r
            print(f"phase3 done: {kind} {key} {r['seconds']:.0f}s", flush=True)
    T["phase3_parallel"] = time.perf_counter() - t_ph3
    res["dense"] = {f"{k[0]}_L{k[1]}": v for k, v in dense.items()}
    res["peaks"] = {f"{k[0]}_L{k[1]}": v for k, v in peaks.items()}
    dump()

    # ---- R1 Pareto set -------------------------------------------------------------------
    ideal = unipolar_limit(sc)
    par = [{"ratio": 1.0, "theta": th1.tolist(), "bv_obj": bv1, "R": r1, "status": int(sl.status),
            "nit": int(sl.nit), "bv_true_min": dense[("final", 0)]["bv_min"],
            "x_argmin_um": dense[("final", 0)]["x_argmin_um"]}]  # fmt: skip
    for rt in PARETO_RATIOS:
        q = dict(pareto_res[rt])
        d = dense[(f"ratio_{rt:g}", 0)]
        q["bv_true_min"], q["x_argmin_um"] = d["bv_min"], d["x_argmin_um"]
        par.append(q)
    par.append({"ratio": None, "theta": TH0.tolist(), "bv_obj": bv0, "R": r0, "note": "initial",
                "bv_true_min": dense[("initial", 0)]["bv_min"],
                "x_argmin_um": dense[("initial", 0)]["x_argmin_um"]})  # fmt: skip
    for q in par:
        for key in ("bv_true_min", "bv_obj"):
            q[f"R_ideal_at_{key}"] = ideal_r_at(ideal, q[key])
            q[f"R_over_ideal_{key}"] = q["R"] / q[f"R_ideal_at_{key}"]
    res["R1_pareto"] = par
    res["R1_unipolar_limit"] = ideal
    res["R1_formula"] = (
        "R_ideal = W_dep/(q mu_n,par(N_D) N_D), W_dep = 2 BV/E_m; triangular_breakdown (I_n = 1, Niwa); stage-13 R1"
    )
    dump()

    # ---- criteria ------------------------------------------------------------------------
    d_ = gc0["directions"]
    u2 = bool(all(r["bv"]["rel_err"] <= U2_TOL_BV for r in d_)
              and all(r["ron"]["rel_err"] <= U2_TOL_R for r in d_))  # fmt: skip
    d1 = gc1["directions"]
    tol_bv = 1e-6 * max(1.0, gc1["grad_bv_norm"])
    tol_r = 1e-6 * max(1.0, gc1["grad_ron_over_R0_norm"])
    tol_r_lit = 1e-6 * max(1.0, gc1["grad_ron_norm"])
    u3_bv = all(r["bv"]["abs_diff"] <= tol_bv for r in d1)
    u3_r = all(r["ron_norm"]["abs_diff"] <= tol_r for r in d1)
    u3_r_lit = all(r["ron"]["abs_diff"] <= tol_r_lit for r in d1)
    u3 = bool(u3_bv and u3_r)
    viol = bounds.violation(th1)
    u4 = bool(r1 <= r0 * (1.0 + FEAS_TOL) and viol <= 1e-9)

    ini = {lev: dense[("initial", lev)] for lev in (0, 1)}
    fin = {lev: dense[("final", lev)] for lev in (0, 1)}
    imp = {lev: fin[lev]["bv_min"] - ini[lev]["bv_min"] for lev in (0, 1)}
    u5 = bool(imp[0] > 0 and imp[1] > 0 and imp[1] >= 0.5 * imp[0])

    u1_vals = [peaks[("initial_u1", lev)]["E_peak_MV_cm"] for lev in U1_LEVELS]
    u1r = u1_richardson(u1_vals)
    u1 = bool(u1r["pass"])

    u6_rec = {}
    for lev in (0, 1):
        pi, pf = peaks[("initial_u6", lev)], peaks[("final_u6", lev)]
        u6_rec[str(lev)] = {"v_ref": v_ref[lev], "E_initial": pi["E_peak_MV_cm"],
                            "E_final": pf["E_peak_MV_cm"],
                            "loc_initial": [pi["x_um"], pi["y_um"]], "loc_final": [pf["x_um"], pf["y_um"]],
                            "E_srcbox_initial": pi["E_source_region_peak_MV_cm"],
                            "E_srcbox_final": pf["E_source_region_peak_MV_cm"]}  # fmt: skip
    u6 = bool(all(v["E_final"] <= v["E_initial"] for v in u6_rec.values()))

    n_failed = ev.stats["failed_evals"] + sum(d["stats"]["failed_evals"] for d in dense.values())
    pareto_failures = {
        str(k): v["failure"] for k, v in pareto_res.items() if v["failure"]
    }  # reported, not gating
    u7 = bool(sl.status == 0 and sl.nit <= args.maxiter and n_failed == 0
              and all(d["resolved_all"] for d in dense.values()))  # fmt: skip
    crit = {"U1": u1, "U2": u2, "U3": u3, "U4": u4, "U5": u5, "U6": u6, "U7": u7}
    crit = {k: bool(v) for k, v in crit.items()}

    # ---- blind-spot analysis (report) ------------------------------------------------------
    cut_x = np.asarray(CUTS_UM)
    blind = {}
    for tag in ("initial", "final"):
        for lev in (0, 1):
            d = dense[(tag, lev)]
            x_am = d["x_argmin_um"]
            bvc = np.asarray(d["bv_at_objective_columns"])
            blind[f"{tag}_L{lev}"] = {
                "dense_min_V": d["bv_min"], "x_argmin_um": x_am,
                "nearest_objective_cut_um": float(cut_x[np.argmin(np.abs(cut_x - x_am))]),
                "distance_to_nearest_cut_um": float(np.min(np.abs(cut_x - x_am))),
                "min_over_objective_columns_V": float(np.nanmin(bvc)),
                "limiting_objective_column_um": float(cut_x[int(np.nanargmin(bvc))]),
                "gap_columns_minus_dense_V": float(np.nanmin(bvc) - d["bv_min"]),
                "n_unbracketed": d["n_unbracketed"],
            }  # fmt: skip
    for lev in (0, 1):
        fi, ii = fin[lev], ini[lev]
        bi, bf = np.asarray(ii["bv"]), np.asarray(fi["bv"])
        both = ~np.isnan(bi) & ~np.isnan(bf)
        blind[f"L{lev}_x_where_final_lower_than_initial_um"] = [
            float(x) for x, a, b in zip(fi["x_um"], bi, bf, strict=True) if not (np.isnan(a) or np.isnan(b)) and b < a
        ]  # fmt: skip
        blind[f"L{lev}_worst_cut_BV_drop_vs_initial_V"] = (
            float(np.max(bi[both] - bf[both])) if both.any() else None
        )
    blind["source_junction_box_peak_E_MV_cm"] = {k: [v["E_srcbox_initial"], v["E_srcbox_final"]]
                                                  for k, v in u6_rec.items()}  # fmt: skip
    blind["global_peak_location_um"] = {
        k: [v["loc_initial"], v["loc_final"]] for k, v in u6_rec.items()
    }
    res["blind_spot_analysis"] = blind

    # ---- R2: figures, stage-13 comparison -------------------------------------------------
    figdir = (REPO / "figures") if out.resolve() == DEFAULT_OUT.resolve() else out / "figures"
    label = f"stage 13b | commit {str(man.get('git_commit'))[:8]} | {man['timestamp_utc'][:19]}Z"
    fig_err = None
    try:
        S13.plot_maps(ev, gam, sc, figdir / "stage13b_sensitivity_maps.png", label)
        plot_per_cut({(k[0], k[1]): v for k, v in dense.items()}, CUTS_UM,
                     dense.get(("stage13_final", 0)), figdir / "stage13b_per_cut_bv.png", label)  # fmt: skip
        dsg = {"initial": TH0, "optimized (R<=R0)": th1}
        for q in par[1:-1]:
            dsg[f"R<={q['ratio']:g}R0"] = np.asarray(q["theta"])
        if s13_theta is not None:
            dsg["stage-13 optimum"] = s13_theta
        plot_profile(ev, sc, dsg, figdir / "stage13b_optimized_profile.png", label)
    except Exception as exc:  # figures are optional
        fig_err = f"{type(exc).__name__}: {exc}"
    res["figures"] = {"dir": str(figdir), "error": fig_err}
    ys = np.linspace(-(ev.prm.d_body + ev.prm.w_drift), -ev.prm.d_body, 7)
    res["R2"] = {
        "node_x_um": (ev.mesh.node_xy[:, 0] * sc.l_d * 1e4).tolist(),
        "node_y_um": (ev.mesh.node_xy[:, 1] * sc.l_d * 1e4).tolist(),
        "gamma_BV": gam["BV"].tolist(), "gamma_Ron": gam["R_on"].tolist(),
        "profile_depth_below_body_um": [(-y - ev.prm.d_body) * sc.l_d * 1e4 for y in ys],
        "N_D_initial_cm3": (ev.dev.drift_profile(TH0, ys) * sc.n0).tolist(),
        "N_D_final_cm3": (ev.dev.drift_profile(th1, ys) * sc.n0).tolist(),
        "jfet_initial_cm3": ev.dev.jfet_level(TH0) * sc.n0, "jfet_final_cm3": ev.dev.jfet_level(th1) * sc.n0,
        "N_D_pareto_cm3": {f"{q['ratio']}": (ev.dev.drift_profile(np.asarray(q["theta"]), ys) * sc.n0).tolist() for q in par if q.get("ratio")},
        "jfet_pareto_cm3": {f"{q['ratio']}": ev.dev.jfet_level(np.asarray(q["theta"])) * sc.n0 for q in par if q.get("ratio")},
        "stage13": s13_summary,
        "stage13_theta": None if s13_theta is None else s13_theta.tolist(),
        "stage13_optimum_on_rounded_device": None if s13_theta is None else {
            "bv_true_min": dense[("stage13_final", 0)]["bv_min"], "x_argmin_um": dense[("stage13_final", 0)]["x_argmin_um"],
            "R": dense[("stage13_final", 0)]["R"], "R_over_R0": dense[("stage13_final", 0)]["R"] / r0,
            "bv_at_objective_columns": dense[("stage13_final", 0)]["bv_at_objective_columns"]},
        "flags_stage12": FLAGS12,
        "note": "qualitative comparison only; sign convention D = N_A - N_D",
    }  # fmt: skip

    res["runtime_seconds"] = time.perf_counter() - t_start
    res["peak_rss_MiB"] = rss_mib()
    res["solver_stats_main"] = dict(ev.stats)
    summary = {
        "criteria": crit,
        "gate_passed": all(crit.values()),
        "values": {
            "bv_obj_initial_V": bv0, "bv_obj_final_V": bv1, "R_initial": r0, "R_final": r1,
            "objective_cuts_um": ev.cut_x_um, "tau_V": TAU_V, "r_um": R_UM,
            "bv_k_initial": o0["bv_k"].tolist(), "bv_k_final": o1["bv_k"].tolist(),
            "U1": {"E_peak_MV_cm_h_h2_h4": u1_vals, "v_volt": v_ref[0], **u1r,
                   "locations_um": [[peaks[("initial_u1", lev)]["x_um"], peaks[("initial_u1", lev)]["y_um"]] for lev in U1_LEVELS]},
            "U2": {"worst_rel_BV_obj": max(r["bv"]["rel_err"] for r in d_),
                   "worst_rel_Ron": max(r["ron"]["rel_err"] for r in d_),
                   "fd_spread_worst": [max(r["bv"]["fd_spread"] / abs(r["bv"]["richardson"]) for r in d_),
                                       max(r["ron"]["fd_spread"] / abs(r["ron"]["richardson"]) for r in d_)]},
            "U3": {"worst_abs_BV_V": max(r["bv"]["abs_diff"] for r in d1), "tol_BV": tol_bv,
                   "worst_abs_Ron_over_R0": max(r["ron_norm"]["abs_diff"] for r in d1), "tol_Ron_over_R0": tol_r,
                   "worst_abs_Ron_ohm_cm2_literal": max(r["ron"]["abs_diff"] for r in d1),
                   "tol_Ron_literal": tol_r_lit, "Ron_literal_pass": bool(u3_r_lit),
                   "grad_norms": [gc1["grad_bv_norm"], gc1["grad_ron_norm"], gc1["grad_ron_over_R0_norm"]]},
            "U4": {"R_over_R0_minus_1": r1 / r0 - 1.0, "bound_violation": viol},
            "U5": {"min_bv_initial_h_h2": [ini[0]["bv_min"], ini[1]["bv_min"]],
                   "min_bv_final_h_h2": [fin[0]["bv_min"], fin[1]["bv_min"]],
                   "improvement_h": imp[0], "improvement_h2": imp[1],
                   "ratio": imp[1] / imp[0] if imp[0] != 0 else None,
                   "x_argmin_initial_um": [ini[0]["x_argmin_um"], ini[1]["x_argmin_um"]],
                   "x_argmin_final_um": [fin[0]["x_argmin_um"], fin[1]["x_argmin_um"]],
                   "R_h": [dense[("initial", 0)]["R"], dense[("final", 0)]["R"]],
                   "R_h2": [dense[("initial", 1)]["R"], dense[("final", 1)]["R"]]},
            "U6": u6_rec,
            "U7": {"status": int(sl.status), "nit": int(sl.nit), "message": sl.message,
                   "failed_evals": n_failed, "pareto_failures": pareto_failures, "pareto_status": {str(k): v["status"] for k, v in pareto_res.items()},
                   "dense_resolved_all": {f"{k[0]}_L{k[1]}": v["resolved_all"] for k, v in dense.items()}},
            "R1": [{k: q.get(k) for k in ("ratio", "bv_obj", "bv_true_min", "x_argmin_um", "R", "R_over_ideal_bv_true_min", "R_over_ideal_bv_obj")} for q in par],
            "R1_formula_check_max_rel": max(r["rel_formula_diff"] for r in ideal),
            "blind_spot": blind,
            "stage13_comparison": res["R2"]["stage13_optimum_on_rounded_device"],
            "flips_max": ev.stats["flips_max"],
        },
        "seconds": T,
        "runtime_seconds": res["runtime_seconds"],
        "peak_rss_mib": res["peak_rss_MiB"],
    }  # fmt: skip
    dump()
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    print(json.dumps(summary["criteria"], indent=2))
    print(f"runtime {res['runtime_seconds']:.1f} s")


if __name__ == "__main__":
    main()
