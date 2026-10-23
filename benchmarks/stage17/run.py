"""Stage-17 benchmark (docs/benchmarks/stage17_channel_aware_opt.md): J1-J7 (gating), R1-R3.

The stage-13d problem (5 design variables, soft-min BV on h, E_ox constraint at 400 V on h) with
R_on from the stage-16 channel-mobility model on h/2 (mesh level 1), using the exact adjoint of
the Picard fixed point (widegap.adjoint_channel).  The off-state (BV, E_ox) is the 13d evaluator
on level 0; the on-state (R_on,ch, its adjoint, the R_on constraint) is a separate level-1 solve.

    maximize BV_obj(theta)  s.t.  R_on,ch(theta) <= R_on,ch(0) (level 1),
                                  E_ox(theta) <= E_ox(0) (400 V, level 0),  bounds,

SLSQP (13d settings), Pareto at R/R0 = 1.05, 1.1, 1.2, J5 dense true-minimum BV (49 lines) on h and
h/2, J6 independent R_on,ch re-check on h/4 (level 2), R3 field-line BV (stage-15b method).

Interpretations (the spec leaves them open; recorded in results.json "interpretations"):
* J1 gates the R_on,ch gradient on the five unit directions (vector norm); the soft-min BV and E_ox
  gradients are re-checked exactly as 13d D2 (seeded directions, relative tolerances 1e-4).
* J2 uses the 13d D3 form for all three gradients in units of the initial values (BV/BV_0, R/R_0,
  E_ox/E_ox,0): worst absolute directional-derivative error <= 1e-6 max(1, ||grad||).  BV and E_ox
  use the seeded directions of 13d, R_on,ch the five unit directions.
* J4 is evaluated at EVERY channel on-state evaluation of every job (FD points, SLSQP, Pareto,
  h/4, 13d final), plus the surface-sign check.
* J7: "every Picard and Newton solve converged" is gated on evaluations: no failed evaluation (any
  job, on- or off-state), no non-converged Picard loop.  Every non-finite SuperLU call (counted
  over all jobs, incl. the 13b/15b helper jobs) must lie inside a failed (rejected) Newton solve:
  the cold-march step halving of the on-state / off-state produces a few such rejected trial
  steps per job; they are REPORTED (failed_newton_calls, nonfinite_calls) and any non-finite call
  elsewhere (converged Newton, adjoint or tangent solves) fails J7.  The Pareto SLSQP exit codes
  are reported, not gated (as in 13d).
  Fd 2 is captured to <out>/stderr.log and scanned for "On entry to".

Writes <out>/{results,summary,manifest}.json and stderr.log (default results/stage17) and figures
stage17_*.png (repo figures/ for the default output, else <out>/figures).  ``--dev`` marks the run
as DEV (refused when --out lies under results/).
"""

# ruff: noqa: E402
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import functools
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

import widegap.adjoint2d as _adj2d
import widegap.linalg as _linalg
from widegap.adjoint2d import solve_transposed
from widegap.adjoint_channel import ChannelRon, mobility_part_theta
from widegap.adjoint_sic import richardson_fd
from widegap.channel_mobility import NOGUCHI_NAMES, NoguchiParams
from widegap.fv2d.mos import MOSDriftDiffusion
from widegap.manifest import build_manifest
from widegap.materials.registry import ProvenanceError, check_entry
from widegap.sic import REGISTRY_PATH

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[2]
DEFAULT_OUT = REPO / "results" / "stage17"
R13D = REPO / "results" / "stage13d"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


# stage 15b first: it re-registers "stage13b_run"; the 13d chain below registers the final one
S15B = _load("stage15b_run", REPO / "benchmarks/stage15b/run.py")
M13D = _load("stage13d_run", REPO / "benchmarks/stage13d/run.py")  # fixed names: fork pickling
M = M13D.M
S13 = M.S13
DESIGN, TH0, Bounds = M13D.DESIGN, M13D.TH0, M13D.Bounds
MEASURES, V_BLK_V = M13D.MEASURES, M13D.V_BLK_V

# ---- study parameters (spec) ----------------------------------------------------------------
LEVEL_ON = 1  # h/2: on-state
LEVEL_OFF = 0  # h: off-state (as 13d)
LEVEL_RECHECK = 2  # h/4 (J6)
J1_TOL = 1e-5
J2_TOL = 1e-6
D2_TOL_REL = 1e-4  # 13d D2 for BV_obj and E_ox
J5_RATIO = 0.5
J6_TOL = 1e-3
FEAS_TOL = M.FEAS_TOL
PARETO_RATIOS = M.PARETO_RATIOS
EPS_FD = M.EPS_FD
N_S_FL, DS_FL_UM = S15B.N_S, S15B.DS_UM  # 400 seeds, 0.00125 um (15b)
Y_BANDS_UM = (0.0, -0.02, -0.05, -0.1, -0.3, -0.8, -1e9)

# ---- instrumentation (J7) -------------------------------------------------------------------
NF = {"adjoint": 0, "converged_newton": 0, "failed_newton": 0, "failed_newton_calls": 0}
_orig_adj_solve = _adj2d.safe_spsolve


def _traced_adjoint_solve(a, b):
    n0 = _linalg.nonfinite_calls
    x = _orig_adj_solve(a, b)
    NF["adjoint"] += _linalg.nonfinite_calls - n0
    return x


_adj2d.safe_spsolve = _traced_adjoint_solve


def _traced_solve_bc(self, X0, bc, *a, **k):
    """Class-level wrapper of ``MOSDriftDiffusion.solve_bc``: counts non-finite SuperLU inputs per
    Newton solve, split by converged / failed.  Applied before any evaluator exists, so every model
    built by the 13b / 15b / 13d helpers is covered."""
    n0 = _linalg.nonfinite_calls
    sol = _orig_solve_bc(self, X0, bc, *a, **k)
    d = _linalg.nonfinite_calls - n0
    if sol.report.converged:
        NF["converged_newton"] += d
    else:
        NF["failed_newton"] += d
        NF["failed_newton_calls"] += 1
    return sol


_orig_solve_bc = MOSDriftDiffusion.solve_bc
MOSDriftDiffusion.solve_bc = _traced_solve_bc


_BASE = {"NF": dict.fromkeys(NF, 0), "nf": 0}


def tracked(fn):
    """Job wrapper: the J7 counters of a job are the deltas over the job (workers are reused)."""

    @functools.wraps(fn)
    def wrapper(*a, **k):
        _BASE["NF"], _BASE["nf"] = dict(NF), int(_linalg.nonfinite_calls)
        return fn(*a, **k)

    return wrapper


def rss_mib():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def _lst(a):
    return np.asarray(a, dtype=float).tolist()


def make_ev_off():
    return M13D.make_ev13d(LEVEL_OFF)


def make_cr(level, params=None):
    ev = M13D.make_ev13d(level, with_cuts=False)
    p = params or NoguchiParams.from_registry(ev.sc.reg)
    return ChannelRon(ev, p)


def guard_summary(guards):
    if not guards:
        return {"n": 0, "n_bad": 0, "all_ok": True}
    bad = [g for g in guards if not g["ok"]]
    return {
        "n": len(guards), "n_bad": len(bad), "all_ok": not bad,
        "bad_examples": [{k: g[k] for k in ("theta", "sign_changes", "gated_near_ns_clamp",
                                            "gated_near_eeff_clamp")} for g in bad[:3]],
        "min_gated_ns_cm2": min(g["gated_ns_min"] for g in guards),
        "min_gated_eeff_raw_MV_cm": min(g["gated_eeff_raw_min"] for g in guards),
        "max_sign_changes": max(g["sign_changes"] for g in guards),
        "max_ungated_near_ns_clamp": max(g["ungated_near_ns_clamp"] for g in guards),
        "max_surface_NA_rel_change": max(g["na_surface_max_rel_change"] for g in guards),
    }  # fmt: skip


def job_diag(crs=(), evs=()):
    """J7/J4 bookkeeping of one job (never hides anything)."""
    stats = {"picard_evals": 0, "picard_iterations": 0, "failed_picard": 0, "failed_evals": 0,
             "max_picard_iterations": 0, "failed_attempts_on": 0, "attempts_on": 0}  # fmt: skip
    guards = []
    for c in crs:
        for k in ("picard_evals", "picard_iterations", "failed_picard"):
            stats[k] += c.stats[{"picard_evals": "evals", "picard_iterations": "picard_iterations",
                                 "failed_picard": "failed_picard"}[k]]  # fmt: skip
        stats["failed_evals"] += c.stats["failed_evals"]
        if c.picard_log:
            stats["max_picard_iterations"] = max(
                stats["max_picard_iterations"], max(p["iterations"] for p in c.picard_log)
            )
        stats["failed_attempts_on"] += c.ev.stats["failed_attempts"]
        stats["attempts_on"] += c.ev.stats["attempts"]
        guards += c.guards
    off = {"failed_evals": 0, "failed_attempts": 0, "attempts": 0, "warm_fallbacks": 0}
    for e in evs:
        off["failed_evals"] += e.stats["failed_evals"]
        off["failed_attempts"] += e.stats["failed_attempts"]
        off["attempts"] += e.stats["attempts"]
        off["warm_fallbacks"] += e.stats["warm_fallbacks"]
    return {"on_state": stats, "off_state": off, "guards": guard_summary(guards),
            "nonfinite_calls": int(_linalg.nonfinite_calls) - _BASE["nf"],
            "NF": {k: NF[k] - _BASE["NF"][k] for k in NF}}  # fmt: skip


# ============================================================ composite evaluator for SLSQP
class Eval17:
    """What the 13d ``Problem`` needs: BV_obj and E_ox from the level-0 off-state evaluator, R_on
    from the level-1 channel model (exact adjoint)."""

    def __init__(self, ev_off, cr):
        self.ev_off, self.cr, self.sc = ev_off, cr, ev_off.sc

    def objective(self, th, grad=True):
        return self.ev_off.objective(th, grad=grad)

    def fields(self, th, grad=True, which=MEASURES, v_volt=None):
        return self.ev_off.fields(th, grad=grad, which=which, v_volt=v_volt)

    def ron(self, th, grad=True):
        o = self.cr.ron(th, grad=grad)
        out = {"R": o["R"]}
        if grad:
            out["grad"] = o["grad"]
        return out


# ============================================================ gradient checks
def fd_ron_unit(cr, th, k, level_tag=""):
    u = np.zeros(5)
    u[k] = 1.0
    rich, f1, f2 = richardson_fd(lambda t: cr.ron(t, grad=False)["R"], th, u, EPS_FD)
    return {"k": k, "richardson": float(rich), "fd_eps": float(f1), "fd_eps_half": float(f2),
            "spread": float(abs(f1 - f2))}  # fmt: skip


@tracked
def job_fd_ron(args):
    """Richardson FD of R_on,ch along one unit direction (4 channel solves) on ``level``."""
    th, k, level, tag = args
    t0 = time.perf_counter()
    cr = make_cr(level)
    out = {"tag": tag, "level": level, **fd_ron_unit(cr, np.asarray(th, float), k)}
    out["diag"] = job_diag([cr])
    out["seconds"] = time.perf_counter() - t0
    return out


@tracked
def job_grad_off(args):
    """13d gradient checks of BV_obj (D2/D3 seeded directions) and the field measures at theta."""
    th, e0, tag = args
    t0 = time.perf_counter()
    ev = make_ev_off()
    th = np.asarray(th, float)
    gc = M.grad_check(ev, th)
    gf = M13D.grad_check_fields(ev, th, e0)
    return {"tag": tag, "grad_check": gc, "grad_check_fields": gf, "diag": job_diag([], [ev]),
            "seconds": time.perf_counter() - t0}  # fmt: skip


@tracked
def job_ron_point(args):
    """R_on,ch (and the exact / lagged gradient) at one theta on ``level`` (level-2: no gradient)."""
    th, level, tag, grad = args
    t0 = time.perf_counter()
    cr = make_cr(level)
    th = np.asarray(th, float)
    o = cr.ron(th, grad=grad, detail=grad)
    out = {"tag": tag, "level": level, "theta": th.tolist(), "R": o["R"],
           "n_nodes": int(cr.ev.mesh.n_nodes), "picard_iterations": o["iterations"],
           "guard": o["guard"]}  # fmt: skip
    if grad:
        out["grad"], out["grad_lagged"] = _lst(o["grad"]), _lst(o["grad_lagged"])
    out["diag"] = job_diag([cr])
    out["seconds"] = time.perf_counter() - t0
    return out


@tracked
def job_r2(args):
    """R2: lagged-gradient error and the JFET-level (theta_5) sensitivity split, channel vs bulk."""
    th, level, tag = args
    t0 = time.perf_counter()
    cr = make_cr(level)
    ev, cmap = cr.ev, cr.cmap
    th = np.asarray(th, float)
    o = cr.ron(th, detail=True)
    g, gl = o["grad"], o["grad_lagged"]
    st = o["state"]
    mu_ch = o["picard"]["mu_edge"].copy()
    # bulk model at the same theta: R, gradient, adjoint weights
    ev.set_theta(th)
    bulk = ev.ron(th)
    s, bc = ev.on_state()
    mod = ev.mod
    rows, _ = mod._all_bc(bc)
    dn = ev.dev.nodes["drain"]
    sel = np.r_[3 * dn + 1, 3 * dn + 2]
    import scipy.sparse as sp

    jb = sp.csr_matrix(mod._jacobian_core(s.X))
    lam = solve_transposed(mod.jacobian_bc(s.X, rows), np.asarray(jb[sel].sum(axis=0)).ravel())
    lam[rows] = 0.0
    w_b = -lam
    w_b[sel] += 1.0
    i_b = mod.contact_current(s, dn)
    fac_b = -bulk["R"] / i_b
    fac_c = -o["R"] / o["I_scaled"]
    xy = ev.mesh.node_xy
    um = float(ev.sc.l_d) * 1e4
    a, b = mod.a, mod.b
    mid_y = 0.5 * (xy[a, 1] + xy[b, 1]) * um
    mid_x = 0.5 * (xy[a, 0] + xy[b, 0]) * um
    k = 4
    res = {}
    for name, w, X, mu, zero, fac in (
        ("channel_exact", o["w"], st["X"], mu_ch, cmap.edges, fac_c),
        ("channel_lagged", o["w_lagged"], st["X"], mu_ch, cmap.edges, fac_c),
        ("bulk", w_b, s.X, cr.cmap.mu_bulk, None, fac_b),
    ):
        mod.mu_edge = np.asarray(mu, float).copy()
        edge, node = mobility_part_theta(ev, zero, w, X, k)
        edge, node = fac * edge, fac * node
        bands_e, bands_n = [], []
        for lo, hi in zip(Y_BANDS_UM[:-1], Y_BANDS_UM[1:], strict=True):
            me = (mid_y <= lo + 1e-12) & (mid_y > hi)
            mn = (xy[:, 1] * um <= lo + 1e-12) & (xy[:, 1] * um > hi)
            bands_e.append(float(edge[me].sum()))
            bands_n.append(float(node[mn].sum()))
        order = np.argsort(-np.abs(edge))[:8]
        res[name] = {
            "total_edge_mobility_part": float(edge.sum()),
            "total_node_ionization_part": float(node.sum()),
            "by_y_band_edge": bands_e, "by_y_band_node": bands_n,
            "x_edges_sum": float(edge[ev.mesh.edge_dir == 0].sum()),
            "y_edges_sum": float(edge[ev.mesh.edge_dir == 1].sum()),
            "top_edges": [{"x_um": float(mid_x[i]), "y_um": float(mid_y[i]),
                           "dir": int(ev.mesh.edge_dir[i]), "contrib": float(edge[i])}
                          for i in order],
        }  # fmt: skip
    mod.mu_edge = cr.cmap.mu_bulk.copy()
    out = {"tag": tag, "level": level, "theta": th.tolist(), "R_channel": o["R"],
           "R_bulk": bulk["R"], "grad_channel_exact": _lst(g), "grad_channel_lagged": _lst(gl),
           "grad_bulk": _lst(bulk["grad"]),
           "lagged_rel_error_norm": float(np.linalg.norm(gl - g) / np.linalg.norm(g)),
           "lagged_rel_error_per_component": _lst(np.abs(gl - g) / np.abs(g)),
           "channel_over_bulk_per_component": _lst(g / bulk["grad"]),
           "y_bands_um": [[lo, hi] for lo, hi in zip(Y_BANDS_UM[:-1], Y_BANDS_UM[1:], strict=True)],
           "theta5": res,
           "theta5_check": {"exact_total": float(res["channel_exact"]["total_edge_mobility_part"]
                                                 + res["channel_exact"]["total_node_ionization_part"]),
                            "grad5": float(g[4])},
           "diag": job_diag([cr], [])}  # fmt: skip
    out["seconds"] = time.perf_counter() - t0
    return out


@tracked
def job_pareto17(args):
    """Independent SLSQP run from ``th_start`` at R_on,ch <= ratio R_on,ch,0 and E_ox <= E_ox,0."""
    th_start, ratio, bv0, r0, e0, maxiter = args
    t0 = time.perf_counter()
    ev_off = make_ev_off()
    cr = make_cr(LEVEL_ON)
    comp = Eval17(ev_off, cr)
    bounds = Bounds(DESIGN, ev_off.sc)
    th, res, pr, its = M13D.run_slsqp(comp, th_start, bv0, r0, e0, ratio, bounds,
                                      maxiter=min(40, maxiter))  # fmt: skip
    b = ev_off.objective(th, grad=False)
    f = ev_off.fields(th, grad=False)
    out = {
        "ratio": ratio, "theta": th.tolist(), "bv_obj": b["BV"], "bv_k": b["bv_k"].tolist(),
        "R": cr.ron(th, grad=False)["R"], "E_ox": f["E_ox"]["S"],
        "fields": {k: f[k]["S"] for k in MEASURES}, "status": int(res.status),
        "message": res.message, "nit": int(res.nit), "bound_violation": bounds.violation(th),
        "failure": pr.failure, "n_evals": pr.n_eval, "diag": job_diag([cr], [ev_off]),
        "seconds": time.perf_counter() - t0,
    }  # fmt: skip
    return out


@tracked
def _with_diag(r):
    """Attach the J7 deltas of the job and its own off-state evaluator statistics (if any)."""
    d = job_diag()
    st = r.get("stats") or {}
    for k in ("failed_evals", "failed_attempts", "attempts", "warm_fallbacks"):
        d["off_state"][k] += int(st.get(k, 0)) if isinstance(st.get(k, 0), (int, float)) else 0
    r["diag"] = d
    return r


def job_dense17(args):
    return _with_diag(M.job_dense(args))


@tracked
def job_fields17(args):
    return _with_diag(M13D.job_fields(args))


@tracked
def job_fl17(args):
    """Field-line BV (stage-15b method) of one design on one mesh."""
    level, theta, tag, n_s, ds_um = args
    r = S15B.job_fl((level, theta, tag, n_s, ds_um, False))
    r["diag"] = job_diag()  # field-line states use the same Newton solver; no EvalError escapes
    return r


@tracked
def job_ron_13d(args):
    """Channel-model R_on of the 13d designs (final + Pareto) on ``level`` (R1)."""
    thetas, level = args
    t0 = time.perf_counter()
    cr = make_cr(level)
    out = {"level": level, "R": {k: cr.ron(np.asarray(t, float), grad=False)["R"]
                                 for k, t in thetas.items()}}  # fmt: skip
    out["diag"] = job_diag([cr])
    out["seconds"] = time.perf_counter() - t0
    return out


# ============================================================ registry check
def registry_check():
    main = {e["name"]: e for e in json.loads(REGISTRY_PATH.read_text())["parameters"]}
    rows = {}
    for n in NOGUCHI_NAMES:
        row = {"in_sic4h_json": n in main}
        try:
            check_entry(main[n])
            row["check_entry"] = "ok"
        except (ProvenanceError, TypeError, KeyError) as exc:
            row["check_entry"] = f"REFUSED: {exc}"
        rows[n] = row
    return {"entries": rows, "all_usable": all(r["check_entry"] == "ok" for r in rows.values())}


# ============================================================ figures
def plot_all(figdir, res, label):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figdir = Path(figdir)
    figdir.mkdir(parents=True, exist_ok=True)
    col = {"initial": "#0072B2", "17": "#D55E00", "13d": "#009E73"}
    foot = lambda fig: fig.text(0.5, 0.005, label, ha="center", fontsize=7, color="0.4")  # noqa: E731

    fig, ax = plt.subplots(figsize=(8.2, 5.0), constrained_layout=True)
    for key, c, name in (("initial_L0", col["initial"], "initial"),
                         ("final_L0", col["17"], "stage 17"),
                         ("d13_L0", col["13d"], "stage-13d optimum")):  # fmt: skip
        d = res["dense"].get(key)
        if d is not None:
            ax.plot(d["x_um"], d["bv"], "-", color=c,
                    label=f"{name}: min {d['bv_min']:.0f} V @ {d['x_argmin_um']:.2f} um")  # fmt: skip
    ax.set_xlabel("x  [um]")
    ax.set_ylabel("avalanche BV of the vertical cut  [V]")
    ax.set_title("per-cut BV on h (49 dense lines)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    foot(fig)
    fig.savefig(figdir / "stage17_per_cut_bv.png", dpi=130)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
    p17 = res["R1"]["pareto_17"]
    ax.plot([q["R_over_R0"] for q in p17], [q["bv_true_min"] for q in p17], "o-", color=col["17"],
            label="stage 17 (R_on,ch on h/2, E_ox <= E_ox,0)")  # fmt: skip
    p13 = res["R1"]["pareto_13d_channel_R"]
    ax.plot([q["R_ch_over_R0"] for q in p13], [q["bv_true_min"] for q in p13], "s--",
            color=col["13d"], label="stage-13d designs re-evaluated with R_on,ch")  # fmt: skip
    ax.set_xlabel("R_on,ch / R_on,ch,0  (h/2)")
    ax.set_ylabel("dense true-minimum BV on h  [V]")
    ax.set_title("Pareto set")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    foot(fig)
    fig.savefig(figdir / "stage17_pareto.png", dpi=130)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 4.0), constrained_layout=True)
    names = [f"theta_{i}" for i in range(4)] + ["theta_J"]
    x = np.arange(5)
    ax.bar(x - 0.2, res["R1"]["theta_13d_final"], 0.4, color=col["13d"], label="stage 13d")
    ax.bar(x + 0.2, res["final"]["theta"], 0.4, color=col["17"], label="stage 17")
    ax.set_xticks(x, names)
    ax.set_ylabel("design variable")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, axis="y")
    ax.set_title("design variables: 13d vs 17")
    foot(fig)
    fig.savefig(figdir / "stage17_theta.png", dpi=130)
    plt.close(fig)

    r2 = res.get("R2", {}).get("initial")
    if r2:
        fig, ax = plt.subplots(figsize=(8.0, 4.2), constrained_layout=True)
        nb = len(r2["y_bands_um"])
        xb = np.arange(nb)
        labels = [f"{b[0]:g}..{max(b[1], -9):g}" for b in r2["y_bands_um"]]
        for j, (k, c) in enumerate((("bulk", "0.5"), ("channel_lagged", "#56B4E9"),
                                    ("channel_exact", col["17"]))):  # fmt: skip
            v = np.array(r2["theta5"][k]["by_y_band_edge"]) + np.array(
                r2["theta5"][k]["by_y_band_node"])  # fmt: skip
            ax.bar(xb + 0.27 * (j - 1), v, 0.26, color=c, label=k)
        ax.set_xticks(xb, labels, fontsize=7)
        ax.set_xlabel("depth band y [um]")
        ax.set_ylabel(
            r"$\lambda$-weighted $\partial F/\partial\theta_J$ contribution to dR/dtheta_J"
        )
        ax.set_title("JFET-level sensitivity by depth (theta = 0, h/2)")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3, axis="y")
        foot(fig)
        fig.savefig(figdir / "stage17_jfet_sensitivity.png", dpi=130)
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
    # capture C-level stderr (OpenBLAS writes to fd 2)
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
    reg_chk = registry_check()
    ev_off = make_ev_off()
    sc = ev_off.sc
    cr = make_cr(LEVEL_ON)
    comp = Eval17(ev_off, cr)
    bounds = Bounds(DESIGN, sc)
    T: dict = {}
    res: dict = {"dev": bool(args.dev)}
    diags: dict = {}
    ctx = mp.get_context("fork")

    def dump(name="results.json"):
        (out / name).write_text(json.dumps(res, indent=2, default=float) + "\n")

    def timed(name, fn, *a, **k):
        t0 = time.perf_counter()
        r = fn(*a, **k)
        T[name] = time.perf_counter() - t0
        print(f"{name}: {T[name]:.1f}s", flush=True)
        return r

    th13 = np.asarray(json.loads((R13D / "results.json").read_text())["final"]["theta"], float)
    r13d = json.loads((R13D / "results.json").read_text())
    th13_pareto = {f"ratio_{q['ratio']:g}": np.asarray(q["theta"], float)
                   for q in r13d["R2"]["pareto"] if q.get("ratio") and q["ratio"] > 1.0}  # fmt: skip
    res["config"] = {
        "levels": {"on_state": LEVEL_ON, "off_state": LEVEL_OFF, "recheck": LEVEL_RECHECK},
        "mesh_nodes": {"off_L0": int(ev_off.mesh.n_nodes), "on_L1": int(cr.ev.mesh.n_nodes)},
        "d_ch_um": 0.020, "picard_tol": 1e-10, "n_channel_columns": int(len(cr.cmap.cols)),
        "n_gated_columns": int(cr.cmap.col_gated.sum()), "n_channel_edges": int(len(cr.cmap.edges)),
        "pareto_ratios": [1.0, *PARETO_RATIOS], "eps_fd": EPS_FD, "max_iter": args.maxiter,
        "workers": args.workers, "V_blk_V": V_BLK_V, "tolerances": {
            "J1": J1_TOL, "J1_bv_E_ox_D2_rel": D2_TOL_REL, "J2": J2_TOL, "J5_ratio": J5_RATIO,
            "J6": J6_TOL, "feas": FEAS_TOL},
        "fieldline": {"n_seeds": N_S_FL, "ds_um": DS_FL_UM},
        "registry": reg_chk,
        "interpretations": __doc__.split("Interpretations", 1)[1].split("Writes", 1)[0],
    }  # fmt: skip

    # ---- initial design (serial) --------------------------------------------------------------
    o0 = timed("initial_objective", lambda: ev_off.objective(TH0, grad=False))
    f0 = timed("initial_fields", lambda: ev_off.fields(TH0, grad=False))
    r0o = timed("initial_ron_channel", lambda: cr.ron(TH0, grad=False))
    bv0, e0, r0 = o0["BV"], f0["E_ox"]["S"], r0o["R"]
    res["initial"] = {
        "theta": TH0.tolist(), "bv_obj": bv0, "bv_k": o0["bv_k"].tolist(), "R_channel_h2": r0,
        "E_ox": e0, "fields": {k: f0[k]["S"] for k in MEASURES},
        "limiting_cut_x_um": ev_off.cut_x_um[o0["limiting_cut"]], "guard": r0o["guard"],
    }  # fmt: skip
    print(f"initial: BV_obj {bv0:.3f} V, R_ch {r0:.6e}, E_ox {e0:.4f}; J4 {r0o['guard']['ok']}")
    dump()

    pool = ProcessPoolExecutor(args.workers, mp_context=ctx)
    fut: dict = {}

    def submit(key, fn, a):
        fut[pool.submit(fn, a)] = key

    # ---- async jobs that only need theta = 0 (run while the main SLSQP proceeds) --------------
    for k in range(5):
        submit(("j1_fd", k), job_fd_ron, (TH0, k, LEVEL_ON, "initial"))
    submit(("gc0", None), job_grad_off, (TH0, e0, "initial"))
    submit(("r2", "initial"), job_r2, (TH0, LEVEL_ON, "initial"))
    submit(("ron_pt", "initial_L1"), job_ron_point, (TH0, LEVEL_ON, "initial", True))
    submit(("ron_pt", "initial_L2"), job_ron_point, (TH0, LEVEL_RECHECK, "initial", False))
    for lev in (0, 1):
        submit(("dense", f"initial_L{lev}"), job_dense17, (lev, TH0, "initial", False))
        submit(("fields", f"initial_L{lev}_blk"), job_fields17,
               (lev, TH0, f"initial_L{lev}_blk", V_BLK_V, lev == 0))  # fmt: skip
    submit(("fields", "initial_L2_blk"), job_fields17, (2, TH0, "initial_L2_blk", V_BLK_V, False))
    for lev in (0, 1):
        submit(("fl", f"initial_L{lev}"), job_fl17, (lev, TH0, "initial", N_S_FL, DS_FL_UM))

    # ---- main SLSQP (serial, here) ------------------------------------------------------------
    th1, sl, pr, its = timed("slsqp_main", M13D.run_slsqp, comp, TH0, bv0, r0, e0, 1.0, bounds,
                             args.maxiter)  # fmt: skip
    o1 = ev_off.objective(th1, grad=False)
    ro1 = cr.ron(th1, grad=False)
    f1 = ev_off.fields(th1, grad=False)
    bv1, r1, e1 = o1["BV"], ro1["R"], f1["E_ox"]["S"]
    print(f"main: BV_obj {bv0:.3f} -> {bv1:.3f} V, R_ch {r0:.5e} -> {r1:.5e}, "
          f"E_ox {e0:.4f} -> {e1:.4f} MV/cm, {sl.message}", flush=True)  # fmt: skip
    res["final"] = {
        "theta": th1.tolist(), "bv_obj": bv1, "bv_k": o1["bv_k"].tolist(), "R_channel_h2": r1,
        "R_over_R0_minus_1": r1 / r0 - 1.0, "E_ox": e1, "E_ox_over_E0_minus_1": e1 / e0 - 1.0,
        "fields": {k: f1[k]["S"] for k in MEASURES},
        "limiting_cut_x_um": ev_off.cut_x_um[o1["limiting_cut"]],
        "bound_violation": bounds.violation(th1), "max_violation_during_run": pr.max_violation,
        "flips": o1["flips"], "guard": ro1["guard"],
    }  # fmt: skip
    res["slsqp_main"] = {"status": int(sl.status), "message": sl.message, "nit": int(sl.nit),
                         "nfev": int(sl.nfev), "njev": int(sl.njev), "iterates": its,
                         "evaluations": pr.history}  # fmt: skip
    diags["main"] = job_diag([cr], [ev_off])
    dump()

    # ---- jobs that need the final theta -------------------------------------------------------
    for k in range(5):
        submit(("j2_fd", k), job_fd_ron, (th1, k, LEVEL_ON, "final"))
    submit(("gc1", None), job_grad_off, (th1, e0, "final"))
    submit(("r2", "final"), job_r2, (th1, LEVEL_ON, "final"))
    submit(("ron_pt", "final_L1"), job_ron_point, (th1, LEVEL_ON, "final", True))
    submit(("ron_pt", "final_L2"), job_ron_point, (th1, LEVEL_RECHECK, "final", False))
    for rt in PARETO_RATIOS:
        submit(("pareto", rt), job_pareto17, (th1, rt, bv0, r0, e0, args.maxiter))
    for lev in (0, 1):
        submit(("dense", f"final_L{lev}"), job_dense17, (lev, th1, "final", False))
        submit(("fields", f"final_L{lev}_blk"), job_fields17,
               (lev, th1, f"final_L{lev}_blk", V_BLK_V, lev == 0))  # fmt: skip
        submit(("fields", f"d13_L{lev}_blk"), job_fields17,
               (lev, th13, f"d13_L{lev}_blk", V_BLK_V, lev == 0))  # fmt: skip
    submit(("fields", "final_L2_blk"), job_fields17, (2, th1, "final_L2_blk", V_BLK_V, False))
    submit(("dense", "d13_L0"), job_dense17, (0, th13, "d13", False))
    submit(("dense", "d13_L1"), job_dense17, (1, th13, "d13", False))
    for lev in (0, 1):
        submit(("fl", f"final_L{lev}"), job_fl17, (lev, th1, "final", N_S_FL, DS_FL_UM))
        submit(("fl", f"d13_L{lev}"), job_fl17, (lev, th13, "d13", N_S_FL, DS_FL_UM))
    r13_thetas = {"d13_final": th13, **th13_pareto}
    submit(("ron13d", 1), job_ron_13d, (r13_thetas, LEVEL_ON))
    submit(("ron13d", 2), job_ron_13d, ({"d13_final": th13}, LEVEL_RECHECK))

    jobs: dict = {"j1_fd": {}, "j2_fd": {}, "pareto": {}, "dense": {}, "fields": {}, "fl": {},
                  "ron_pt": {}, "r2": {}, "ron13d": {}, "gc0": None, "gc1": None}  # fmt: skip
    t_par = time.perf_counter()
    pending = set(fut)
    while pending:
        for f in as_completed(list(pending)):
            pending.discard(f)
            kind, key = fut[f]
            try:
                r = f.result()
            except Exception as exc:  # recorded; the corresponding criterion then fails
                r = {"error": f"{type(exc).__name__}: {exc}", "seconds": float("nan")}
            if kind in ("gc0", "gc1"):
                jobs[kind] = r
            else:
                jobs[kind][key] = r
            if isinstance(r, dict) and "diag" in r:
                diags[f"{kind}|{key}"] = r["diag"]
            print(f"done {kind} {key} {r.get('seconds', float('nan')):.0f}s"
                  f"{' ERROR ' + r['error'] if 'error' in r else ''}", flush=True)  # fmt: skip
            if kind == "pareto" and "error" not in r:
                nf = pool.submit(job_dense17, (0, np.asarray(r["theta"]), f"ratio_{key:g}", False))
                fut[nf] = ("dense", f"ratio_{key:g}_L0")
                pending.add(nf)
                nf2 = pool.submit(job_ron_point, (np.asarray(r["theta"]), LEVEL_RECHECK,
                                                  f"ratio_{key:g}", False))  # fmt: skip
                fut[nf2] = ("ron_pt", f"ratio_{key:g}_L2")
                pending.add(nf2)
            break
    pool.shutdown()
    T["phase2_parallel"] = time.perf_counter() - t_par
    dense, fields, fl = jobs["dense"], jobs["fields"], jobs["fl"]
    res["jobs"] = {k: v for k, v in jobs.items() if k in ("j1_fd", "j2_fd", "ron_pt", "ron13d")}
    res["dense"] = dense
    res["fields"] = fields
    res["fieldline"] = fl
    res["pareto_runs"] = jobs["pareto"]
    res["grad_checks"] = {"initial": jobs["gc0"], "final": jobs["gc1"]}
    dump()

    # ============================================================ criteria
    crit, vals = {}, {}
    # ---- J1
    g0 = jobs["ron_pt"]["initial_L1"]
    g_adj = np.array(g0["grad"])
    fdv = np.array([jobs["j1_fd"][k]["richardson"] for k in range(5)])
    j1_err = float(np.linalg.norm(g_adj - fdv) / np.linalg.norm(fdv))
    gc0, gf0 = jobs["gc0"]["grad_check"], jobs["gc0"]["grad_check_fields"]
    d2_bv = max(r["bv"]["rel_err"] for r in gc0["directions"])
    d2_e = max(r["E_ox"]["rel_err"] for r in gf0["directions"])
    crit["J1"] = bool(j1_err <= J1_TOL and d2_bv <= D2_TOL_REL and d2_e <= D2_TOL_REL)
    vals["J1"] = {"rel_err_norm_R_on_ch": j1_err, "tol": J1_TOL, "g_adj": _lst(g_adj),
                  "g_fd_richardson": _lst(fdv),
                  "fd_spread_rel_max": max(jobs["j1_fd"][k]["spread"] / abs(fdv[k])
                                           for k in range(5)),
                  "worst_rel_BV_obj_D2": d2_bv, "worst_rel_E_ox_D2": d2_e, "D2_tol": D2_TOL_REL,
                  "R_on_ch_h2_initial": g0["R"]}  # fmt: skip
    # ---- J2 (units of the initial values)
    gc1, gf1 = jobs["gc1"]["grad_check"], jobs["gc1"]["grad_check_fields"]
    g1 = np.array(jobs["ron_pt"]["final_L1"]["grad"])
    r1_h2 = jobs["ron_pt"]["final_L1"]["R"]
    fd1 = np.array([jobs["j2_fd"][k]["richardson"] for k in range(5)])
    a_r = float(np.max(np.abs(g1 - fd1)) / r0)
    n_r = float(np.linalg.norm(g1) / r0)
    n_bv = gc1["grad_bv_norm"] / bv0
    n_e = gf1["grad_E_over_E0_norm"]
    a_bv = max(r["bv"]["abs_diff"] for r in gc1["directions"]) / bv0
    a_e = max(r["E_ox"]["abs_diff"] for r in gf1["directions"]) / e0
    t_ = {k: J2_TOL * max(1.0, n) for k, n in (("bv", n_bv), ("r", n_r), ("e", n_e))}
    crit["J2"] = bool(a_bv <= t_["bv"] and a_r <= t_["r"] and a_e <= t_["e"])
    vals["J2"] = {"worst_abs_BV_over_BV0": a_bv, "tol_BV": t_["bv"], "worst_abs_R_over_R0": a_r,
                  "tol_R": t_["r"], "worst_abs_E_over_E0": a_e, "tol_E": t_["e"],
                  "grad_norms_normalized": [n_bv, n_r, n_e],
                  "R_unit_dir_rel_err_norm": float(np.linalg.norm(g1 - fd1) / np.linalg.norm(fd1))}  # fmt: skip
    # ---- J3
    viol = bounds.violation(th1)
    r1_check = r1_h2
    crit["J3"] = bool(r1_check <= r0 * (1.0 + 1e-8) and r1 <= r0 * (1.0 + 1e-8)
                      and e1 <= e0 * (1.0 + 1e-8) and viol <= 1e-9)  # fmt: skip
    vals["J3"] = {"R_over_R0_minus_1": r1 / r0 - 1.0, "E_over_E0_minus_1": e1 / e0 - 1.0,
                  "bound_violation": viol, "R_final_independent_job": r1_check}  # fmt: skip
    # ---- J4 (every channel on-state evaluation of every job + main)
    g_all = [d["guards"] for d in diags.values() if "guards" in d and d["guards"]["n"]]
    n_g = sum(g["n"] for g in g_all)
    n_bad = sum(g["n_bad"] for g in g_all)
    crit["J4"] = bool(n_g > 0 and n_bad == 0)
    vals["J4"] = {"n_evaluations": n_g, "n_violations": n_bad,
                  "min_gated_ns_cm2": min(g["min_gated_ns_cm2"] for g in g_all),
                  "min_gated_eeff_raw_MV_cm": min(g["min_gated_eeff_raw_MV_cm"]
                                                  for g in g_all),
                  "max_ungated_near_ns_clamp": max(g["max_ungated_near_ns_clamp"]
                                                   for g in g_all),
                  "max_surface_NA_rel_change": max(g["max_surface_NA_rel_change"]
                                                   for g in g_all)}  # fmt: skip
    # ---- J5
    ini = {lev: dense[f"initial_L{lev}"] for lev in (0, 1)}
    fin = {lev: dense[f"final_L{lev}"] for lev in (0, 1)}
    imp = {lev: fin[lev]["bv_min"] - ini[lev]["bv_min"] for lev in (0, 1)}
    crit["J5"] = bool(imp[0] > 0 and imp[1] > 0 and imp[1] >= J5_RATIO * imp[0]
                      and len(ini[0]["x_um"]) >= 40)  # fmt: skip
    vals["J5"] = {"min_bv_initial_h_h2": [ini[0]["bv_min"], ini[1]["bv_min"]],
                  "min_bv_final_h_h2": [fin[0]["bv_min"], fin[1]["bv_min"]],
                  "improvement_h": imp[0], "improvement_h2": imp[1],
                  "ratio": imp[1] / imp[0] if imp[0] != 0 else None,
                  "x_argmin_initial_um": [ini[0]["x_argmin_um"], ini[1]["x_argmin_um"]],
                  "x_argmin_final_um": [fin[0]["x_argmin_um"], fin[1]["x_argmin_um"]]}  # fmt: skip
    # ---- J6
    r_i4, r_f4 = jobs["ron_pt"]["initial_L2"]["R"], jobs["ron_pt"]["final_L2"]["R"]
    crit["J6"] = bool(r_f4 <= r_i4 * (1.0 + J6_TOL))
    vals["J6"] = {"R_ch_h4_initial": r_i4, "R_ch_h4_final": r_f4,
                  "final_over_initial_minus_1": r_f4 / r_i4 - 1.0, "tol": J6_TOL,
                  "R_ch_h2_initial_final": [r0, r1_h2],
                  "h2_to_h4_change_initial": r_i4 / r0 - 1.0, "h2_to_h4_change_final": r_f4 / r1_h2 - 1.0}  # fmt: skip

    # ---- J7
    def agg(path):
        return sum(d_[path[0]].get(path[1], 0) for d_ in diags.values() if path[0] in d_)

    nf_tot = sum(d_["nonfinite_calls"] for d_ in diags.values())
    stderr_text = log_path.read_text(errors="replace") if log_path.exists() else ""
    entry_lines = [ln for ln in stderr_text.splitlines() if "On entry to" in ln]
    failed_picard = agg(("on_state", "failed_picard"))
    failed_evals_on = agg(("on_state", "failed_evals"))
    failed_evals_off = agg(("off_state", "failed_evals"))
    job_errors = [f"{k}" for kind in ("j1_fd", "j2_fd", "pareto", "dense", "fields", "fl",
                                      "ron_pt", "r2", "ron13d") for k, v in jobs[kind].items()
                  if isinstance(v, dict) and "error" in v]  # fmt: skip
    job_errors += [k for k in ("gc0", "gc1") if jobs[k] is None or "error" in jobs[k]]
    nf_adj = sum(d_["NF"]["adjoint"] for d_ in diags.values())
    nf_conv = sum(d_["NF"]["converged_newton"] for d_ in diags.values())
    nf_failed = sum(d_["NF"]["failed_newton"] for d_ in diags.values())
    crit["J7"] = bool(
        nf_tot == nf_failed and
        sl.status == 0 and sl.nit <= args.maxiter and failed_picard == 0 and failed_evals_on == 0
        and failed_evals_off == 0 and nf_adj == 0 and nf_conv == 0 and not entry_lines
        and not job_errors and all(d["resolved_all"] for d in dense.values())
        and all(v.get("converged", True) for v in fl.values())
    )  # fmt: skip
    vals["J7"] = {
        "slsqp_status": int(sl.status), "slsqp_message": sl.message, "nit": int(sl.nit),
        "failed_picard": failed_picard, "failed_evals_on_state": failed_evals_on,
        "failed_evals_off_state": failed_evals_off,
        "nonfinite_calls_total_over_jobs": nf_tot, "nonfinite_in_adjoint_solves": nf_adj,
        "nonfinite_in_converged_newton": nf_conv, "nonfinite_in_failed_newton": nf_failed,
        "nonfinite_outside_newton_trial_steps": nf_tot - nf_failed,
        "failed_newton_calls_rejected_trials": agg(("NF", "failed_newton_calls")),
        "on_state_failed_attempts": agg(("on_state", "failed_attempts_on")),
        "off_state_failed_attempts": agg(("off_state", "failed_attempts")),
        "off_state_warm_fallbacks": agg(("off_state", "warm_fallbacks")),
        "max_picard_iterations": max((d_["on_state"]["max_picard_iterations"]
                                      for d_ in diags.values()), default=0),
        "on_entry_to_lines": len(entry_lines), "stderr_bytes": len(stderr_text),
        "job_errors": job_errors,
        "pareto_status": {str(k): v["status"] for k, v in jobs["pareto"].items()},
        "dense_resolved_all": {k: v["resolved_all"] for k, v in dense.items()},
        "fieldline_converged": {k: v.get("converged", True) for k, v in fl.items()},
        "note": "every non-finite SuperLU call must lie inside a failed (rejected) Newton solve; "
                "those are reported, any other one fails J7",
    }  # fmt: skip

    # ============================================================ reports R1-R3
    # R1
    ch13 = jobs["ron13d"].get(1, {}).get("R", {}) if jobs["ron13d"].get(1) else {}
    ch13_4 = jobs["ron13d"].get(2, {}).get("R", {}) if jobs["ron13d"].get(2) else {}
    par17 = [{"ratio": 1.0, "theta": th1.tolist(), "bv_obj": bv1, "R": r1, "E_ox": e1,
              "status": int(sl.status), "nit": int(sl.nit), "bv_true_min": fin[0]["bv_min"],
              "x_argmin_um": fin[0]["x_argmin_um"]}]  # fmt: skip
    for rt in PARETO_RATIOS:
        q = dict(jobs["pareto"][rt])
        d = dense.get(f"ratio_{rt:g}_L0", {})
        q["bv_true_min"], q["x_argmin_um"] = d.get("bv_min"), d.get("x_argmin_um")
        par17.append(q)
    par17.insert(0, {"ratio": None, "note": "initial", "theta": TH0.tolist(), "bv_obj": bv0,
                     "R": r0, "E_ox": e0, "bv_true_min": ini[0]["bv_min"],
                     "x_argmin_um": ini[0]["x_argmin_um"]})  # fmt: skip
    for q in par17:
        q["R_over_R0"] = q["R"] / r0
        q["E_ox_over_E0"] = q["E_ox"] / e0
    par13 = [{"ratio": 1.0, "label": "d13_final", "R_ch": ch13.get("d13_final"),
              "R_ch_over_R0": ch13.get("d13_final", np.nan) / r0,
              "bv_true_min": dense["d13_L0"]["bv_min"]}]  # fmt: skip
    for q in r13d["R2"]["pareto"]:
        if q.get("ratio") and q["ratio"] > 1.0:
            key = f"ratio_{q['ratio']:g}"
            par13.append(
                {
                    "ratio": q["ratio"],
                    "label": key,
                    "R_ch": ch13.get(key),
                    "R_ch_over_R0": ch13.get(key, np.nan) / r0,
                    "bv_true_min": q["bv_true_min"],
                    "R_bulk_over_R0_bulk": q["R_over_R0"],
                }
            )
    res["R1"] = {
        "theta_17_final": th1.tolist(), "theta_13d_final": th13.tolist(),
        "theta_diff_17_minus_13d": (th1 - th13).tolist(),
        "theta_J_17_vs_13d": [float(th1[4]), float(th13[4])],
        "bv_true_min_h": {"initial": ini[0]["bv_min"], "17": fin[0]["bv_min"],
                          "13d": dense["d13_L0"]["bv_min"]},
        "bv_true_min_h2": {"initial": ini[1]["bv_min"], "17": fin[1]["bv_min"],
                           "13d": dense["d13_L1"]["bv_min"]},
        "bv_gain_h_17": imp[0], "bv_gain_h_13d": dense["d13_L0"]["bv_min"] - ini[0]["bv_min"],
        "R_ch_13d_final_h2": ch13.get("d13_final"), "R_ch_13d_final_h4": ch13_4.get("d13_final"),
        "R_ch_13d_final_over_R0_minus_1_h2": ch13.get("d13_final", np.nan) / r0 - 1.0,
        "R_ch_17_final_h2": r1_h2, "pareto_17": par17, "pareto_13d_channel_R": par13,
        "bv_obj_13d_final_on_level0": r13d["final"]["bv_obj"],
    }  # fmt: skip
    # R2
    r2i, r2f = jobs["r2"].get("initial"), jobs["r2"].get("final")
    res["R2"] = {"initial": r2i, "final": r2f}

    # R3
    def fld(key):
        return fields[key]["fields"]

    r3 = {"fields_V_blk": {f"{t}_L{lv}": fld(f"{t}_L{lv}_blk") for t in ("initial", "d13", "final")
                           for lv in (0, 1)},
          "fields_V_blk_h4": {t: fld(f"{t}_L2_blk") for t in ("initial", "final") if
                              f"{t}_L2_blk" in fields},  # fmt: skip
          "same_mesh_rel_change_17_vs_initial": {
              f"L{lv}": {k: fld(f"final_L{lv}_blk")[k] / fld(f"initial_L{lv}_blk")[k] - 1.0
                         for k in MEASURES} for lv in (0, 1)},
          "same_mesh_rel_change_17_vs_13d": {
              f"L{lv}": {k: fld(f"final_L{lv}_blk")[k] / fld(f"d13_L{lv}_blk")[k] - 1.0
                         for k in MEASURES} for lv in (0, 1)},
          "fieldline_bv_V": {k: v.get("bv") for k, v in fl.items()},
          "fieldline_vs_dense_cut_min": {
              k: {"fieldline": fl[k].get("bv"),
                  "dense_cut_min": dense[k]["bv_min"] if k in dense else None} for k in fl}}  # fmt: skip
    res["R3"] = r3
    res["diagnostics"] = diags
    res["runtime_seconds"] = time.perf_counter() - t_start
    res["seconds"] = T
    res["rss_mib_main"] = rss_mib()

    figdir = (REPO / "figures") if out == DEFAULT_OUT.resolve() else out / "figures"
    label = f"stage 17 | commit {str(man.get('git_commit'))[:8]} | {man['timestamp_utc'][:19]}Z"
    fig_err = None
    try:
        plot_all(figdir, res, label)
    except Exception as exc:  # figures are optional
        fig_err = f"{type(exc).__name__}: {exc}"
    res["figures"] = {"dir": str(figdir), "error": fig_err}
    summary = {
        "criteria": crit, "gate_passed": bool(all(crit.values())), "dev": bool(args.dev),
        "values": {
            "bv_obj_initial_V": bv0, "bv_obj_final_V": bv1,
            "R_ch_h2_initial": r0, "R_ch_h2_final": r1_h2, "E_ox_initial_MV_cm": e0,
            "E_ox_final_MV_cm": e1, **vals, "R1": res["R1"],
            "R2_lagged_error": {t: (jobs["r2"][t]["lagged_rel_error_norm"]
                                    if t in jobs["r2"] and "lagged_rel_error_norm" in jobs["r2"][t]
                                    else None) for t in ("initial", "final")},
            "R3": r3,
        },
        "seconds": T, "runtime_seconds": res["runtime_seconds"],
    }  # fmt: skip
    dump()
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    print(json.dumps(crit, indent=2))
    print(f"runtime {res['runtime_seconds']:.1f} s")


if __name__ == "__main__":
    main()
