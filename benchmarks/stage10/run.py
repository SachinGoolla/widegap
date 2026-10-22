"""Stage-10 benchmark (docs/benchmarks/stage10_mosfet_optimization.md): O1-O7 (gating), R1-R3.

2D adjoint doping sensitivities of the drain current I_D and of a smoothed field-based
breakdown voltage BV_s of a vertical DMOSFET half-cell, and SLSQP optimization of the
drift-layer profile + JFET level (max BV_s s.t. I_D >= I_D,0).
Writes <out>/{results,summary,manifest}.json (default results/stage10) and the figures
stage10_sensitivity_maps.png, stage10_optimized_profile.png (repo figures/ for the default
output, else <out>/figures).
"""

from __future__ import annotations

import argparse
import json
import resource
import sys
import time
import warnings
from pathlib import Path

import numpy as np
from numpy.polynomial import legendre
from scipy.optimize import minimize

from widegap.adjoint2d import (
    FieldIndex,
    breakdown_voltage,
    bv_gradient_c,
    drain_current_grad_c,
)
from widegap.devices.dmosfet import DMOSFET, DMOSParams
from widegap.fv2d.mos import continue_bias
from widegap.manifest import build_manifest

warnings.filterwarnings("ignore", category=RuntimeWarning)  # transient overflow in damped Newton

REPO = Path(__file__).resolve().parents[2]

# ---- study parameters (all recorded in results.json) -------------------------------------
PRM = DMOSParams(n0_drift=16.0, n_jfet=18.0)
E_CRIT, KAPPA = 75.0, 200.0  # field-index scale [V_T/L_D] and stage-6 smoothing
V_ON, V_LIN = 200.0, 2.0  # on-state bias (V_GS, V_DS), scaled V_T
V_STEP, V_MAX = 25.0, 600.0  # drain-bias marching step and ceiling for the BV search
N_DRIFT_BOUNDS = (4.0, 64.0)  # pointwise bounds on N_D(y)/N0 over the drift depth
N_JFET_BOUNDS = (6.0, 72.0)  # JFET donor level bounds
THETA_BOX = 3.0  # |theta_k| <= THETA_BOX for the Legendre coefficients
PARETO_RATIOS = (0.98, 0.95, 0.90)  # looser I_D >= ratio * I_D,0 levels
SEEDS = (0, 1, 2)
EPS_FD = 1e-3
FEAS_TOL = 1e-8  # O4: I_D >= I_D0 (1 - 1e-8)
MAX_ITER = 60
TH0 = np.zeros(5)


def rss_mib():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


# ============================================================ evaluator with warm starts
class Evaluator:
    """BV_s and I_D (+ adjoint gradients) of the DMOSFET at theta, on mesh ``level``.
    Newton starting points are reused from the previously evaluated (nearby) design."""

    def __init__(self, level=0, prm=PRM):
        self.dev = DMOSFET(prm, level)
        mask = self.dev.blocking_edge_mask()
        self.idx = FieldIndex(E_CRIT, mask, KAPPA)
        self.idx_exact = FieldIndex(E_CRIT, mask, None)
        self.X_eq = None
        self.X_on = None
        self.bv_cache: list = []
        self.n_bv = self.n_id = 0
        self.log: list = []

    def eq(self, th):
        mod, c = self.dev.model(th)
        bc = self.dev.bc(c, 0.0, 0.0, 0.0)
        eqbc = mod.equilibrium_bc(bc)
        for guess in (self.X_eq, None):
            X0 = mod.initial_guess(bc) if guess is None else guess
            s = mod.solve_bc(X0, eqbc)
            if s.report.converged:
                self.X_eq = s.X.copy()
                return mod, c, s
        raise RuntimeError("equilibrium solve did not converge")

    def i_d(self, th, grad=True):
        mod, c, s_eq = self.eq(th)
        dev = self.dev
        bc = dev.bc(c, V_ON, 0.0, V_LIN)
        s = mod.solve_bc(self.X_on, bc) if self.X_on is not None else None
        if s is None or not s.report.converged:
            out = continue_bias(
                mod,
                s_eq.X,
                lambda p: dev.bc(c, p[0], 0.0, p[1]),
                [0.0, 0.0],
                [[0.0, V_LIN], [V_ON, V_LIN]],
            )
            s = out[-1][1]
            if not s.report.converged:
                raise RuntimeError("on-state solve did not converge")
        self.X_on = s.X.copy()
        self.n_id += 1
        cur = mod.contact_current(s, dev.nodes["drain"])
        if not grad:
            return cur, None, None
        _, g_c = drain_current_grad_c(mod, s.X, bc, dev.nodes["drain"])
        return cur, dev.doping(th)[1].T @ g_c, g_c

    def bv_s(self, th, grad=True, exact=False, warm=True):
        mod, c, s_eq = self.eq(th)
        dev = self.dev
        idx = self.idx_exact if exact else self.idx
        r = breakdown_voltage(
            mod,
            lambda v: dev.bc(c, 0.0, 0.0, v),
            idx,
            s_eq.X,
            guesses=self.bv_cache if warm else (),
            v_step=V_STEP,
            vmax=V_MAX,
        )
        if not exact and warm:
            self.bv_cache = r.cache
        self.n_bv += 1
        if not grad:
            return r.bv, None, None
        g_c, _, _ = bv_gradient_c(mod, r.X, r.bc, dev.nodes["drain"], idx)
        return r.bv, dev.doping(th)[1].T @ g_c, g_c

    def both(self, th):
        bv, gbv, _ = self.bv_s(th)
        i, gi, _ = self.i_d(th)
        self.log.append({"theta": np.asarray(th).tolist(), "bv_s": bv, "i_d": i})
        return bv, gbv, i, gi


# ============================================================ design-variable constraints
def drift_basis(n=33):
    eta = np.linspace(-1.0, 1.0, n)
    return np.stack([legendre.legval(eta, np.eye(4)[k]) for k in range(4)], axis=1)  # (n, 4)


def profile_constraints():
    """Linear inequalities A theta + b >= 0 for N_lo <= N_D(y) <= N_hi (log scale)."""
    P = drift_basis()
    lo, hi = np.log(np.array(N_DRIFT_BOUNDS) / PRM.n0_drift)
    A = np.vstack([np.hstack([P, np.zeros((len(P), 1))]), np.hstack([-P, np.zeros((len(P), 1))])])
    b = np.r_[-lo * np.ones(len(P)), hi * np.ones(len(P))]
    return A, b


def bounds_list():
    lo = np.log(N_JFET_BOUNDS[0] / PRM.n_jfet)
    hi = np.log(N_JFET_BOUNDS[1] / PRM.n_jfet)
    return [(-THETA_BOX, THETA_BOX)] * 4 + [(lo, hi)]


def bound_violation(th):
    """Largest violation of the box + profile bounds (<= 0 means feasible)."""
    th = np.asarray(th, float)
    A, b = profile_constraints()
    v = [float(np.max(-(A @ th + b)))]
    for t, (lo, hi) in zip(th, bounds_list(), strict=True):
        v += [lo - t, t - hi]
    return max(v)


# ============================================================ optimization
class Problem:
    """Shares one (BV_s, I_D, gradients) evaluation between objective and constraints."""

    def __init__(self, ev: Evaluator, bv0, i0, ratio):
        self.ev, self.bv0, self.i0, self.ratio = ev, bv0, i0, ratio
        self.last = (None, None)
        self.history: list = []
        self.max_bound_violation = -np.inf

    def _get(self, th):
        key = np.asarray(th, float).tobytes()
        if self.last[0] != key:
            self.max_bound_violation = max(self.max_bound_violation, bound_violation(th))
            r = self.ev.both(np.asarray(th, float))
            self.last = (key, r)
        return self.last[1]

    def fun(self, th):
        bv, gbv, _, _ = self._get(th)
        return -bv / self.bv0

    def jac(self, th):
        bv, gbv, _, _ = self._get(th)
        return -gbv / self.bv0

    def c_fun(self, th):
        _, _, i, _ = self._get(th)
        return i / self.i0 - self.ratio

    def c_jac(self, th):
        _, _, i, gi = self._get(th)
        return gi / self.i0


def run_slsqp(ev, th_start, bv0, i0, ratio, maxiter=MAX_ITER, ftol=1e-9):
    A, b = profile_constraints()
    pr = Problem(ev, bv0, i0, ratio)
    cons = [
        {"type": "ineq", "fun": pr.c_fun, "jac": pr.c_jac},
        {"type": "ineq", "fun": lambda t: A @ t + b, "jac": lambda t: A},
    ]
    it_log = []

    def cb(xk):
        it_log.append(np.asarray(xk).tolist())

    res = minimize(
        pr.fun,
        np.asarray(th_start, float),
        jac=pr.jac,
        method="SLSQP",
        bounds=bounds_list(),
        constraints=cons,
        options={"maxiter": maxiter, "ftol": ftol},
        callback=cb,
    )
    th = np.asarray(res.x, float)
    return th, res, pr, it_log


# ============================================================ checks
def directions():
    out = []
    for s in SEEDS:
        u = np.random.default_rng(s).standard_normal(5)
        out.append(u / np.linalg.norm(u))
    return out


def grad_check(ev: Evaluator, th):
    """Adjoint vs central FD along the seeded directions.  The gated FD uses EPS_FD = 1e-3;
    report-only extras: FD at eps/2, eps/10, eps/100 and the Richardson value
    (4 FD(eps/2) - FD(eps)) / 3, to separate adjoint bugs from FD truncation error."""
    th = np.asarray(th, float)
    bv, gbv, _ = ev.bv_s(th)
    i, gi, _ = ev.i_d(th)
    rows = []
    epss = (EPS_FD, EPS_FD / 2, EPS_FD / 10, EPS_FD / 100)
    for s, u in zip(SEEDS, directions(), strict=True):
        fd_b, fd_i = {}, {}
        for e in epss:
            fb = [ev.bv_s(th + sg * e * u, grad=False)[0] for sg in (1, -1)]
            fi = [ev.i_d(th + sg * e * u, grad=False)[0] for sg in (1, -1)]
            fd_b[e], fd_i[e] = (fb[0] - fb[1]) / (2 * e), (fi[0] - fi[1]) / (2 * e)
        rb = (4 * fd_b[epss[1]] - fd_b[epss[0]]) / 3
        ri = (4 * fd_i[epss[1]] - fd_i[epss[0]]) / 3
        ab, ai = float(gbv @ u), float(gi @ u)
        rows.append(
            {
                "seed": s,
                "bv_adjoint": ab, "bv_fd": fd_b[EPS_FD],
                "bv_abs_diff": abs(ab - fd_b[EPS_FD]), "bv_rel_err": abs(ab - fd_b[EPS_FD]) / abs(fd_b[EPS_FD]),
                "id_adjoint": ai, "id_fd": fd_i[EPS_FD],
                "id_abs_diff": abs(ai - fd_i[EPS_FD]), "id_rel_err": abs(ai - fd_i[EPS_FD]) / abs(fd_i[EPS_FD]),
                "report_only": {
                    "fd_by_eps_bv": {str(e): v for e, v in fd_b.items()},
                    "fd_by_eps_id": {str(e): v for e, v in fd_i.items()},
                    "bv_richardson": rb, "bv_abs_diff_richardson": abs(ab - rb),
                    "bv_abs_diff_smallest_eps": abs(ab - fd_b[epss[-1]]),
                    "id_richardson": ri, "id_abs_diff_richardson": abs(ai - ri),
                    "id_abs_diff_smallest_eps": abs(ai - fd_i[epss[-1]]),
                },
            }
        )  # fmt: skip
    # re-establish the base-point warm state
    ev.bv_s(th, grad=False)
    ev.i_d(th, grad=False)
    return {
        "bv_s": bv, "i_d": i,
        "grad_bv_norm": float(np.linalg.norm(gbv)), "grad_id_norm": float(np.linalg.norm(gi)),
        "grad_bv": gbv.tolist(), "grad_id": gi.tolist(), "directions": rows,
    }  # fmt: skip


def uniform_grid(ev: Evaluator, i0):
    """O6 family: theta_{1..3} = 0, theta_J = 0, 21 points in theta_0 over the allowed range."""
    lim = float(np.log(N_DRIFT_BOUNDS[1] / PRM.n0_drift))
    rows = []
    for t0 in np.linspace(-lim, lim, 21):
        th = np.zeros(5)
        th[0] = t0
        entry = {"theta0": float(t0), "N_D": float(PRM.n0_drift * np.exp(t0))}
        try:
            i = ev.i_d(th, grad=False)[0]
            entry["i_d"] = i
            entry["feasible"] = bool(i >= i0 * (1.0 - FEAS_TOL))
            if entry["feasible"]:
                entry["bv_s"] = ev.bv_s(th, grad=False)[0]
            else:
                entry["bv_s"] = None
                entry["skipped"] = "infeasible (I_D < I_D,0); BV_s not needed for O6"
        except (RuntimeError, ValueError) as exc:
            entry["error"] = f"{type(exc).__name__}: {exc}"
            entry["feasible"] = False
        rows.append(entry)
    return rows


def reverify(th0, th1, levels=(0, 1)):
    """Fresh, cold-started evaluators on mesh h (level 0) and h/2 (level 1): exact BV, I_D."""
    out = {}
    for lev in levels:
        ev = Evaluator(lev)
        t0 = time.perf_counter()
        b0 = ev.bv_s(th0, grad=False, exact=True, warm=False)[0]
        b1 = ev.bv_s(th1, grad=False, exact=True, warm=False)[0]
        i0 = ev.i_d(th0, grad=False)[0]
        i1 = ev.i_d(th1, grad=False)[0]
        s0 = ev.bv_s(th0, grad=False, warm=False)[0]
        s1 = ev.bv_s(th1, grad=False, warm=False)[0]
        out[str(lev)] = {
            "level": lev, "cells": [ev.dev.mesh.nx, ev.dev.mesh.ny],
            "bv_exact_initial": b0, "bv_exact_final": b1, "improvement_exact": b1 - b0,
            "bv_s_initial": s0, "bv_s_final": s1, "i_d_initial": i0, "i_d_final": i1,
            "seconds": time.perf_counter() - t0,
        }  # fmt: skip
    return out


# ============================================================ figures
def sensitivity_maps(ev: Evaluator, th):
    """gamma(x, y) = dJ/dC_i / vol_s,i (per unit area) for BV_s and I_D, ZA19/BA26 convention
    D = N_A - N_D, i.e. gamma_D = -gamma_C."""
    _, _, g_bv = ev.bv_s(th)
    _, _, g_id = ev.i_d(th)
    vol = ev.dev.reg.vol_s
    ok = vol > 0
    gam = {}
    for name, g in (("BV_s", g_bv), ("I_D", g_id)):
        gd = np.zeros_like(g)
        gd[ok] = -g[ok] / vol[ok]
        gam[name] = gd
    return gam


def plot_maps(ev, gam, out_png, man):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm

    m, reg = ev.dev.mesh, ev.dev.reg
    shape = (m.nx + 1, m.ny + 1)
    X, Y = np.meshgrid(m.x, m.y, indexing="ij")
    fig, ax = plt.subplots(1, 2, figsize=(12, 5.2), constrained_layout=True)
    for a, (name, g) in zip(ax, gam.items(), strict=True):
        gm = np.where(reg.semi_node, g, np.nan)
        lin = np.nanpercentile(np.abs(gm[np.isfinite(gm)]), 20) + 1e-300
        vmax = np.nanmax(np.abs(gm))
        mp = a.pcolormesh(
            X, Y, np.ma.masked_invalid(gm.reshape(shape)), shading="gouraud", cmap="RdBu_r",
            norm=SymLogNorm(linthresh=lin, vmin=-vmax, vmax=vmax),
        )  # fmt: skip
        fig.colorbar(mp, ax=a, label=f"gamma_{name}  (d{name}/dD per unit area, D = N_A - N_D)")
        a.set_ylim(m.y[0] if name == "BV_s" else -14.0, 0.0)
        a.set_title(f"{name} sensitivity to D = N_A - N_D, initial design")
        a.set_xlabel("x  [L_D]")
        a.set_ylabel("y  [L_D]")
    fig.text(0.5, 0.005, man, ha="center", fontsize=7, color="0.4")
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


def plot_profile(ev, th0, th1, pareto, out_png, man):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    dev = ev.dev
    p = dev.prm
    ys = np.linspace(-(p.d_body + p.w_drift), -p.d_body, 200)
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2), constrained_layout=True)
    for th, lab, col in ((th0, "initial", "#0072B2"), (th1, "optimized", "#D55E00")):
        ax[0].plot(
            dev.drift_profile(th, ys),
            ys,
            color=col,
            label=f"{lab}  (JFET {dev.jfet_level(th):.1f})",
        )
        ev.bv_s(th, grad=False)  # warm the caches (also used for the field profile below)
        mod, c, s = ev.eq(th)
        r = breakdown_voltage(
            mod, lambda v, c=c: dev.bc(c, 0.0, 0.0, v), ev.idx, s.X, guesses=ev.bv_cache,
            v_step=V_STEP, vmax=V_MAX,
        )  # fmt: skip
        e = ev.idx.edge_fields(mod, r.X)
        xy = dev.mesh.node_xy
        sel = (dev.mesh.edge_dir == 1) & (xy[mod.a, 0] == p.x_max) & (xy[mod.a, 1] < -p.d_body + 3)
        yc = 0.5 * (xy[mod.a, 1] + xy[mod.b, 1])
        order = np.argsort(yc[sel])
        ax[1].plot(
            np.abs(e[sel])[order], yc[sel][order], color=col, label=f"{lab}: BV_s = {r.bv:.1f}"
        )
    ax[0].set_xscale("log")
    ax[0].set_xlabel("drift donor N_D(y)/N0")
    ax[0].set_ylabel("y  [L_D]")
    ax[0].set_title("Drift profile (body junction at top)")
    ax[0].legend(fontsize=8)
    ax[1].axvline(E_CRIT, color="0.5", ls="--", lw=0.8)
    ax[1].set_xlabel("|E| at BV_s along x = x_max  [V_T/L_D]")
    ax[1].set_ylabel("y  [L_D]")
    ax[1].set_title("Field at the breakdown state")
    ax[1].legend(fontsize=8)
    pb = [q["bv_s"] for q in pareto]
    pi = [q["i_d"] for q in pareto]
    ax[2].plot(pi, pb, "o-", color="#009E73")
    ax[2].set_xlabel(f"I_D  (V_GS = {V_ON:g}, V_DS = {V_LIN:g})")
    ax[2].set_ylabel("BV_s  [V_T]")
    ax[2].set_title("Pareto set (R1)")
    fig.text(0.5, 0.005, man, ha="center", fontsize=7, color="0.4")
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


# ============================================================ main
def main() -> None:
    default_out = REPO / "results" / "stage10"
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=default_out)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    figdir = REPO / "figures" if args.out == default_out else args.out / "figures"
    manifest = build_manifest(REPO, sys.argv, list(SEEDS))
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    stamp = (
        f"widegap | commit {str(manifest['git_commit'])[:7]} | {manifest['cpu_model']} | CPU only"
    )

    T = {}

    def timed(name, fn, *a, **k):
        t0 = time.perf_counter()
        r = fn(*a, **k)
        T[name] = time.perf_counter() - t0
        print(f"{name}: {T[name]:.1f}s", flush=True)
        return r

    ev = Evaluator(0)
    dev = ev.dev
    config = {
        "device": {k: (list(v) if isinstance(v, tuple) else v) for k, v in vars(PRM).items()},
        "mesh_level0": [dev.mesh.nx, dev.mesh.ny],
        "mesh_level1": [2 * dev.mesh.nx, 2 * dev.mesh.ny],
        "E_crit": E_CRIT, "kappa": KAPPA, "V_on": V_ON, "V_lin": V_LIN,
        "v_step": V_STEP, "v_max": V_MAX,
        "N_drift_bounds": N_DRIFT_BOUNDS, "N_jfet_bounds": N_JFET_BOUNDS, "theta_box": THETA_BOX,
        "eps_fd": EPS_FD, "seeds": list(SEEDS), "pareto_ratios": PARETO_RATIOS,
        "n_blocking_edges": int(ev.idx.mask.sum()), "n_semiconductor_edges": int((dev.reg.w_s > 0).sum()),
        "index_note": "B_s over semiconductor edges except the zero-bias source-junction box "
        "(x in [x_src0-3.5, x_src1+2], y > -(d_src+3.5)); see DMOSFET.blocking_edge_mask",
    }  # fmt: skip

    # ---- initial design: BV_s, I_D, gradients, O1/O2
    bv0, _, _ = timed("initial_bv_s", ev.bv_s, TH0, grad=False)
    i0, _, _ = timed("initial_i_d", ev.i_d, TH0, grad=False)
    gc0 = timed("grad_check_initial", grad_check, ev, TH0)

    # ---- optimization (main problem, ratio 1)
    th1, res, pr, it_log = timed("slsqp_main", run_slsqp, ev, TH0, bv0, i0, 1.0)
    bv1 = ev.bv_s(th1, grad=False)[0]
    i1 = ev.i_d(th1, grad=False)[0]
    print(f"main: BV_s {bv0:.4f} -> {bv1:.4f}, I_D {i0:.6f} -> {i1:.6f}, {res.message}", flush=True)
    gc1 = timed("grad_check_final", grad_check, ev, th1)

    # ---- O6 nested uniform family (before the Pareto runs; they re-warm anyway)
    grid = timed("uniform_grid", uniform_grid, ev, i0)
    feas = [g["bv_s"] for g in grid if g.get("feasible") and "bv_s" in g]
    best_uniform = max(feas) if feas else float("nan")

    # ---- R1 Pareto
    pareto = [
        {"ratio": 1.0, "bv_s": bv1, "i_d": i1, "theta": th1.tolist(), "status": int(res.status),
         "nit": int(res.nit)},
    ]  # fmt: skip
    th_prev = th1
    for r in PARETO_RATIOS:
        th_r, rs, prr, _ = timed(f"slsqp_ratio_{r}", run_slsqp, ev, th_prev, bv0, i0, r, 40)
        b_r, i_r = ev.bv_s(th_r, grad=False)[0], ev.i_d(th_r, grad=False)[0]
        pareto.append(
            {"ratio": r, "bv_s": b_r, "i_d": i_r, "theta": th_r.tolist(), "status": int(rs.status),
             "nit": int(rs.nit), "message": rs.message}
        )  # fmt: skip
        th_prev = th_r
    pareto.append({"ratio": None, "bv_s": bv0, "i_d": i0, "theta": TH0.tolist(), "note": "initial"})

    # ---- O5 independent re-verification on h and h/2
    rv = timed("reverify", reverify, TH0, th1)
    pareto_exact = []
    for q in pareto:
        if q.get("theta") is not None:
            pareto_exact.append(float(ev.bv_s(np.asarray(q["theta"]), grad=False, exact=True)[0]))
    for q, be in zip([q for q in pareto if q.get("theta") is not None], pareto_exact, strict=True):
        q["bv_exact_level0"] = be

    # ---- R2 maps and R3 trends at the initial design
    gam = sensitivity_maps(ev, TH0)
    fig_err = None
    try:
        plot_maps(ev, gam, figdir / "stage10_sensitivity_maps.png", stamp)
        plot_profile(ev, TH0, th1, [q for q in pareto if q.get("i_d") is not None],
                     figdir / "stage10_optimized_profile.png", stamp)  # fmt: skip
    except Exception as exc:  # figures are optional
        fig_err = f"skipped: {type(exc).__name__}: {exc}"
    ys = np.linspace(-(PRM.d_body + PRM.w_drift), -PRM.d_body, 7)
    r3 = {
        "y_samples": ys.tolist(),
        "N_D_initial": dev.drift_profile(TH0, ys).tolist(),
        "N_D_final": dev.drift_profile(th1, ys).tolist(),
        "jfet_initial": dev.jfet_level(TH0),
        "jfet_final": dev.jfet_level(th1),
        "pareto_jfet": [dev.jfet_level(np.asarray(q["theta"])) for q in pareto if q.get("theta")],
        "gamma_BV_drift_mean_near_junction": float(
            np.mean(
                gam["BV_s"][
                    (dev.mesh.node_xy[:, 1] < -PRM.d_body) & (dev.mesh.node_xy[:, 1] > -PRM.d_body - 3)
                    & (dev.mesh.node_xy[:, 0] > PRM.x_jfet)
                ]
            )
        ),
        "note": "sign convention D = N_A - N_D: gamma_BV < 0 near the body junction means lowering "
        "N_D (raising D) there raises BV_s.",
    }  # fmt: skip

    # ---- criteria
    o1 = all(r["id_rel_err"] <= 1e-5 for r in gc0["directions"])
    o2 = all(r["bv_rel_err"] <= 1e-4 for r in gc0["directions"])
    # O3' (spec section 3): adjoint vs Richardson-extrapolated central FD
    o3 = all(
        r["report_only"]["id_abs_diff_richardson"] <= 1e-6 * max(1.0, gc1["grad_id_norm"])
        and r["report_only"]["bv_abs_diff_richardson"] <= 1e-6 * max(1.0, gc1["grad_bv_norm"])
        for r in gc1["directions"]
    )
    viol = bound_violation(th1)
    o4 = bool(i1 >= i0 * (1.0 - FEAS_TOL) and viol <= 1e-9)
    l0, l1 = rv["0"], rv["1"]
    o5 = bool(
        l0["improvement_exact"] > 0
        and l1["improvement_exact"] > 0
        and l1["improvement_exact"] >= 0.5 * l0["improvement_exact"]
        and abs(l1["i_d_final"] - l0["i_d_final"]) <= 0.02 * abs(l0["i_d_final"])
    )
    o6 = bool(bv1 >= best_uniform - 1e-6 * bv1)
    o7 = bool(res.status == 0 and res.nit <= MAX_ITER)
    crit = {"O1": o1, "O2": o2, "O3": o3, "O4": o4, "O5": o5, "O6": o6, "O7": o7}
    crit = {k: bool(v) for k, v in crit.items()}

    results = {
        "config": config,
        "initial": {"theta": TH0.tolist(), "bv_s": bv0, "i_d": i0},
        "final": {
            "theta": th1.tolist(),
            "bv_s": bv1,
            "i_d": i1,
            "bound_violation": viol,
            "max_bound_violation_during_run": pr.max_bound_violation,
        },  # fmt: skip
        "grad_check_initial": gc0,
        "grad_check_final": gc1,
        "slsqp_main": {
            "status": int(res.status),
            "message": res.message,
            "nit": int(res.nit),
            "nfev": int(res.nfev),
            "njev": int(res.njev),
            "iterates": it_log,
        },  # fmt: skip
        "evaluations_main": ev.log,
        "uniform_grid": grid,
        "best_uniform_bv_s": best_uniform,
        "reverify": rv,
        "R1_pareto": pareto,
        "R2": {
            "figure": str(figdir / "stage10_sensitivity_maps.png"),
            "figure_error": fig_err,
            "node_x": dev.mesh.node_xy[:, 0].tolist(),
            "node_y": dev.mesh.node_xy[:, 1].tolist(),
            "gamma_BV_s_D": gam["BV_s"].tolist(),
            "gamma_I_D_D": gam["I_D"].tolist(),
        },  # fmt: skip
        "R3": r3,
    }
    summary = {
        "criteria": crit,
        "gate_passed": all(crit.values()),
        "values": {
            "bv_s_initial": bv0, "bv_s_final": bv1, "i_d_initial": i0, "i_d_final": i1,
            "O1_worst_rel_err_I_D": max(r["id_rel_err"] for r in gc0["directions"]),
            "O2_worst_rel_err_BV_s": max(r["bv_rel_err"] for r in gc0["directions"]),
            "O3_worst_abs_diff_I_D": max(r["id_abs_diff"] for r in gc1["directions"]),
            "O3_worst_abs_diff_BV_s": max(r["bv_abs_diff"] for r in gc1["directions"]),
            "O3_grad_norms": [gc1["grad_id_norm"], gc1["grad_bv_norm"]],
            "O3_report_only_worst_abs_diff_vs_richardson": [
                max(r["report_only"]["id_abs_diff_richardson"] for r in gc1["directions"]),
                max(r["report_only"]["bv_abs_diff_richardson"] for r in gc1["directions"]),
            ],
            "O3_report_only_worst_abs_diff_vs_eps_1e-5": [
                max(r["report_only"]["id_abs_diff_smallest_eps"] for r in gc1["directions"]),
                max(r["report_only"]["bv_abs_diff_smallest_eps"] for r in gc1["directions"]),
            ],
            "O4_i_d_over_i_d0_minus_1": i1 / i0 - 1.0, "O4_bound_violation": viol,
            "O5": {"improvement_h": l0["improvement_exact"], "improvement_h2": l1["improvement_exact"],
                   "ratio": l1["improvement_exact"] / l0["improvement_exact"],
                   "i_d_final_h": l0["i_d_final"], "i_d_final_h2": l1["i_d_final"],
                   "i_d_rel_diff": abs(l1["i_d_final"] - l0["i_d_final"]) / abs(l0["i_d_final"])},
            "O6": {"bv_s_final": bv1, "best_uniform_bv_s": best_uniform,
                   "margin": bv1 - best_uniform},
            "O7": {"status": int(res.status), "nit": int(res.nit), "message": res.message},
            "R1": [{k: q.get(k) for k in ("ratio", "bv_s", "i_d", "bv_exact_level0")} for q in pareto],
            "R3": {k: r3[k] for k in ("N_D_initial", "N_D_final", "jfet_initial", "jfet_final",
                                      "gamma_BV_drift_mean_near_junction")},
        },
        "evaluations": {"bv": ev.n_bv, "i_d": ev.n_id},
        "seconds": T,
        "peak_rss_mib": rss_mib(),
    }  # fmt: skip
    (args.out / "results.json").write_text(json.dumps(results, indent=2, default=float) + "\n")
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    print(json.dumps(summary, indent=2, default=float))


if __name__ == "__main__":
    main()
