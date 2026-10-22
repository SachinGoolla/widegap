"""Stage-13 benchmark (docs/benchmarks/stage13_sic_dmosfet_optimization.md): T1-T7 (gating), R1-R3.

Full-physics adjoint optimization of the stage-12 4H-SiC DMOSFET half-cell (sourced registry
parameters, incomplete ionization, Caughey-Thomas mobility incl. doping dependence and field
dependence, interface charge): maximize the avalanche BV (ionization integral along the FIXED
vertical column x = x_bc, Niwa coefficients) subject to R_on,sp <= R_on,sp,0 at V_GS = 15 V,
V_DS = 0.1 V, SLSQP on adjoint gradients (widegap.adjoint_sic).

Writes <out>/{results,summary,manifest}.json (default results/stage13) and the figures
stage13_sensitivity_maps.png, stage13_pareto.png, stage13_optimized_profile.png (repo figures/
for the default output, else <out>/figures).

Study choices (recorded in results.json; none is a material parameter):

* device: the stage-12 12b design (imported from benchmarks/stage12/run.py), mesh level 0 for
  the optimization (1107 nodes) and fresh cold evaluators on levels 0 and 1 for T5.
* theta bounds: drift N_D(y) in [1e15, 1e17] cm^-3 imposed pointwise at 65 equispaced eta values
  plus the cell-centre eta of both meshes (dense linear inequalities on theta); JFET N_D in
  [1e16, 5e17] (box on theta_J); |theta_k| <= 6 (never active, only bounds SLSQP).
* Pareto levels: R_on <= {1.00, 1.05, 1.10, 1.20} R_on,0 (spec: "three looser levels"; the
  levels are a study choice), each warm-started from the previous solution, <= 40 iterations
  for the looser ones.
* gradient units: BV in V, R_on,sp in ohm cm^2; the T3 tolerance max(1, ||grad||) is
  unit-dependent, so T3 for R_on is gated on R/R_0 (dimensionless) and the literal
  ohm cm^2 numbers are reported too (spec issue, see the report).
* R1 unipolar limit: for each N_D the ideal non-punch-through (abrupt one-sided, triangular
  field) drift layer has E_m, W_dep, BV from ``triangular_breakdown`` (ionization integral I_n=1
  with the registry Niwa coefficients), and
      R_on,sp,ideal = W_dep / (q mu_n,par(N_D) N_D),   W_dep = 2 BV / E_m,
  with mu_n,par the registry Ishikawa Caughey-Thomas fit (Hall, E || c).  Algebraically equal
  to the textbook 4 BV^2 / (eps_s mu_n E_c^3) with E_c = E_m (checked).  No fitting.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import resource
import sys
import time
import warnings
from pathlib import Path

import numpy as np
from numpy.polynomial import legendre
from scipy.optimize import minimize

from widegap.adjoint_sic import EvalError, SiCAdjoint, richardson_fd
from widegap.manifest import build_manifest
from widegap.materials.ionization import Q
from widegap.sic import SiC4H, SiCDMOS, TracedRegistry, triangular_breakdown

warnings.filterwarnings("ignore")  # transient overflow in damped Newton

REPO = Path(__file__).resolve().parents[2]
DEFAULT_OUT = REPO / "results" / "stage13"

_spec = importlib.util.spec_from_file_location("stage12_run", REPO / "benchmarks/stage12/run.py")
S12 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(S12)
DESIGN, FLAGS12 = S12.DESIGN, S12.FLAGS

# ---- study parameters ----------------------------------------------------------------------
ND_DRIFT_BOUNDS = (1e15, 1e17)  # cm^-3 (spec)
ND_JFET_BOUNDS = (1e16, 5e17)
THETA_BOX = 6.0
PARETO_RATIOS = (1.05, 1.10, 1.20)  # looser R_on levels (R_on <= ratio R_on,0)
SEEDS = (0, 1, 2)
EPS_FD = 1e-3
FEAS_TOL = 1e-8
MAX_ITER = 60
N_GRID = 21
TH0 = np.zeros(5)


def rss_mib():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


# ============================================================ constraints
class Bounds:
    """Linear inequalities A theta + b >= 0 (drift pointwise bounds) and box bounds."""

    def __init__(self, design=DESIGN, sc: SiC4H | None = None, levels=(0, 1)):
        n0 = design.n_drift
        self.lo = float(np.log(ND_DRIFT_BOUNDS[0] / n0))
        self.hi = float(np.log(ND_DRIFT_BOUNDS[1] / n0))
        etas = [np.linspace(-1.0, 1.0, 65)]
        if sc is not None:
            for lev in levels:
                d = SiCDMOS(sc, design, lev)
                etas.append(np.unique(np.round(d.dev.eta[d.dev.m_drift], 14)))
        eta = np.unique(np.concatenate(etas))
        self.eta = eta
        P = np.stack([legendre.legval(eta, np.eye(4)[k]) for k in range(4)], axis=1)
        z = np.zeros((len(eta), 1))
        self.A = np.vstack([np.hstack([P, z]), np.hstack([-P, z])])
        self.b = np.r_[-self.lo * np.ones(len(eta)), self.hi * np.ones(len(eta))]
        jlo = float(np.log(ND_JFET_BOUNDS[0] / design.n_jfet))
        jhi = float(np.log(ND_JFET_BOUNDS[1] / design.n_jfet))
        self.box = [(-THETA_BOX, THETA_BOX)] * 4 + [(jlo, jhi)]

    def violation(self, th):
        th = np.asarray(th, float)
        v = [float(np.max(-(self.A @ th + self.b)))]
        for t, (lo, hi) in zip(th, self.box, strict=True):
            v += [lo - t, t - hi]
        return max(v)


# ============================================================ evaluator with cached states
class Problem:
    """Shares one (BV, R_on, gradients) evaluation between objective and constraint."""

    def __init__(self, ev: SiCAdjoint, bv0, r0, ratio, bounds: Bounds):
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
            b = self.ev.bv(th)
            r = self.ev.ron(th)
            self.n_eval += 1
            self.history.append({"theta": th.tolist(), "bv": b["BV"], "R": r["R"]})
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
    res = minimize(
        pr.fun, np.asarray(th_start, float), jac=pr.jac, method="SLSQP", bounds=bounds.box,
        constraints=cons, options={"maxiter": maxiter, "ftol": ftol},
        callback=lambda xk: its.append(np.asarray(xk).tolist()),
    )  # fmt: skip
    return np.asarray(res.x, float), res, pr, its


# ============================================================ gradient checks
def directions():
    out = []
    for s in SEEDS:
        u = np.random.default_rng(s).standard_normal(5)
        out.append(u / np.linalg.norm(u))
    return out


def grad_check(ev: SiCAdjoint, th):
    """Adjoint vs Richardson-extrapolated central FD (eps = 1e-3 and 1e-3/2) along 3 seeded
    directions for BV [V] and R_on [ohm cm^2].  Also recorded: both plain FDs, the FD spread
    |FD(eps) - FD(eps/2)| and a round-off floor estimate (|f| 1e-13 / eps_half)."""
    th = np.asarray(th, float)
    b = ev.bv(th)
    r = ev.ron(th)
    gb, gr = b["grad"], r["grad"]
    r0 = r["R"]
    rows = []
    for s, u in zip(SEEDS, directions(), strict=True):
        rec = {"seed": s}
        for name, fun, g, f0 in (
            ("bv", lambda t: ev.bv(t, grad=False)["BV"], gb, b["BV"]),
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
            "adjoint": rec["ron"]["adjoint"] / r0,
            "richardson": rec["ron"]["richardson"] / r0,
            "abs_diff": rec["ron"]["abs_diff"] / r0,
        }
        rows.append(rec)
    ev.bv(th)
    ev.ron(th)  # re-establish warm states
    return {
        "bv_V": b["BV"], "R_ohm_cm2": r0,
        "grad_bv_norm": float(np.linalg.norm(gb)), "grad_ron_norm": float(np.linalg.norm(gr)),
        "grad_ron_over_R0_norm": float(np.linalg.norm(gr) / r0),
        "grad_bv": gb.tolist(), "grad_ron": gr.tolist(), "directions": rows,
    }  # fmt: skip


# ============================================================ nested uniform family (T6)
def uniform_grid(ev: SiCAdjoint, r0, bounds: Bounds):
    rows = []
    for t0 in np.linspace(bounds.lo, bounds.hi, N_GRID):
        th = np.zeros(5)
        th[0] = t0
        entry = {"theta0": float(t0), "N_D_cm3": float(DESIGN.n_drift * np.exp(t0))}
        try:
            r = ev.ron(th, grad=False)["R"]
            entry["R"] = r
            entry["feasible"] = bool(r <= r0 * (1.0 + FEAS_TOL))
            if entry["feasible"]:
                entry["bv"] = ev.bv(th, grad=False)["BV"]
        except EvalError as exc:
            entry["error"] = str(exc)
            entry["feasible"] = False
        rows.append(entry)
    return rows


# ============================================================ re-verification (T5)
def reverify(sc, designs, levels=(0, 1)):
    """Fresh, cold-started evaluators on mesh h (level 0) and h/2 (level 1)."""
    out = {}
    for lev in levels:
        ev = SiCAdjoint(sc, DESIGN, lev)
        t0 = time.perf_counter()
        rec = {"level": lev, "cells": [int(ev.mesh.nx), int(ev.mesh.ny)], "designs": {}}
        for name, th in designs.items():
            ev.bv_cache = None
            ev.X_on = None  # cold: no warm state from another design
            ev.X_eq = None
            rec["designs"][name] = {
                "bv": ev.bv(th, grad=False, warm=False)["BV"],
                "R": ev.ron(th, grad=False)["R"],
                "cut_x_um": float(sc.length_cm(ev.col_x)) * 1e4,
            }
        rec["seconds"] = time.perf_counter() - t0
        rec["stats"] = dict(ev.stats)
        out[str(lev)] = rec
    return out


def dfdtheta_accuracy(ev: SiCAdjoint, th, step=1e-6):
    """Analytic dF/dtheta vs a residual-only central difference in theta (no re-solve, step 1e-6)
    contracted with random weights, at the converged on-state (all rows) and off-state at the BV
    (Poisson rows only: the continuity rows are Dirichlet there)."""
    th = np.asarray(th, float)
    states = {"on": (ev.ron(th)["X"], False), "off": (ev.bv(th)["state"].X, True)}
    out = {"step": step}
    for name, (X, poisson_only) in states.items():
        errs = []
        for s in SEEDS:
            rng = np.random.default_rng(100 + s)
            w = rng.standard_normal(len(X))
            if poisson_only:
                w[1::3] = 0.0
                w[2::3] = 0.0
            u = directions()[s]
            ev.set_theta(th)
            g = float((ev.dnd_dtheta.T @ ev.dF_dnd(w, X)) @ u)

            def f(t, w=w, X=X):
                ev.set_theta(t)
                F = ev.mod._residual_core(X)
                return w[w != 0] @ F[w != 0]  # (rows with zero weight may be inf in the off-state)

            fd = (f(th + step * u) - f(th - step * u)) / (2 * step)
            errs.append(abs(g - fd) / max(abs(fd), 1e-300))
        out[name] = {"rel_err": errs, "worst": max(errs)}
    ev.set_theta(th)
    return out


def corner_diag(sc, ev: SiCAdjoint, designs):
    """Report-only: the global peak field at each design's BV state and the BV along the column
    through the initial design's global-peak edge (the stage-12 'literal peak' / JFET-corner
    cut), to expose what the fixed planar-cut objective does not see."""
    out = {"designs": {}}
    peaks = {}
    for name, th in designs.items():
        st = ev.bv(th, grad=False)["state"]
        pk = ev.base.peak_field(st.X[0::3])
        peaks[name] = pk
        out["designs"][name] = {
            "E_peak_MV_cm": float(sc.field_phys(pk[0])) * 1e-6,
            "x_um": float(sc.length_cm(pk[1])) * 1e4, "y_um": float(sc.length_cm(pk[2])) * 1e4,
            "dir": pk[3],
        }  # fmt: skip
    x_corner = peaks["initial"][1]
    ec = SiCAdjoint(sc, DESIGN, 0, cut_x=x_corner)
    out["corner_column_x_um"] = float(sc.length_cm(ec.col_x)) * 1e4
    for name, th in designs.items():
        try:
            out["designs"][name]["BV_on_initial_peak_column_V"] = ec.bv(th, grad=False)["BV"]
        except EvalError as exc:
            out["designs"][name]["BV_on_initial_peak_column_V"] = f"error: {exc}"
    out["note"] = "report-only; the optimized objective is the planar cut x = x_bc"
    return out


def live_class_diag(ev: SiCAdjoint, designs):
    """Report-only: BV with the n-type node set (phi_n = V_D) taken from the current doping
    instead of the fixed theta = 0 classification (cold marches, then restored)."""
    out = {}
    saved = (ev.bv_cache, ev.live_class)
    ev.live_class = True
    try:
        for name, th in designs.items():
            ev.bv_cache = None
            try:
                ev.set_theta(th)
                out[name] = {"bv_live": ev.bv(th, grad=False, warm=False)["BV"], "flips": ev.flips}
            except EvalError as exc:
                out[name] = {"error": str(exc), "flips": ev.flips}
    finally:
        ev.live_class = False
        ev.bv_cache = None
    for name, th in designs.items():
        out[name]["bv_fixed"] = ev.bv(th, grad=False, warm=False)["BV"]
        out[name]["rel_diff"] = (
            abs(out[name].get("bv_live", np.nan) - out[name]["bv_fixed"]) / out[name]["bv_fixed"]
        )
    del saved
    return out


# ============================================================ R1 unipolar limit
def unipolar_limit(sc: SiC4H, nd_cm=None):
    nd_cm = np.logspace(15.0, 17.7, 28) if nd_cm is None else nd_cm
    ch = sc.niwa()
    eps_cm = 8.8541878128e-12 * sc.eps_r * 1e-2  # F/cm
    rows = []
    for n in nd_cm:
        e_m, w_d, v_bd = triangular_breakdown(float(sc.dens_s(n)), ch)
        bv = float(sc.volt_phys(v_bd))
        ec = float(sc.field_phys(e_m))
        w_cm = float(sc.length_cm(w_d))
        mu = float(sc.mobility_cm2(n, "n", "par"))
        r = w_cm / (Q * n * mu)
        r_baliga = 4.0 * bv**2 / (eps_cm * mu * ec**3)
        rows.append({"N_D_cm3": float(n), "BV_V": bv, "E_c_V_cm": ec, "W_dep_um": w_cm * 1e4,
                     "mu_n_par_cm2Vs": mu, "R_on_sp_ohm_cm2": r,
                     "R_baliga_ohm_cm2": r_baliga, "rel_formula_diff": abs(r - r_baliga) / r})  # fmt: skip
    return rows


def ideal_r_at(rows, bv):
    b = np.log([r["BV_V"] for r in rows])
    r = np.log([r["R_on_sp_ohm_cm2"] for r in rows])
    order = np.argsort(b)
    return float(np.exp(np.interp(np.log(bv), b[order], r[order])))


# ============================================================ R2 maps
def sensitivity_maps(ev: SiCAdjoint, th):
    """gamma(x, y) = dJ/dN_D,i / vol_s,i per unit area, ZA19/BA26 sign (D = N_A - N_D, so
    gamma_D = -gamma_ND), in units of [BV in V | R in R_0] per (N_0 L_D^2)."""
    b = ev.bv(th)
    r = ev.ron(th)
    vol = ev.dev.reg.vol_s
    ok = vol > 0
    out = {}
    for name, g in (("BV", b["g_nd"]), ("R_on", r["g_nd"] / r["R"])):
        gd = np.zeros_like(g)
        gd[ok] = -g[ok] / vol[ok]
        out[name] = gd
    return out


def plot_maps(ev, gam, sc, out_png, man):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm

    m, reg = ev.mesh, ev.dev.reg
    shape = (m.nx + 1, m.ny + 1)
    um = sc.l_d * 1e4
    X, Y = np.meshgrid(m.x * um, m.y * um, indexing="ij")
    fig, ax = plt.subplots(1, 2, figsize=(12, 5.0), constrained_layout=True)
    labels = {"BV": "gamma_BV  [V per (N_0 L_D^2)], D = N_A - N_D",
              "R_on": "gamma_Ron  [(dR/R_0) per (N_0 L_D^2)], D = N_A - N_D"}  # fmt: skip
    for a, (name, g) in zip(ax, gam.items(), strict=True):
        gm = np.where(reg.semi_node, g, np.nan)
        vmax = float(np.nanpercentile(np.abs(gm), 99.5))  # clip the contact-node outliers
        mp = a.pcolormesh(
            X, Y, np.ma.masked_invalid(gm.reshape(shape)), shading="gouraud", cmap="RdBu_r",
            norm=SymLogNorm(linthresh=1e-3 * vmax, vmin=-vmax, vmax=vmax),
        )  # fmt: skip
        fig.colorbar(mp, ax=a, label=labels[name])
        a.set_ylim(-(sc.length_cm(-ev.prm.y_bottom) * 1e4), 0.0)
        a.set_title(f"{name} sensitivity, initial design")
        a.set_xlabel("x  [um]")
        a.set_ylabel("y  [um]")
    fig.text(0.5, 0.005, man, ha="center", fontsize=7, color="0.4")
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


def plot_pareto(pareto, ideal, rv1, out_png, man):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, a = plt.subplots(1, 1, figsize=(6.4, 4.8), constrained_layout=True)
    bi = np.array([r["BV_V"] for r in ideal])
    ri = np.array([r["R_on_sp_ohm_cm2"] for r in ideal])
    a.loglog(bi, ri, "-", color="0.35", label="unipolar 4H-SiC limit (registry, no fit)")
    pts = [q for q in pareto if q.get("theta") is not None]
    opt = [q for q in pts if q.get("ratio")]
    a.loglog([q["bv"] for q in opt], [q["R"] for q in opt], "o-", color="#D55E00", lw=0.8,
             label="optimized Pareto set (mesh h)")  # fmt: skip
    ini = [q for q in pts if not q.get("ratio")]
    a.loglog(
        [q["bv"] for q in ini], [q["R"] for q in ini], "s", color="#0072B2", label="initial design"
    )
    for q in pts:
        a.annotate("initial" if not q.get("ratio") else f"{q['ratio']:g}",
                   (q["bv"], q["R"]), fontsize=7, xytext=(4, -10), textcoords="offset points")  # fmt: skip
    a.set_xlim(300.0, 1500.0)
    a.set_ylim(2e-5, 1e-3)
    if rv1:
        a.loglog([v["bv"] for v in rv1], [v["R"] for v in rv1], "x", color="#CC79A7",
                 label="same designs re-verified on h/2")  # fmt: skip
    a.set_xlabel("avalanche BV  [V] (fixed vertical cut)")
    a.set_ylabel("R_on,sp  [ohm cm$^2$]  (V_GS = 15 V, V_DS = 0.1 V)")
    a.set_title("Pareto set vs unipolar limit")
    a.legend(fontsize=7)
    a.grid(True, which="both", alpha=0.3)
    fig.text(0.5, 0.005, man, ha="center", fontsize=7, color="0.4")
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


def plot_profile(ev: SiCAdjoint, sc, designs, out_png, man):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    dev, prm = ev.dev, ev.prm
    ys = np.linspace(-(prm.d_body + prm.w_drift), -prm.d_body, 200)
    um = sc.l_d * 1e4
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6), constrained_layout=True)
    cols = ["#0072B2", "#D55E00", "#E69F00", "#CC79A7", "#009E73"]
    for (name, th), col in zip(designs.items(), cols, strict=False):
        nd = dev.drift_profile(th, ys) * sc.n0
        ax[0].plot(nd, -ys * um, color=col, label=f"{name} (JFET {dev.jfet_level(th) * sc.n0:.2e})")
        b = ev.bv(th)
        st = b["state"]
        nodes = ev.cut_nodes
        depth = -ev.mesh.node_xy[nodes, 1]
        e = np.abs(np.diff(st.X[0::3][nodes]) / np.diff(depth)) * sc.e_unit * 1e-6
        ax[1].plot(e, 0.5 * (depth[1:] + depth[:-1]) * um, color=col,
                   label=f"{name}: BV = {b['BV']:.0f} V")  # fmt: skip
    ax[0].set_xscale("log")
    ax[0].set_ylim(ax[0].get_ylim()[::-1])
    ax[0].set_xlabel("drift donor density [cm$^{-3}$]")
    ax[0].set_ylabel("depth below the surface [um]")
    ax[0].set_title("drift profile (body junction at top)")
    ax[0].legend(fontsize=7)
    ax[1].set_ylim(ax[1].get_ylim()[::-1])
    ax[1].set_xlabel("|E| on the cut at BV  [MV/cm]")
    ax[1].set_ylabel("depth [um]")
    ax[1].set_title("field along x = x_bc at the BV state")
    ax[1].legend(fontsize=7)
    fig.text(0.5, 0.005, man, ha="center", fontsize=7, color="0.4")
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


# ============================================================ main
def main(argv=None, maxiter=MAX_ITER, quick=False) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args(argv)
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    man = build_manifest(REPO, sys.argv, list(SEEDS))
    (out / "manifest.json").write_text(json.dumps(man, indent=2) + "\n")
    t_start = time.perf_counter()
    reg = TracedRegistry.from_json()
    sc = SiC4H(reg)
    ev = SiCAdjoint(sc, DESIGN, 0)
    bounds = Bounds(DESIGN, sc)
    T: dict = {}
    res: dict = {}

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
        "mesh_level0": [int(ev.mesh.nx), int(ev.mesh.ny)], "n_nodes_level0": int(ev.mesh.n_nodes),
        "cut_column_x_um": float(sc.length_cm(ev.col_x)) * 1e4, "n_cut_nodes": int(len(ev.cut_nodes)),
        "drift_bounds_cm3": ND_DRIFT_BOUNDS, "jfet_bounds_cm3": ND_JFET_BOUNDS,
        "theta_box": THETA_BOX, "n_profile_constraint_points": int(len(bounds.eta)),
        "pareto_ratios": [1.0, *PARETO_RATIOS], "eps_fd": EPS_FD, "seeds": list(SEEDS),
        "max_iter": MAX_ITER, "V_GS_V": ev.V_GS, "V_DS_V": ev.V_DS,
        "coefficients": "Niwa 2015 (Okuto-Crowell, electron-initiated I_n)",
        "dF_dtheta": "analytic (Poisson via ionized fraction, continuity via mobility factor)",
        "offstate_node_classification": "fixed to the theta = 0 device (n-type nodes carry phi_n = V_D)",
        "x_half_cell_um": ev.x_cm * 1e4,
    }  # fmt: skip

    bv0 = timed("initial_bv", lambda: ev.bv(TH0, grad=False)["BV"])
    r0 = timed("initial_ron", lambda: ev.ron(TH0, grad=False)["R"])
    res["initial"] = {"theta": TH0.tolist(), "bv": bv0, "R": r0}
    gc0 = timed("grad_check_initial", grad_check, ev, TH0)
    res["grad_check_initial"] = gc0
    dump()

    th1, sl, pr, its = timed("slsqp_main", run_slsqp, ev, TH0, bv0, r0, 1.0, bounds, maxiter)
    b1, rr1 = ev.bv(th1, grad=False), ev.ron(th1, grad=False)
    bv1, r1 = b1["BV"], rr1["R"]
    print(f"main: BV {bv0:.3f} -> {bv1:.3f} V, R {r0:.5e} -> {r1:.5e}, {sl.message}", flush=True)
    res["final"] = {"theta": th1.tolist(), "bv": bv1, "R": r1, "R_over_R0_minus_1": r1 / r0 - 1.0,
                    "bound_violation": bounds.violation(th1), "max_violation_during_run": pr.max_violation,
                    "flips": b1["flips"]}  # fmt: skip
    res["slsqp_main"] = {"status": int(sl.status), "message": sl.message, "nit": int(sl.nit),
                         "nfev": int(sl.nfev), "njev": int(sl.njev), "iterates": its,
                         "evaluations": pr.history}  # fmt: skip
    gc1 = timed("grad_check_final", grad_check, ev, th1)
    res["grad_check_final"] = gc1
    dump()

    grid = timed("uniform_grid", uniform_grid, ev, r0, bounds)
    feas = [g["bv"] for g in grid if g.get("feasible") and "bv" in g]
    best_uniform = max(feas) if feas else float("nan")
    res["uniform_grid"], res["best_uniform_bv"] = grid, best_uniform
    dump()

    # ---- R1 Pareto set
    pareto = [
        {
            "ratio": 1.0,
            "theta": th1.tolist(),
            "bv": bv1,
            "R": r1,
            "status": int(sl.status),
            "nit": int(sl.nit),
        }
    ]
    th_prev = th1
    for rt in PARETO_RATIOS:
        th_r, rs, _, _ = timed(
            f"slsqp_ratio_{rt}", run_slsqp, ev, th_prev, bv0, r0, rt, bounds, min(40, maxiter)
        )
        pareto.append({"ratio": rt, "theta": th_r.tolist(), "bv": ev.bv(th_r, grad=False)["BV"],
                       "R": ev.ron(th_r, grad=False)["R"], "status": int(rs.status),
                       "nit": int(rs.nit), "message": rs.message,
                       "bound_violation": bounds.violation(th_r)})  # fmt: skip
        th_prev = th_r
    pareto.append({"ratio": None, "theta": TH0.tolist(), "bv": bv0, "R": r0, "note": "initial"})
    res["R1_pareto"] = pareto
    dump()

    # ---- T5 independent re-verification (fresh cold evaluators) on h and h/2
    designs = {"initial": TH0, "final": th1}
    for q in pareto[1:-1]:
        designs[f"ratio_{q['ratio']:g}"] = np.asarray(q["theta"])
    rv = timed("reverify", reverify, sc, designs)
    res["reverify"] = rv
    dump()

    res["dFdtheta_accuracy"] = timed("dfdtheta_accuracy", dfdtheta_accuracy, ev, th1)
    res["corner_diag"] = timed("corner_diag", corner_diag, sc, ev, {"initial": TH0, "final": th1})
    dump()
    res["live_classification_diag"] = timed(
        "live_class", live_class_diag, ev, {"initial": TH0, "final": th1}
    )
    dump()

    # ---- R1 unipolar limit
    ideal = timed("unipolar_limit", unipolar_limit, sc)
    for q in pareto:
        q["R_ideal_at_BV"] = ideal_r_at(ideal, q["bv"])
        q["R_over_ideal"] = q["R"] / q["R_ideal_at_BV"]
    res["R1_unipolar_limit"] = ideal
    res["R1_pareto"] = pareto
    res["R1_formula"] = (
        "R_ideal = W_dep/(q mu_n,par(N_D) N_D), W_dep = 2 BV/E_m; triangular_breakdown (I_n = 1, Niwa)"
    )
    dump()

    # ---- R2 sensitivity maps (initial design), R3 profile trend
    gam = timed("sensitivity", sensitivity_maps, ev, TH0)
    ev.bv(th1, grad=False)
    st1 = ev.bv(th1, grad=False)["state"]
    psi = st1.X[0::3]
    pk = ev.base.peak_field(psi)
    ys = np.linspace(-(ev.prm.d_body + ev.prm.w_drift), -ev.prm.d_body, 7)
    y_um = (-ys * sc.l_d * 1e4).tolist()
    ycut = ev.dev.mesh.node_xy[ev.cut_nodes, 1]
    near = (ycut < -ev.prm.d_body) & (ycut > -ev.prm.d_body - float(sc.length_s(0.5e-4)))  # 0.5 um
    res["R2"] = {
        "node_x_um": (ev.mesh.node_xy[:, 0] * sc.l_d * 1e4).tolist(),
        "node_y_um": (ev.mesh.node_xy[:, 1] * sc.l_d * 1e4).tolist(),
        "gamma_BV": gam["BV"].tolist(), "gamma_Ron": gam["R_on"].tolist(),
        "units": "BV [V] and R_on/R_0 per (N_0 L_D^2); D = N_A - N_D sign (gamma_D = -gamma_ND)",
        "gamma_BV_cut_mean_below_body": float(np.mean(gam["BV"][ev.cut_nodes][near])),
    }  # fmt: skip
    res["R3"] = {
        "depth_below_body_um": [(-y - ev.prm.d_body) * sc.l_d * 1e4 for y in ys],
        "y_samples_um": y_um,
        "N_D_initial_cm3": (ev.dev.drift_profile(TH0, ys) * sc.n0).tolist(),
        "N_D_final_cm3": (ev.dev.drift_profile(th1, ys) * sc.n0).tolist(),
        "N_D_pareto_cm3": {f"{q['ratio']}": (ev.dev.drift_profile(np.asarray(q['theta']), ys) * sc.n0).tolist() for q in pareto if q.get("ratio")},
        "jfet_initial_cm3": ev.dev.jfet_level(TH0) * sc.n0, "jfet_final_cm3": ev.dev.jfet_level(th1) * sc.n0,
        "jfet_pareto_cm3": {f"{q['ratio']}": ev.dev.jfet_level(np.asarray(q['theta'])) * sc.n0 for q in pareto if q.get("ratio")},
        "global_peak_field_at_final_BV": {"E_MV_cm": float(sc.field_phys(pk[0])) * 1e-6, "x_um": float(sc.length_cm(pk[1])) * 1e4,
                                          "y_um": float(sc.length_cm(pk[2])) * 1e4, "dir": pk[3]},
        "flags_stage12": FLAGS12,
        "note": "qualitative comparison only (BA26 parameters are unpublished); sign convention D = N_A - N_D",
    }  # fmt: skip
    figdir = (REPO / "figures") if out.resolve() == DEFAULT_OUT.resolve() else out / "figures"
    label = f"stage 13 | commit {str(man.get('git_commit'))[:8]} | {man['timestamp_utc'][:19]}Z"
    rv1 = [{"bv": d["bv"], "R": d["R"]} for k, d in rv["1"]["designs"].items()]
    fig_err = None
    try:
        plot_maps(ev, gam, sc, figdir / "stage13_sensitivity_maps.png", label)
        plot_pareto(pareto, ideal, rv1, figdir / "stage13_pareto.png", label)
        dsg = {"initial": TH0, "optimized (R<=R0)": th1}
        for q in pareto[1:-1]:
            dsg[f"R<={q['ratio']:g}R0"] = np.asarray(q["theta"])
        plot_profile(ev, sc, dsg, figdir / "stage13_optimized_profile.png", label)
    except Exception as exc:  # figures are optional
        fig_err = f"{type(exc).__name__}: {exc}"
    res["figures"] = {"dir": str(figdir), "error": fig_err}

    # ---- criteria
    d_ = gc0["directions"]
    t1 = all(r["ron"]["rel_err"] <= 1e-5 for r in d_)
    t2 = all(r["bv"]["rel_err"] <= 1e-4 for r in d_)
    d1 = gc1["directions"]
    tol_bv = 1e-6 * max(1.0, gc1["grad_bv_norm"])
    tol_r = 1e-6 * max(1.0, gc1["grad_ron_over_R0_norm"])
    tol_r_literal = 1e-6 * max(1.0, gc1["grad_ron_norm"])
    t3_bv = all(r["bv"]["abs_diff"] <= tol_bv for r in d1)
    t3_r = all(r["ron_norm"]["abs_diff"] <= tol_r for r in d1)
    t3_r_lit = all(r["ron"]["abs_diff"] <= tol_r_literal for r in d1)
    t3 = bool(t3_bv and t3_r)
    viol = bounds.violation(th1)
    t4 = bool(r1 <= r0 * (1.0 + FEAS_TOL) and viol <= 1e-9)
    l0, l1 = rv["0"]["designs"], rv["1"]["designs"]
    imp0 = l0["final"]["bv"] - l0["initial"]["bv"]
    imp1 = l1["final"]["bv"] - l1["initial"]["bv"]
    r_rel = abs(l1["final"]["R"] - l0["final"]["R"]) / abs(l0["final"]["R"])
    t5 = bool(imp0 > 0 and imp1 > 0 and imp1 >= 0.5 * imp0 and r_rel <= 0.03)
    t6 = bool(bv1 >= best_uniform - 1e-6 * bv1)
    n_failed = ev.stats["failed_evals"] + sum(v["stats"]["failed_evals"] for v in rv.values())
    t7 = bool(sl.status == 0 and sl.nit <= MAX_ITER and n_failed == 0)
    crit = {"T1": t1, "T2": t2, "T3": t3, "T4": t4, "T5": t5, "T6": t6, "T7": t7}
    crit = {k: bool(v) for k, v in crit.items()}
    res["runtime_seconds"] = time.perf_counter() - t_start
    res["peak_rss_MiB"] = rss_mib()
    res["solver_stats_main"] = dict(ev.stats)
    summary = {
        "criteria": crit,
        "gate_passed": all(crit.values()),
        "values": {
            "bv_initial_V": bv0, "bv_final_V": bv1, "R_initial": r0, "R_final": r1,
            "T1_worst_rel_err_Ron": max(r["ron"]["rel_err"] for r in d_),
            "T2_worst_rel_err_BV": max(r["bv"]["rel_err"] for r in d_),
            "T1_T2_fd_spread_worst": [max(r["ron"]["fd_spread"] / abs(r["ron"]["richardson"]) for r in d_),
                                       max(r["bv"]["fd_spread"] / abs(r["bv"]["richardson"]) for r in d_)],
            "T3_worst_abs_BV_V": max(r["bv"]["abs_diff"] for r in d1), "T3_tol_BV": tol_bv,
            "T3_worst_abs_Ron_over_R0": max(r["ron_norm"]["abs_diff"] for r in d1), "T3_tol_Ron_over_R0": tol_r,
            "T3_worst_abs_Ron_ohm_cm2_literal": max(r["ron"]["abs_diff"] for r in d1),
            "T3_tol_Ron_literal": tol_r_literal, "T3_Ron_literal_pass": bool(t3_r_lit),
            "T3_grad_norms": [gc1["grad_bv_norm"], gc1["grad_ron_norm"], gc1["grad_ron_over_R0_norm"]],
            "T3_worst_rel_BV": max(r["bv"]["rel_err"] for r in d1),
            "T3_worst_rel_Ron": max(r["ron"]["rel_err"] for r in d1),
            "T4": {"R_over_R0_minus_1": r1 / r0 - 1.0, "bound_violation": viol},
            "T5": {"improvement_h": imp0, "improvement_h2": imp1, "ratio": imp1 / imp0,
                   "R_final_h": l0["final"]["R"], "R_final_h2": l1["final"]["R"], "R_rel_diff": r_rel,
                   "R_initial_h": l0["initial"]["R"], "R_initial_h2": l1["initial"]["R"],
                   "bv_h": [l0["initial"]["bv"], l0["final"]["bv"]], "bv_h2": [l1["initial"]["bv"], l1["final"]["bv"]]},
            "T6": {"bv_final": bv1, "best_uniform_bv": best_uniform, "margin": bv1 - best_uniform},
            "T7": {"status": int(sl.status), "nit": int(sl.nit), "message": sl.message,
                   "failed_evals": n_failed, "failed_newton_attempts_main": ev.stats["failed_attempts"],
                   "newton_attempts_main": ev.stats["attempts"]},
            "R1": [{k: q.get(k) for k in ("ratio", "bv", "R", "R_ideal_at_BV", "R_over_ideal")} for q in pareto],
            "R1_rel_formula_diff_max": max(r["rel_formula_diff"] for r in ideal),
            "R3": {k: res["R3"][k] for k in ("N_D_initial_cm3", "N_D_final_cm3", "jfet_initial_cm3", "jfet_final_cm3", "global_peak_field_at_final_BV")},
            "flips_max": ev.stats["flips_max"],
            "dFdtheta_vs_residual_fd_worst": [res["dFdtheta_accuracy"]["on"]["worst"], res["dFdtheta_accuracy"]["off"]["worst"]],
            "corner_diag": res["corner_diag"],
            "live_classification_diag": res["live_classification_diag"],
            "warm_fallbacks": ev.stats["warm_fallbacks"],
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
