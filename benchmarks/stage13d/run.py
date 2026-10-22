"""Stage-13d benchmark (docs/benchmarks/stage13d_oxide_constrained.md): D1-D7 (gating), R1-R2.

Stage 13c (corrected off-state model, rounded device r = 0.1 um, the same K = 10 cuts, tau = 5 V,
49 dense lines, theta bounds, SLSQP) plus an oxide-field constraint:

    maximize  BV_obj(theta)   s.t.  R_on <= R_on,0,   E_ox(theta) <= E_ox(theta_0),   bounds,

E_ox = smooth maximum (log-sum-exp, kappa = 20 per MV/cm) of the field magnitude on the oxide edges
that touch the SiC/SiO2 interface, at the FIXED drain voltage V_blk = 400 V (below the initial true
BV of 517 V; gate and source at 0), with its adjoint gradient (widegap.adjoint_sic_oxide, built on
the corrected off-state electrostatics of widegap.offstate).  Three looser R_on levels
(1.05, 1.10, 1.20 R_on,0), each with the same oxide constraint, give the Pareto set.

Defined local SiC field measures W1 (JFET/body corner), W2 (under the source edge), W3 (JFET-centre
surface): smooth maxima over the semiconductor edges whose midpoints lie in the fixed boxes
WINDOWS_UM of widegap.adjoint_sic_oxide, recorded in that module (and in results.json) BEFORE any
stage-13d run; they are diagnostics, not constraints.

Study choices (nothing tuned after a run): the device, mesh levels, r, cuts, tau, dense set,
bounds, Pareto ratios and SLSQP settings are those of benchmarks/stage13b/run.py (imported through
benchmarks/stage13c/run.py, which replaces its CornerAdjoint by the corrected off-state one).
D3 is evaluated on the quantities in units of their initial values (BV_obj/BV_obj,0, R/R_0,
E_ox/E_ox,0) with tolerance 1e-6 max(1, ||grad||); the 13c-style literal BV [V] and E_ox [MV/cm]
checks are reported alongside.

Writes <out>/{results,summary,comparison,manifest}.json (default results/stage13d) and figures
stage13d_*.png (repo figures/ for the default output, else <out>/figures).
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

from widegap.adjoint_sic import EvalError
from widegap.adjoint_sic_oxide import (
    KAPPA_PER_MV_CM,
    MEASURES,
    V_BLK_V,
    WINDOWS_UM,
    OxideAdjoint,
)
from widegap.manifest import build_manifest
from widegap.sic import SiC4H, TracedRegistry

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[2]
DEFAULT_OUT = REPO / "results" / "stage13d"
R13C = REPO / "results" / "stage13c"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


M13 = _load("stage13c_run", REPO / "benchmarks/stage13c/run.py")  # fixed name: fork-pool pickling
M = M13.M  # stage-13b module with CornerAdjoint -> CorrectedCornerAdjoint
S13 = M.S13
DESIGN, TH0, Bounds = M.DESIGN, M.TH0, M.Bounds

# ---- study parameters ------------------------------------------------------------------------
GUIDELINE_MV_CM = 3.0  # unsourced design guideline, reported only
FEAS_TOL = M.FEAS_TOL
D1_TOL = 0.03
D2_TOL = {"bv": 1e-4, "ron": 1e-5, "E_ox": 1e-4}
D3_TOL = 1e-6
D6_RTOL = 1e-6
PARETO_RATIOS = M.PARETO_RATIOS
EPS_FD = M.EPS_FD
LEVELS_R1 = (0, 1, 2)
SRC_ZONE_UM = ((1.3, 2.9), (-0.6, 0.0))  # n+ source junction zone (built-in field, V_D independent)


def make_ev13d(level, with_cuts=True):
    sc = SiC4H(TracedRegistry.from_json())
    return OxideAdjoint(sc, DESIGN, level, r_um=M.R_UM, cuts_um=M.CUTS_UM if with_cuts else (),
                        tau_v=M.TAU_V)  # fmt: skip


# ============================================================ optimization problem
class Problem:
    """Soft-min BV objective with the R_on and E_ox inequality constraints (shared evaluation)."""

    failure = None

    def __init__(self, ev: OxideAdjoint, bv0, r0, e0, ratio, bounds: Bounds):
        self.ev, self.bv0, self.r0, self.e0, self.ratio, self.bounds = (
            ev,
            bv0,
            r0,
            e0,
            ratio,
            bounds,
        )
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
                f = self.ev.fields(th, which=("E_ox",))["E_ox"]
            except EvalError as exc:
                self.failure = {"theta": th.tolist(), "error": str(exc), "n_eval": self.n_eval}
                raise M.EvalAbort(str(exc)) from exc
            self.n_eval += 1
            self.history.append({"theta": th.tolist(), "bv_obj": b["BV"], "bv_k": b["bv_k"].tolist(),
                                 "R": r["R"], "E_ox": f["S"]})  # fmt: skip
            self.last = (key, (b, r, f))
        return self.last[1]

    def fun(self, th):
        return -self._get(th)[0]["BV"] / self.bv0

    def jac(self, th):
        return -self._get(th)[0]["grad"] / self.bv0

    def c_ron(self, th):
        return self.ratio - self._get(th)[1]["R"] / self.r0

    def c_ron_jac(self, th):
        return -self._get(th)[1]["grad"] / self.r0

    def c_ox(self, th):
        return 1.0 - self._get(th)[2]["S"] / self.e0

    def c_ox_jac(self, th):
        return -self._get(th)[2]["grad"] / self.e0


def run_slsqp(ev, th_start, bv0, r0, e0, ratio, bounds, maxiter=M.MAX_ITER, ftol=1e-10):
    pr = Problem(ev, bv0, r0, e0, ratio, bounds)
    cons = [
        {"type": "ineq", "fun": pr.c_ron, "jac": pr.c_ron_jac},
        {"type": "ineq", "fun": pr.c_ox, "jac": pr.c_ox_jac},
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
    except M.EvalAbort as exc:
        best, bv_best = np.asarray(th_start, float), -np.inf
        for h in pr.history:
            t = np.asarray(h["theta"])
            if (h["R"] <= ratio * r0 * (1.0 + FEAS_TOL) and h["E_ox"] <= e0 * (1.0 + FEAS_TOL)
                    and bounds.violation(t) <= 1e-9 and h["bv_obj"] > bv_best):  # fmt: skip
                best, bv_best = t, h["bv_obj"]
        res = SimpleNamespace(status=-1, message=f"ABORTED (forward solve failed): {exc}",
                              nit=len(its), nfev=pr.n_eval, njev=pr.n_eval, x=best)  # fmt: skip
        return best, res, pr, its


# ============================================================ gradient checks (D2, D3)
def grad_check_fields(ev: OxideAdjoint, th, e0):
    """Adjoint vs Richardson central FD (eps = 1e-3, 5e-4) of the four field measures at V_blk along
    the three seeded directions (E_ox gates; W1-W3 reported)."""
    th = np.asarray(th, float)
    o = ev.fields(th)
    rows = [{"seed": s} for s in M.SEEDS]
    for rec, u in zip(rows, S13.directions(), strict=True):

        def fd(e, u=u):
            fp = ev.fields(th + e * u, grad=False)
            fm = ev.fields(th - e * u, grad=False)
            return {k: (fp[k]["S"] - fm[k]["S"]) / (2.0 * e) for k in MEASURES}

        f1, f2 = fd(EPS_FD), fd(EPS_FD / 2.0)
        for k in MEASURES:
            rich = (4.0 * f2[k] - f1[k]) / 3.0
            a = float(o[k]["grad"] @ u)
            rec[k] = {"adjoint": a, "richardson": rich, "fd_eps": f1[k], "fd_eps_half": f2[k],
                      "abs_diff": abs(a - rich), "rel_err": abs(a - rich) / abs(rich),
                      "fd_spread": abs(f1[k] - f2[k]),
                      "roundoff_floor_est": abs(o[k]["S"]) * 1e-13 / (EPS_FD / 2.0)}  # fmt: skip
    ev.fields(th)
    return {"E_ox_MV_cm": o["E_ox"]["S"], "fields": {k: o[k]["S"] for k in MEASURES},
            "grad_E_norm": float(np.linalg.norm(o["E_ox"]["grad"])),
            "grad_E_over_E0_norm": float(np.linalg.norm(o["E_ox"]["grad"]) / e0),
            "grad_E": o["E_ox"]["grad"].tolist(), "directions": rows}  # fmt: skip


# ============================================================ parallel jobs (module level)
def job_pareto(args):
    """Independent SLSQP run (level 0) from ``th_start`` at R_on <= ratio R_on,0 and E_ox <= E_ox,0."""
    th_start, ratio, bv0, r0, e0, maxiter = args
    t0 = time.perf_counter()
    ev = make_ev13d(0)
    bounds = Bounds(DESIGN, ev.sc)
    th, res, pr, its = run_slsqp(ev, th_start, bv0, r0, e0, ratio, bounds, maxiter=min(40, maxiter))
    b = ev.objective(th, grad=False)
    f = ev.fields(th, grad=False)
    return {
        "ratio": ratio, "theta": th.tolist(), "bv_obj": b["BV"], "bv_k": b["bv_k"].tolist(),
        "R": ev.ron(th, grad=False)["R"], "E_ox": f["E_ox"]["S"],
        "fields": {k: f[k]["S"] for k in MEASURES}, "status": int(res.status),
        "message": res.message, "nit": int(res.nit), "bound_violation": bounds.violation(th),
        "failed_evals": ev.stats["failed_evals"], "seconds": time.perf_counter() - t0,
        "failure": pr.failure, "n_evals": pr.n_eval,
    }  # fmt: skip


def job_fields(args):
    """Field measures (smooth and hard maxima) at drain voltage ``v_volt`` of one design on one mesh."""
    level, theta, tag, v_volt, want_edges = args
    t0 = time.perf_counter()
    ev = make_ev13d(level, with_cuts=False)
    th = np.asarray(theta, float)
    o = ev.fields(th, grad=False, v_volt=v_volt)
    out = {"tag": tag, "level": level, "v_volt": float(v_volt),
           "fields": {k: o[k]["S"] for k in MEASURES},
           "hard_max": {k: o[k]["hard_max"] for k in MEASURES},
           "windows": ev.window_info(), "n_nodes": int(ev.mesh.n_nodes),
           "stats": dict(ev.stats)}  # fmt: skip
    if want_edges:
        ef = ev.edge_fields(o["X"])
        out["edges"] = {k: np.asarray(v).tolist() for k, v in ef.items()}
    out["seconds"] = time.perf_counter() - t0
    return out


# ============================================================ derived analyses
def blind_spot_fields(ef_i, ef_f, windows_um=WINDOWS_UM):
    """Compare per-edge fields of two designs on the same mesh: the largest SiC field outside
    W1-W3 and the n+ source-junction zone, the largest increase over all SiC edges (and outside the
    windows), and the oxide-edge maxima."""
    x, y = np.asarray(ef_i["semi_x_um"]), np.asarray(ef_i["semi_y_um"])
    ei, ef = np.asarray(ef_i["semi_E"]), np.asarray(ef_f["semi_E"])
    inw = np.zeros(len(x), dtype=bool)
    for (x0, x1), (y0, y1) in windows_um.values():
        inw |= (x >= x0 - 1e-9) & (x <= x1 + 1e-9) & (y >= y0 - 1e-9) & (y <= y1 + 1e-9)
    (sx0, sx1), (sy0, sy1) = SRC_ZONE_UM
    zone = (x >= sx0) & (x <= sx1) & (y >= sy0) & (y <= sy1)
    out_ = ~inw & ~zone

    def loc(arr, mask):
        a = np.where(mask, arr, -np.inf)
        k = int(np.argmax(a))
        return {"value": float(a[k]), "x_um": float(x[k]), "y_um": float(y[k]),
                "dir": int(ef_i["semi_dir"][k])}  # fmt: skip

    d = ef - ei
    oi, of = np.asarray(ef_i["ox_E"]), np.asarray(ef_f["ox_E"])
    ox_x = np.asarray(ef_i["ox_x_um"])
    return {
        "max_outside_windows_initial": loc(ei, out_), "max_outside_windows_final": loc(ef, out_),
        "max_all_sic_initial": loc(ei, np.ones(len(x), dtype=bool)),
        "max_all_sic_final": loc(ef, np.ones(len(x), dtype=bool)),
        "max_increase_outside_windows": loc(d, out_), "max_increase_all_sic": loc(d, np.ones(len(x), dtype=bool)),
        "n_edges_outside_windows_with_increase": int(np.sum(d[out_] > 1e-9)),
        "oxide_max_initial": {"value": float(oi.max()), "x_um": float(ox_x[np.argmax(oi)])},
        "oxide_max_final": {"value": float(of.max()), "x_um": float(ox_x[np.argmax(of)])},
        "oxide_max_increase": {"value": float((of - oi).max()), "x_um": float(ox_x[np.argmax(of - oi)])},
    }  # fmt: skip


def plot_all(figdir, res, man_label):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figdir = Path(figdir)
    figdir.mkdir(parents=True, exist_ok=True)
    col = {"initial": "#0072B2", "13d": "#D55E00", "13c": "#009E73"}
    foot = lambda fig: fig.text(0.5, 0.005, man_label, ha="center", fontsize=7, color="0.4")  # noqa: E731

    # 1: per-cut BV
    fig, ax = plt.subplots(figsize=(8.2, 5.0), constrained_layout=True)
    for key, c, name in (
        ("initial_L0", col["initial"], "initial"),
        ("final_L0", col["13d"], "stage 13d"),
        ("c13_L0", col["13c"], "stage-13c optimum"),
    ):
        d = res["dense"].get(key)
        if d is not None:
            ax.plot(
                d["x_um"],
                d["bv"],
                "-",
                color=c,
                label=f"{name}: min {d['bv_min']:.0f} V @ {d['x_argmin_um']:.2f} um",
            )
    ax.set_xlabel("x  [um]")
    ax.set_ylabel("avalanche BV of the vertical cut  [V]")
    ax.set_title("per-cut BV on h (49 dense lines)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    foot(fig)
    fig.savefig(figdir / "stage13d_per_cut_bv.png", dpi=130)
    plt.close(fig)

    # 2: field measures at V_blk
    fig, ax = plt.subplots(figsize=(8.2, 4.6), constrained_layout=True)
    names = ["E_ox", "W1", "W2", "W3"]
    for j, (key, c, name) in enumerate(
        (
            ("initial", col["initial"], "initial"),
            ("final", col["13d"], "stage 13d"),
            ("c13", col["13c"], "stage-13c optimum"),
        )
    ):
        f = res["fields"].get(f"{key}_L0_blk")
        if f is not None:
            ax.bar(
                np.arange(4) + 0.27 * (j - 1),
                [f["fields"][k] for k in names],
                0.25,
                color=c,
                label=name,
            )
    ax.axhline(GUIDELINE_MV_CM, color="0.4", ls="--", lw=1)
    ax.text(
        3.45,
        GUIDELINE_MV_CM,
        "3 MV/cm (unsourced guideline, E_ox only)",
        ha="right",
        va="bottom",
        fontsize=7,
    )
    ax.set_xticks(range(4), ["E_ox (oxide)", "W1 corner", "W2 source edge", "W3 JFET surface"])
    ax.set_ylabel("smooth-max field at V_blk = 400 V  [MV/cm]")
    ax.legend(fontsize=8)
    ax.set_title("oxide and local SiC field measures (h)")
    foot(fig)
    fig.savefig(figdir / "stage13d_fields.png", dpi=130)
    plt.close(fig)

    # 3: oxide field along x
    fig, ax = plt.subplots(figsize=(8.2, 4.6), constrained_layout=True)
    for key, c, name in (
        ("initial", col["initial"], "initial"),
        ("final", col["13d"], "stage 13d"),
        ("c13", col["13c"], "stage-13c optimum"),
    ):
        f = res["fields"].get(f"{key}_L0_blk")
        if f is not None and "edges" in f:
            ax.plot(f["edges"]["ox_x_um"], f["edges"]["ox_E"], "o-", ms=3, color=c, label=name)
    ax.axhline(GUIDELINE_MV_CM, color="0.4", ls="--", lw=1)
    ax.set_xlabel("x  [um]  (JFET centre 0, gate edge 1.6)")
    ax.set_ylabel("oxide field next to the interface  [MV/cm]")
    ax.set_title("gate-oxide field at V_blk = 400 V (h)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    foot(fig)
    fig.savefig(figdir / "stage13d_oxide_field.png", dpi=130)
    plt.close(fig)

    # 4: Pareto
    fig, ax = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
    p13 = res["R2"]["pareto"]
    ax.plot(
        [q["R_over_R0"] for q in p13],
        [q["bv_true_min"] for q in p13],
        "o-",
        color=col["13d"],
        label="stage 13d (E_ox <= E_ox,0)",
    )
    for q in p13:
        ax.annotate(
            f"{q['E_ox_over_E0']:.3f}",
            (q["R_over_R0"], q["bv_true_min"]),
            fontsize=7,
            textcoords="offset points",
            xytext=(3, -9),
        )
    c13 = res["R2"]["pareto_13c"]
    if c13:
        ax.plot(
            [q["R_over_R0"] for q in c13],
            [q["bv_true_min"] for q in c13],
            "s--",
            color=col["13c"],
            label="stage 13c (no oxide constraint)",
        )
    ax.set_xlabel("R_on / R_on,0")
    ax.set_ylabel("dense true-minimum BV on h  [V]")
    ax.set_title("Pareto set (labels: E_ox/E_ox,0 at 400 V)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    foot(fig)
    fig.savefig(figdir / "stage13d_pareto.png", dpi=130)
    plt.close(fig)


# ============================================================ main
def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--maxiter", type=int, default=M.MAX_ITER)
    ap.add_argument("--workers", type=int, default=min(10, os.cpu_count() or 1))
    args = ap.parse_args(argv)
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    man = build_manifest(REPO, sys.argv, list(M.SEEDS))
    (out / "manifest.json").write_text(json.dumps(man, indent=2) + "\n")
    t_start = time.perf_counter()
    ev = make_ev13d(0)
    sc = ev.sc
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
        "rounding_length_um": M.R_UM, "objective_cuts_um_requested": list(M.CUTS_UM),
        "objective_cuts_um_snapped": ev.cut_x_um, "K": len(M.CUTS_UM), "tau_V": M.TAU_V,
        "n_dense_cuts": M.N_DENSE, "mesh_level0": [int(ev.mesh.nx), int(ev.mesh.ny)],
        "n_nodes_level0": int(ev.mesh.n_nodes), "pareto_ratios": [1.0, *PARETO_RATIOS],
        "eps_fd": EPS_FD, "seeds": list(M.SEEDS), "max_iter": args.maxiter, "workers": args.workers,
        "V_blk_V": V_BLK_V, "kappa_per_MV_cm": KAPPA_PER_MV_CM,
        "windows_um": {k: [list(v[0]), list(v[1])] for k, v in WINDOWS_UM.items()},
        "window_edges_level0": ev.window_info(), "src_zone_um_blind_spot_only": [list(SRC_ZONE_UM[0]), list(SRC_ZONE_UM[1])],
        "E_ox_definition": "log-sum-exp (kappa) of |psi_b - psi_a|/h (physical, MV/cm) over the vertical "
                           "oxide edges touching y = 0 (first oxide layer), all mesh columns; gate and source at 0, "
                           "gate Dirichlet psi = V_G - vfb_code (the model's flat-band offset is inside the field)",
        "guideline_MV_cm_unsourced_not_gated": GUIDELINE_MV_CM,
        "offstate_model": "stage 12c (phi_n = V_D w harmonic, phi_p = 0), theta = 0 classification",
    }  # fmt: skip

    # ---- initial design (serial) -----------------------------------------------------------
    o0 = timed("initial_objective", lambda: ev.objective(TH0, grad=False))
    r0 = timed("initial_ron", lambda: ev.ron(TH0, grad=False)["R"])
    f0 = timed("initial_fields", lambda: ev.fields(TH0, grad=False))
    bv0, e0 = o0["BV"], f0["E_ox"]["S"]
    res["initial"] = {"theta": TH0.tolist(), "bv_obj": bv0, "bv_k": o0["bv_k"].tolist(), "R": r0,
                      "E_ox": e0, "fields": {k: f0[k]["S"] for k in MEASURES},
                      "hard_max": {k: f0[k]["hard_max"] for k in MEASURES},
                      "limiting_cut_x_um": ev.cut_x_um[o0["limiting_cut"]]}  # fmt: skip
    gc0 = timed("grad_check_initial", M.grad_check, ev, TH0)
    gf0 = timed("grad_check_fields_initial", grad_check_fields, ev, TH0, e0)
    res["grad_check_initial"], res["grad_check_fields_initial"] = gc0, gf0
    dump()

    th1, sl, pr, its = timed(
        "slsqp_main", run_slsqp, ev, TH0, bv0, r0, e0, 1.0, bounds, args.maxiter
    )
    o1, rr1 = ev.objective(th1, grad=False), ev.ron(th1, grad=False)
    f1 = ev.fields(th1, grad=False)
    bv1, r1, e1 = o1["BV"], rr1["R"], f1["E_ox"]["S"]
    print(
        f"main: BV_obj {bv0:.3f} -> {bv1:.3f} V, R {r0:.5e} -> {r1:.5e}, E_ox {e0:.4f} -> {e1:.4f} "
        f"MV/cm, {sl.message}",
        flush=True,
    )
    res["final"] = {
        "theta": th1.tolist(), "bv_obj": bv1, "bv_k": o1["bv_k"].tolist(), "R": r1,
        "R_over_R0_minus_1": r1 / r0 - 1.0, "E_ox": e1, "E_ox_over_E0_minus_1": e1 / e0 - 1.0,
        "fields": {k: f1[k]["S"] for k in MEASURES},
        "hard_max": {k: f1[k]["hard_max"] for k in MEASURES},
        "limiting_cut_x_um": ev.cut_x_um[o1["limiting_cut"]], "bound_violation": bounds.violation(th1),
        "max_violation_during_run": pr.max_violation, "flips": o1["flips"],
    }  # fmt: skip
    res["slsqp_main"] = {"status": int(sl.status), "message": sl.message, "nit": int(sl.nit),
                         "nfev": int(sl.nfev), "njev": int(sl.njev), "iterates": its,
                         "evaluations": pr.history}  # fmt: skip
    gc1 = timed("grad_check_final", M.grad_check, ev, th1)
    gf1 = timed("grad_check_fields_final", grad_check_fields, ev, th1, e0)
    res["grad_check_final"], res["grad_check_fields_final"] = gc1, gf1
    dump()

    th13 = None
    p13 = R13C / "results.json"
    r13 = json.loads(p13.read_text()) if p13.exists() else None
    if r13 is not None:
        th13 = np.asarray(r13["final"]["theta"], float)

    # ---- phase 2 (parallel) ----------------------------------------------------------------
    jobs_dense = [(0, TH0, "initial", True), (0, th1, "final", True), (1, TH0, "initial", True),
                  (1, th1, "final", True)]  # fmt: skip
    if th13 is not None:
        jobs_dense.append((0, th13, "c13", True))
    jobs_f = []
    for lev in LEVELS_R1:
        jobs_f += [(lev, TH0, f"initial_L{lev}_blk", V_BLK_V, lev == 0),
                   (lev, th1, f"final_L{lev}_blk", V_BLK_V, lev == 0)]  # fmt: skip
    if th13 is not None:
        jobs_f += [(lev, th13, f"c13_L{lev}_blk", V_BLK_V, lev == 0) for lev in (0, 1)]
    t_ph2 = time.perf_counter()
    pareto_res: dict = {}
    dense: dict = {}
    fields: dict = {}
    with ProcessPoolExecutor(args.workers, mp_context=ctx) as pool:
        futs = {}
        for rt in PARETO_RATIOS:
            futs[pool.submit(job_pareto, (th1, rt, bv0, r0, e0, args.maxiter))] = ("pareto", rt)
        for a in jobs_dense:
            futs[pool.submit(M.job_dense, a)] = ("dense", f"{a[2]}_L{a[0]}")
        for a in jobs_f:
            futs[pool.submit(job_fields, a)] = ("fields", a[2])
        for f in as_completed(futs):
            kind, key = futs[f]
            r = f.result()
            {"pareto": pareto_res, "dense": dense, "fields": fields}[kind][key] = r
            print(f"phase2 done: {kind} {key} {r['seconds']:.0f}s", flush=True)
    T["phase2_parallel"] = time.perf_counter() - t_ph2
    res["dense"] = dense
    res["fields"] = fields
    dump()

    # ---- phase 3 (parallel): Pareto designs' dense BV and fields; fields at BV voltages -----------
    v_ref = {lev: dense[f"initial_L{lev}"]["bv_min"] for lev in (0, 1)}
    jobs_d3 = [
        (0, np.asarray(pareto_res[rt]["theta"]), f"ratio_{rt:g}", True) for rt in PARETO_RATIOS
    ]
    jobs_f3 = [(0, TH0, "initial_L0_atBVinit", v_ref[0], True), (0, th1, "final_L0_atBVinit", v_ref[0], True),
               (0, th1, "final_L0_atBVfinal", dense["final_L0"]["bv_min"], True)]  # fmt: skip
    if th13 is not None:
        jobs_f3 += [(0, th13, "c13_L0_atBVinit", v_ref[0], True),
                    (0, th13, "c13_L0_atBVown", dense["c13_L0"]["bv_min"], True)]  # fmt: skip
    jobs_f3 += [(0, np.asarray(pareto_res[rt]["theta"]), f"ratio_{rt:g}_L0_blk", V_BLK_V, True)
                for rt in PARETO_RATIOS]  # fmt: skip
    t_ph3 = time.perf_counter()
    with ProcessPoolExecutor(args.workers, mp_context=ctx) as pool:
        futs = {pool.submit(M.job_dense, a): ("dense", f"{a[2]}_L{a[0]}") for a in jobs_d3}
        futs.update({pool.submit(job_fields, a): ("fields", a[2]) for a in jobs_f3})
        for f in as_completed(futs):
            kind, key = futs[f]
            r = f.result()
            (dense if kind == "dense" else fields)[key] = r
            print(f"phase3 done: {kind} {key} {r['seconds']:.0f}s", flush=True)
    T["phase3_parallel"] = time.perf_counter() - t_ph3
    res["dense"], res["fields"] = dense, fields
    dump()

    # ---- criteria --------------------------------------------------------------------------
    # D1 mesh adequacy of the objective
    b_h, b_h2 = dense["initial_L0"]["bv_min"], dense["initial_L1"]["bv_min"]
    d1_rel = abs(b_h - b_h2) / b_h2
    d1 = bool(d1_rel <= D1_TOL and len(dense["initial_L0"]["x_um"]) >= 40)
    # D2 initial gradients
    d_ = gc0["directions"]
    g_ = gf0["directions"]
    d2_bv = max(r["bv"]["rel_err"] for r in d_)
    d2_r = max(r["ron"]["rel_err"] for r in d_)
    d2_e = max(r["E_ox"]["rel_err"] for r in g_)
    d2 = bool(d2_bv <= D2_TOL["bv"] and d2_r <= D2_TOL["ron"] and d2_e <= D2_TOL["E_ox"])
    # D3 final gradients, units of the initial values
    dd, ge = gc1["directions"], gf1["directions"]
    n_bv, n_r, n_e = (
        gc1["grad_bv_norm"] / bv0,
        gc1["grad_ron_over_R0_norm"],
        gf1["grad_E_over_E0_norm"],
    )
    a_bv = max(r["bv"]["abs_diff"] for r in dd) / bv0
    a_r = max(r["ron_norm"]["abs_diff"] for r in dd)
    a_e = max(r["E_ox"]["abs_diff"] for r in ge) / e0
    t_bv, t_r, t_e = (D3_TOL * max(1.0, n) for n in (n_bv, n_r, n_e))
    d3 = bool(a_bv <= t_bv and a_r <= t_r and a_e <= t_e)
    lit_bv_tol = D3_TOL * max(1.0, gc1["grad_bv_norm"])
    lit_e_tol = D3_TOL * max(1.0, gf1["grad_E_norm"])
    lit_bv = max(r["bv"]["abs_diff"] for r in dd)
    lit_e = max(r["E_ox"]["abs_diff"] for r in ge)
    # D4 feasibility
    viol = bounds.violation(th1)
    d4 = bool(r1 <= r0 * (1.0 + FEAS_TOL) and e1 <= e0 * (1.0 + FEAS_TOL) and viol <= 1e-9)
    # D5 dense true-minimum BV
    ini = {lev: dense[f"initial_L{lev}"] for lev in (0, 1)}
    fin = {lev: dense[f"final_L{lev}"] for lev in (0, 1)}
    imp = {lev: fin[lev]["bv_min"] - ini[lev]["bv_min"] for lev in (0, 1)}
    d5 = bool(imp[0] > 0 and imp[1] > 0 and imp[1] >= 0.5 * imp[0])
    # D6 no local measure increases at V_blk
    d6_rec = {}
    for lev in (0, 1):
        fi, ff = fields[f"initial_L{lev}_blk"]["fields"], fields[f"final_L{lev}_blk"]["fields"]
        d6_rec[str(lev)] = {k: {"initial": fi[k], "final": ff[k], "rel_change": ff[k] / fi[k] - 1.0,
                                "ok": bool(ff[k] <= fi[k] * (1.0 + D6_RTOL))} for k in MEASURES}  # fmt: skip
    d6 = bool(all(v["ok"] for lv in d6_rec.values() for v in lv.values()))
    # D7
    n_failed = (
        ev.stats["failed_evals"]
        + sum(d["stats"]["failed_evals"] for d in dense.values())
        + sum(d["stats"]["failed_evals"] for d in fields.values())
    )
    d7 = bool(sl.status == 0 and sl.nit <= args.maxiter and n_failed == 0
              and all(d["resolved_all"] for d in dense.values()))  # fmt: skip
    crit = {
        k: bool(v)
        for k, v in {"D1": d1, "D2": d2, "D3": d3, "D4": d4, "D5": d5, "D6": d6, "D7": d7}.items()
    }

    # ---- R1: fields on h, h/2, h/4 ---------------------------------------------------------
    r1_rec = {}
    for tag in ("initial", "final"):
        r1_rec[tag] = {}
        for k in MEASURES:
            vals = [fields[f"{tag}_L{lev}_blk"]["fields"][k] for lev in LEVELS_R1]
            r1_rec[tag][k] = {"values_MV_cm_h_h2_h4": vals, **M.u1_richardson(vals)}
    # ---- R2 --------------------------------------------------------------------------------
    ideal = M.unipolar_limit(sc)
    par = [{"ratio": 1.0, "theta": th1.tolist(), "bv_obj": bv1, "R": r1, "E_ox": e1,
            "fields": res["final"]["fields"], "status": int(sl.status), "nit": int(sl.nit),
            "bv_true_min": fin[0]["bv_min"], "x_argmin_um": fin[0]["x_argmin_um"]}]  # fmt: skip
    for rt in PARETO_RATIOS:
        q = dict(pareto_res[rt])
        d = dense[f"ratio_{rt:g}_L0"]
        q["bv_true_min"], q["x_argmin_um"] = d["bv_min"], d["x_argmin_um"]
        par.append(q)
    for q in par:
        q["R_over_R0"] = q["R"] / r0
        q["E_ox_over_E0"] = q["E_ox"] / e0
        q["R_ideal_at_bv_true_min"] = M.ideal_r_at(ideal, q["bv_true_min"])
        q["R_over_ideal_bv_true_min"] = q["R"] / q["R_ideal_at_bv_true_min"]
    par.insert(0, {"ratio": None, "note": "initial", "theta": TH0.tolist(), "bv_obj": bv0, "R": r0,
                   "E_ox": e0, "R_over_R0": 1.0, "E_ox_over_E0": 1.0, "fields": res["initial"]["fields"],
                   "bv_true_min": ini[0]["bv_min"], "x_argmin_um": ini[0]["x_argmin_um"]})  # fmt: skip
    par13 = []
    if r13 is not None:
        for q in r13["R1_pareto"]:
            if q.get("ratio"):
                par13.append(
                    {
                        "ratio": q["ratio"],
                        "R_over_R0": q["R"] / r13["initial"]["R"],
                        "bv_true_min": q["bv_true_min"],
                        "bv_obj": q["bv_obj"],
                    }
                )
    f_ini = fields["initial_L0_blk"]["fields"]
    ab = lambda k: fields[k]["fields"]["E_ox"]  # noqa: E731
    r2 = {
        "E_ox_MV_cm": {"initial_h": ab("initial_L0_blk"), "final_h": ab("final_L0_blk"),
                       "initial_h2": ab("initial_L1_blk"), "final_h2": ab("final_L1_blk"),
                       "initial_h4": ab("initial_L2_blk"), "guideline": GUIDELINE_MV_CM,
                       "initial_exceeds_guideline": bool(ab("initial_L0_blk") > GUIDELINE_MV_CM),
                       "final_exceeds_guideline": bool(ab("final_L0_blk") > GUIDELINE_MV_CM),
                       "initial_hard_max_h": fields["initial_L0_blk"]["hard_max"]["E_ox"],
                       "final_hard_max_h": fields["final_L0_blk"]["hard_max"]["E_ox"],
                       "E_ox_at_initial_true_BV_initial": ab("initial_L0_atBVinit"),
                       "E_ox_at_initial_true_BV_final": ab("final_L0_atBVinit"),
                       "E_ox_final_at_own_true_BV": ab("final_L0_atBVfinal")},
        "BV_true_min_V": {"initial_h_h2": [ini[0]["bv_min"], ini[1]["bv_min"]],
                          "final_h_h2": [fin[0]["bv_min"], fin[1]["bv_min"]]},
        "R_on": {"initial": r0, "final": r1}, "windows_V_blk_h": {"initial": f_ini, "final": fields["final_L0_blk"]["fields"]},
        "pareto": par, "pareto_13c": par13,
    }  # fmt: skip
    if th13 is not None:
        r2["stage13c_optimum"] = {
            "theta": th13.tolist(), "E_ox_MV_cm_h": ab("c13_L0_blk"), "E_ox_MV_cm_h2": ab("c13_L1_blk"),
            "E_ox_over_E0_h": ab("c13_L0_blk") / ab("initial_L0_blk"),
            "fields_h": fields["c13_L0_blk"]["fields"], "bv_true_min_h": dense["c13_L0"]["bv_min"],
            "R_over_R0": dense["c13_L0"]["R"] / r0, "x_argmin_um": dense["c13_L0"]["x_argmin_um"],
            "bv_true_min_13c_reported": [r13["dense"]["final_L0"]["bv_min"], r13["dense"]["final_L1"]["bv_min"]],
        }  # fmt: skip
    res["R1"] = r1_rec
    res["R2"] = r2

    # ---- blind-spot analysis (report) --------------------------------------------------------
    blind = {
        "final_vs_initial_h": blind_spot_fields(
            fields["initial_L0_blk"]["edges"], fields["final_L0_blk"]["edges"]
        )
    }
    if th13 is not None:
        blind["c13_vs_initial_h"] = blind_spot_fields(
            fields["initial_L0_blk"]["edges"], fields["c13_L0_blk"]["edges"]
        )
    for rt in PARETO_RATIOS:
        blind[f"ratio_{rt:g}_vs_initial_h"] = blind_spot_fields(
            fields["initial_L0_blk"]["edges"], fields[f"ratio_{rt:g}_L0_blk"]["edges"]
        )
    cut_x = np.asarray(M.CUTS_UM)
    for tag in ("initial", "final"):
        for lev in (0, 1):
            d = dense[f"{tag}_L{lev}"]
            blind[f"dense_BV_{tag}_L{lev}"] = {
                "dense_min_V": d["bv_min"], "x_argmin_um": d["x_argmin_um"],
                "distance_to_nearest_objective_cut_um": float(np.min(np.abs(cut_x - d["x_argmin_um"])))}  # fmt: skip
    for lev in (0, 1):
        bi, bf = np.asarray(ini[lev]["bv"]), np.asarray(fin[lev]["bv"])
        both = ~np.isnan(bi) & ~np.isnan(bf)
        blind[f"L{lev}_x_where_final_BV_lower_than_initial_um"] = [
            float(x) for x, a, b in zip(fin[lev]["x_um"], bi, bf, strict=True)
            if not (np.isnan(a) or np.isnan(b)) and b < a]  # fmt: skip
        blind[f"L{lev}_worst_cut_BV_drop_V"] = float(np.max(bi[both] - bf[both]))
    res["blind_spot_analysis"] = blind

    # ---- figures ----------------------------------------------------------------------------
    figdir = (REPO / "figures") if out.resolve() == DEFAULT_OUT.resolve() else out / "figures"
    label = f"stage 13d | commit {str(man.get('git_commit'))[:8]} | {man['timestamp_utc'][:19]}Z"
    fig_err = None
    try:
        plot_all(figdir, res, label)
    except Exception as exc:  # figures are optional
        fig_err = f"{type(exc).__name__}: {exc}"
    res["figures"] = {"dir": str(figdir), "error": fig_err}

    res["runtime_seconds"] = time.perf_counter() - t_start
    res["solver_stats_main"] = dict(ev.stats)
    pareto_failures = {str(k): v["failure"] for k, v in pareto_res.items() if v["failure"]}
    summary = {
        "criteria": crit, "gate_passed": all(crit.values()),
        "values": {
            "bv_obj_initial_V": bv0, "bv_obj_final_V": bv1, "R_initial": r0, "R_final": r1,
            "E_ox_initial_MV_cm": e0, "E_ox_final_MV_cm": e1,
            "D1": {"bv_min_h_h2": [b_h, b_h2], "rel_diff": d1_rel, "tol": D1_TOL,
                   "n_dense": len(dense["initial_L0"]["x_um"])},
            "D2": {"worst_rel_BV_obj": d2_bv, "worst_rel_Ron": d2_r, "worst_rel_E_ox": d2_e,
                   "worst_rel_W": {k: max(r[k]["rel_err"] for r in g_) for k in ("W1", "W2", "W3")},
                   "fd_spread_worst_E_ox": max(r["E_ox"]["fd_spread"] / abs(r["E_ox"]["richardson"]) for r in g_),
                   "roundoff_floor_E_ox_abs": max(r["E_ox"]["roundoff_floor_est"] for r in g_)},
            "D3": {"worst_abs_BV_over_BV0": a_bv, "tol_BV": t_bv, "worst_abs_Ron_over_R0": a_r, "tol_Ron": t_r,
                   "worst_abs_E_over_E0": a_e, "tol_E": t_e, "grad_norms_normalized": [n_bv, n_r, n_e],
                   "literal_BV_V": {"worst_abs": lit_bv, "tol": lit_bv_tol, "pass": bool(lit_bv <= lit_bv_tol)},
                   "literal_E_ox_MV_cm": {"worst_abs": lit_e, "tol": lit_e_tol, "pass": bool(lit_e <= lit_e_tol)},
                   "W_worst_abs_MV_cm": {k: max(r[k]["abs_diff"] for r in ge) for k in ("W1", "W2", "W3")}},
            "D4": {"R_over_R0_minus_1": r1 / r0 - 1.0, "E_over_E0_minus_1": e1 / e0 - 1.0, "bound_violation": viol},
            "D5": {"min_bv_initial_h_h2": [ini[0]["bv_min"], ini[1]["bv_min"]],
                   "min_bv_final_h_h2": [fin[0]["bv_min"], fin[1]["bv_min"]],
                   "improvement_h": imp[0], "improvement_h2": imp[1],
                   "ratio": imp[1] / imp[0] if imp[0] != 0 else None,
                   "x_argmin_initial_um": [ini[0]["x_argmin_um"], ini[1]["x_argmin_um"]],
                   "x_argmin_final_um": [fin[0]["x_argmin_um"], fin[1]["x_argmin_um"]]},
            "D6": d6_rec,
            "D7": {"status": int(sl.status), "nit": int(sl.nit), "message": sl.message, "failed_evals": n_failed,
                   "pareto_failures": pareto_failures,
                   "pareto_status": {str(k): v["status"] for k, v in pareto_res.items()},
                   "dense_resolved_all": {k: v["resolved_all"] for k, v in dense.items()}},
            "R1": r1_rec, "R2": r2, "blind_spot": blind, "flips_max": ev.stats["flips_max"],
        },
        "seconds": T, "runtime_seconds": res["runtime_seconds"],
    }  # fmt: skip
    dump()
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    cmp_ = {"stage13c_optimum": r2.get("stage13c_optimum"), "pareto_13d": par, "pareto_13c": par13,
            "criteria_13c": json.loads((R13C / "summary.json").read_text())["criteria"]
            if (R13C / "summary.json").exists() else None}  # fmt: skip
    (out / "comparison.json").write_text(json.dumps(cmp_, indent=2, default=float) + "\n")
    print(json.dumps(summary["criteria"], indent=2))
    print(f"runtime {res['runtime_seconds']:.1f} s")


if __name__ == "__main__":
    main()
