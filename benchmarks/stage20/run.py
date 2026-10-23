"""Stage-20 benchmark (docs/benchmarks/stage20_body_redesign.md): G1-G3 (gating), R1-R2.

p-body depth family d_body in {0.8, 1.0, 1.2, 1.4, 1.6} um on the stage-13d device (rounded
body r = 0.1 um, theta = 0).  Writes <out>/{results,summary,manifest}.json, <out>/stderr.log
(default results/stage20) and figures stage20_*.png (repo figures/ for the default output, else
<out>/figures).  ``--dev`` refuses an output under results/.

Choices fixed BEFORE any official run (recorded in results.json["config"]):

* Family: ``widegap.bodyredesign.body_design`` changes only d_body; the body y-cell count scales
  with the body thickness (``ny_body = ceil(11 (d_body - d_src)/0.5)``; 11 at 0.8 um) so the
  vertical cell size equals the 0.8 um design's.  Fixed: surface/channel doping, N_A, d_src,
  drift thickness w_drift and doping (the drift starts below the body, so the device deepens by
  the same amount), JFET width and doping (the JFET layer extends to the body depth), oxide.
* Mesh levels h = 0 and h/2 = 1 of the 13d device.
* BV: dense true minimum of the Niwa ionization integral over 49 vertical lines (x = 0..x_max,
  stage-13b ``N_DENSE``), stage-12c/13d off-state, ``dense_bv`` tolerances of stage 13b (0.02 V
  near the minimum, 1/3 V loose at levels 0/1).
* V_pt,es: smallest V_D at which the minimum over columns of psi_surface - min_window psi falls
  below 10 V_T (window: surface to d_body + 0.5 um).  Columns: mesh columns with
  x_src0 <= x <= x_bc whose surface node is net n-type (the n+ source incl. its rounded tail).
  The literal "up to the body-contact edge" would include p-type-surface columns, whose
  psi_surface - min psi is ~0 by construction (no n+ -> p barrier); they are excluded.  March in
  5 V steps (10 V steps above 200 V) from 0, then bisection of the final bracket to 0.5 V and
  linear interpolation of the barrier in the bracket; ceiling 1200 V (reported as ">1200").
* DD confirmation: stage-19 1D DD (``Column1D``, body tied, level 1) at the worst column (the
  column attaining the minimum barrier at V_pt,es on level 1), grid 2 V steps to 150 V, 5 V to
  400 V, 10 V to 1000 V, stopped at the first non-converged target; J_n onsets at 1e-3, 1, 1e3
  A/cm^2.  The DD run is NOT part of G1 (its reach is a numerical limit, reported).
* On state (level 1; E_ox on both levels): R_on,ch (stage-16 ``ChannelRon``), electrostatic V_T
  (N_S* = 1e12 cm^-2 sheet definition at the mid gated p column, V_DS = 0), E_ox = smooth max of
  the oxide field at 400 V (stage 13d).
* G1: every dense BV resolved with no unbracketed line and no failed evaluation, every V_pt,es
  march and on-state solve completed, no "On entry to" line in stderr.log.
* G2: |BV(h/2)/BV(h) - 1| <= 3 % and |V_pt,es(h/2)/V_pt,es(h) - 1| <= 5 % for every design
  (both ">ceiling" counts as equal).
* G3: V_pt,es(h/2) of the 0.8 um design in [78, 141] V.
* R2: smallest d_body with V_pt,es >= 1.2 BV on both meshes.
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

import numpy as np
from scipy.optimize import brentq

from widegap.adjoint_channel import ChannelRon
from widegap.adjoint_sic_corner import DenseCuts, dense_bv
from widegap.adjoint_sic_oxide import OxideAdjoint
from widegap.bodyredesign import OffMarch, body_design, column_index_matrix, source_columns
from widegap.channel_mobility import ChannelMap, NoguchiParams
from widegap.fv2d.mos import continue_bias
from widegap.manifest import build_manifest
from widegap.punchthrough import Column1D, first_crossing
from widegap.sic import SiC4H, TracedRegistry

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[2]
DEFAULT_OUT = REPO / "results" / "stage20"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


M13D = _load("stage13d_run", REPO / "benchmarks/stage13d/run.py")  # fixed name: fork-pool pickling
M13B = M13D.M
DESIGN, TH0 = M13D.DESIGN, M13D.TH0

D_BODY_UM = (0.8, 1.0, 1.2, 1.4, 1.6)
LEVELS = (0, 1)
N_DENSE = M13B.N_DENSE
LOOSE_TOL_V = {0: 1.0, 1: 3.0}
BARRIER_LEVEL_VT = 10.0
NS_STAR = 1e12
D_CH = 0.020
DD_GRID_V = (
    [float(v) for v in range(2, 151, 2)]
    + [float(v) for v in range(155, 401, 5)]
    + [float(v) for v in range(410, 1001, 10)]
)
J_LEVELS = (1e-3, 1.0, 1e3)
TOL = {"G2_bv": 0.03, "G2_vpt": 0.05, "R2_ratio": 1.2}
G3_RANGE_V = (78.0, 141.0)
V_CEIL = 1200.0


def rss_mib():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def make_ev(d_body, level):
    sc = SiC4H(TracedRegistry.from_json())
    return OxideAdjoint(sc, body_design(DESIGN, d_body), level, r_um=M13B.R_UM, cuts_um=(),
                        tau_v=M13B.TAU_V)  # fmt: skip


# ============================================================ jobs (module level: fork pool)
def job_dense(args):
    d_body, level = args
    t0 = time.perf_counter()
    ev = make_ev(d_body, level)
    sc = ev.sc
    x_max_um = float(sc.length_cm(ev.prm.x_max)) * 1e4
    xs = np.linspace(0.0, x_max_um, N_DENSE)
    d = dense_bv(ev, TH0, DenseCuts(ev, xs), tol_v=0.02, loose_tol_v=LOOSE_TOL_V[level])
    bv = np.asarray(d["bv"], dtype=float)
    jm = int(np.nanargmin(bv))
    return {
        "kind": "dense", "d_body": d_body, "level": level, "bv_min": float(bv[jm]),
        "x_argmin_um": float(xs[jm]), "bv": bv.tolist(), "x_um": xs.tolist(),
        "n_unbracketed": int(np.isnan(bv).sum()), "resolved_all": bool(d["resolved_all"]),
        "stats": dict(ev.stats), "seconds": time.perf_counter() - t0,
    }  # fmt: skip


def job_vpt(args):
    """V_pt,es (zero-current barrier) and, on level 1, the DD confirmation at the worst column."""
    d_body, level = args
    t0 = time.perf_counter()
    ev = make_ev(d_body, level)
    sc = ev.sc
    ev.set_theta(TH0)
    cols = source_columns(ev)
    idx = column_index_matrix(ev, cols)
    xs_um = [float(sc.length_cm(ev.dev.mesh.x[i])) * 1e4 for i in cols]
    r = OffMarch(ev, TH0).vpt_es(idx, BARRIER_LEVEL_VT, v_max_v=V_CEIL)
    out = {
        "kind": "vpt", "d_body": d_body, "level": level, "V_pt_es": r["V_pt_es"],
        "V_bracket": r.get("V_bracket"), "columns_um": xs_um, "n_solves": r["n_solves"],
        "curve": [[v, b, xs_um[k]] for v, b, k in r["curve"]],
    }  # fmt: skip
    if r["V_pt_es"] is not None:
        out["worst_x_um"] = xs_um[r["worst_column_row"]]
        out["barriers_by_column_at_bracket_hi"] = r["barriers_at_hi"]
    else:
        out["worst_x_um"] = xs_um[r["curve"][-1][2]]
    if level == 1:
        out["dd"] = dd_confirmation(ev, out["worst_x_um"])
    out["seconds"] = time.perf_counter() - t0
    return out


def dd_confirmation(ev, x_um):
    col = Column1D(ev, x_um, TH0, "tied")
    eq = col.equilibrium()
    if not eq.report.converged:
        return {"error": "equilibrium did not converge"}
    res = col.march(eq.X, 0.0, DD_GRID_V)
    v, j, bar, fail = [], [], [], None
    for vv, sol, _, _ in res:
        if not sol.report.converged:
            fail = float(vv)
            break
        o = col.observe(sol)
        v.append(float(vv))
        j.append(abs(o["Jn_A_cm2"]))
        bar.append(o["barrier_VT"])
    ln = np.log(np.maximum(np.array(j), 1e-300))
    return {
        "x_um_column": col.x_um, "V": v, "Jn_A_cm2": j, "barrier_VT": bar,
        "V_reached": v[-1] if v else 0.0, "first_failed_target_V": fail,
        "onset_V": {f"{lv:g}": first_crossing(v, ln, np.log(lv), True) for lv in J_LEVELS},
        "V_barrier_below_10VT": first_crossing(v, bar, BARRIER_LEVEL_VT, above=False),
    }  # fmt: skip


def vt_es(d_body, level=1):
    """Electrostatic V_T [V]: gate voltage with N_S(mid gated p column) = N_S*, V_DS = 0."""
    ev = make_ev(d_body, level)
    ev.set_theta(TH0)
    sc, base, mod = ev.sc, ev.base, ev.mod
    cm = ChannelMap(sc, mod, ev.nd, ev.na, D_CH, gate_x_max=ev.prm.gate_x1)
    imid = int(np.argmin(np.abs(ev.mesh.x - ev.prm.x_ch_mid)))
    if not cm.col_p[imid]:
        raise RuntimeError("mid column is not p-type")
    eq = ev._eq()
    last = {}

    def ns(vg):
        vgs = float(sc.volt_s(vg))
        if "X" in last:
            s = mod.solve_bc(last["X"], base.bc(vgs, 0.0, 0.0), max_step=20.0)
            if s.report.converged:
                last["X"] = s.X.copy()
                return cm.sheet_density(s.n)[imid]
        out = continue_bias(mod, eq.X, lambda p: base.bc(p[0], 0.0, p[1]), [0.0, 0.0],
                            [[0.0, 0.0], [vgs, 0.0]], max_step=20.0)  # fmt: skip
        s = out[-1][1]
        if not s.report.converged:
            raise RuntimeError("V_T solve did not converge")
        last["X"] = s.X.copy()
        return cm.sheet_density(s.n)[imid]

    f = lambda v: np.log(ns(v) / NS_STAR)  # noqa: E731
    a = b = 10.5
    fa = fb = f(a)
    k = 0
    while fa * fb > 0:
        k += 1
        if fa < 0:
            b += 0.25 * k
            fb = f(b)
        else:
            a -= 0.25 * k
            fa = f(a)
        if k > 60:
            raise RuntimeError("no V_T bracket")
    return float(brentq(f, a, b, xtol=1e-10, rtol=1e-12))


def job_onstate(args):
    (d_body,) = args
    t0 = time.perf_counter()
    out = {"kind": "onstate", "d_body": d_body}
    out["V_T_es_V"] = vt_es(d_body, 1)
    ev = make_ev(d_body, 1)
    cr = ChannelRon(ev, NoguchiParams.from_registry(ev.sc.reg))
    o = cr.ron(TH0, grad=False)
    out["R_on_ch_ohm_cm2"] = float(o["R"])
    out["R_on_picard_iterations"] = int(o["iterations"])
    out["R_on_guard_ok"] = bool(o["guard"]["ok"])
    for lv in LEVELS:
        e = make_ev(d_body, lv)
        out[f"E_ox_400V_MV_cm_L{lv}"] = float(
            e.fields(TH0, grad=False, which=("E_ox",))["E_ox"]["S"]
        )
    out["seconds"] = time.perf_counter() - t0
    return out


def run_job(spec):
    fn = {"dense": job_dense, "vpt": job_vpt, "onstate": job_onstate}[spec[0]]
    try:
        return fn(spec[1:])
    except Exception as exc:
        import traceback

        return {"kind": spec[0], "d_body": spec[1], "level": spec[2] if len(spec) > 2 else None,
                "error": repr(exc), "traceback": traceback.format_exc()}  # fmt: skip


# ============================================================ analysis
def rel(a, b):
    return float(abs(a / b - 1.0))


def evaluate(res, n_blas):
    rows, errors = {}, []
    for d in D_BODY_UM:
        r = {"d_body_um": d}
        for lv in LEVELS:
            jd, jv = res[("dense", d, lv)], res[("vpt", d, lv)]
            for j in (jd, jv):
                if "error" in j:
                    errors.append(f"{j['kind']}|{d}|L{lv}: {j['error']}")
            r[f"BV_L{lv}"] = jd.get("bv_min")
            r[f"x_BV_L{lv}_um"] = jd.get("x_argmin_um")
            r[f"V_pt_es_L{lv}"] = jv.get("V_pt_es", "error")
            r[f"worst_x_um_L{lv}"] = jv.get("worst_x_um")
            r[f"dense_resolved_L{lv}"] = bool(
                jd.get("resolved_all") and jd.get("n_unbracketed", 1) == 0
                and jd.get("stats", {}).get("failed_evals", 1) == 0
            )  # fmt: skip
        on = res[("onstate", d)]
        if "error" in on:
            errors.append(f"onstate|{d}: {on['error']}")
        r.update({k: on.get(k) for k in ("V_T_es_V", "R_on_ch_ohm_cm2", "E_ox_400V_MV_cm_L0",
                                          "E_ox_400V_MV_cm_L1")})  # fmt: skip
        r["R_on_guard_ok"] = on.get("R_on_guard_ok")
        dd = res[("vpt", d, 1)].get("dd", {})
        r["DD"] = {"x_um": dd.get("x_um_column"), "onset_V": dd.get("onset_V"),
                   "V_reached": dd.get("V_reached"), "first_failed_target_V": dd.get("first_failed_target_V"),
                   "V_barrier_below_10VT": dd.get("V_barrier_below_10VT")}  # fmt: skip
        for lv in LEVELS:
            v, bv = r[f"V_pt_es_L{lv}"], r[f"BV_L{lv}"]
            r[f"ratio_L{lv}"] = None if (v in (None, "error") or bv is None) else v / bv
            r[f"ratio_L{lv}_at_least"] = (V_CEIL / bv) if (v is None and bv) else None
        rows[d] = r

    def vpt_ok(a, b):
        if a is None and b is None:
            return True
        if a is None or b is None or "error" in (a, b):
            return False
        return rel(b, a) <= TOL["G2_vpt"]

    g2 = {}
    for d, r in rows.items():
        bv_ok = (
            r["BV_L0"] is not None
            and r["BV_L1"] is not None
            and rel(r["BV_L1"], r["BV_L0"]) <= TOL["G2_bv"]
        )
        g2[d] = {
            "bv_rel_change": None if not (r["BV_L0"] and r["BV_L1"]) else rel(r["BV_L1"], r["BV_L0"]),
            "bv_ok": bool(bv_ok),
            "vpt_rel_change": None
            if (r["V_pt_es_L0"] in (None, "error") or r["V_pt_es_L1"] in (None, "error"))
            else rel(r["V_pt_es_L1"], r["V_pt_es_L0"]),
            "vpt_ok": bool(vpt_ok(r["V_pt_es_L0"], r["V_pt_es_L1"])),
        }  # fmt: skip
    g1_ok = (
        not errors and n_blas == 0
        and all(rows[d][f"dense_resolved_L{lv}"] for d in rows for lv in LEVELS)
        and all(rows[d]["R_on_guard_ok"] for d in rows)
    )  # fmt: skip
    g3v = rows[0.8]["V_pt_es_L1"]
    g3_ok = bool(isinstance(g3v, float) and G3_RANGE_V[0] <= g3v <= G3_RANGE_V[1])
    # R2
    safe = None
    for d in D_BODY_UM:
        r = rows[d]
        ok = True
        for lv in LEVELS:
            v, bv = r[f"V_pt_es_L{lv}"], r[f"BV_L{lv}"]
            if bv is None or v == "error" or not (v is None or v >= TOL["R2_ratio"] * bv):
                ok = False
        if ok:
            safe = d
            break
    r2 = {"rule": "smallest d_body with V_pt,es >= 1.2 BV on both meshes", "safe_design_um": safe}
    if safe is not None:
        b, s = rows[0.8], rows[safe]
        r2.update(
            R_on_ch_penalty=s["R_on_ch_ohm_cm2"] / b["R_on_ch_ohm_cm2"] - 1.0,
            BV_change_L0=s["BV_L0"] / b["BV_L0"] - 1.0,
            BV_change_L1=s["BV_L1"] / b["BV_L1"] - 1.0,
            V_T_change_V=s["V_T_es_V"] - b["V_T_es_V"],
            E_ox_change_L1=s["E_ox_400V_MV_cm_L1"] / b["E_ox_400V_MV_cm_L1"] - 1.0,
        )  # fmt: skip
    crit = {
        "G1": bool(g1_ok),
        "G2": bool(all(g["bv_ok"] and g["vpt_ok"] for g in g2.values())),
        "G3": g3_ok,
    }
    return rows, {"G1_errors": errors, "G2": g2, "G3_value_L1": g3v, "R2": r2}, crit


# ============================================================ figures
def plot_all(figdir: Path, rows: dict, res: dict, label: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figdir.mkdir(parents=True, exist_ok=True)
    ds = list(rows)
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    for lv, mk in ((0, "o--"), (1, "s-")):
        bv = [rows[d][f"BV_L{lv}"] for d in ds]
        vp = [
            rows[d][f"V_pt_es_L{lv}"] if isinstance(rows[d][f"V_pt_es_L{lv}"], float) else V_CEIL
            for d in ds
        ]
        ax[0].plot(ds, bv, mk, label=f"avalanche BV L{lv}", color="C3")
        ax[0].plot(ds, vp, mk, label=f"V_pt,es L{lv}", color="C0")
    ax[0].plot(ds, [1.2 * rows[d]["BV_L1"] for d in ds], ":", color="gray", label="1.2 BV (L1)")
    ax[0].set_xlabel("$d_{body}$ [um]")
    ax[0].set_ylabel("voltage [V]")
    ax[0].legend(fontsize=7)
    ax[1].plot(ds, [rows[d]["R_on_ch_ohm_cm2"] * 1e3 for d in ds], "s-")
    ax[1].set_xlabel("$d_{body}$ [um]")
    ax[1].set_ylabel("R_on,ch [m$\\Omega$ cm$^2$]")
    for d in ds:
        c = res[("vpt", d, 1)]["curve"]
        ax[2].plot([p[0] for p in c], [p[1] for p in c], label=f"{d} um")
    ax[2].axhline(BARRIER_LEVEL_VT, color="k", lw=0.6)
    ax[2].set_xlabel("V_D [V]")
    ax[2].set_ylabel("min barrier [$V_T$]")
    ax[2].legend(fontsize=7)
    fig.suptitle(label, fontsize=8)
    fig.tight_layout()
    fig.savefig(figdir / "stage20_summary.png", dpi=130)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for d in ds:
        dd = res[("vpt", d, 1)].get("dd", {})
        if dd.get("V"):
            ax.semilogy(dd["V"], np.maximum(dd["Jn_A_cm2"], 1e-60), label=f"{d} um")
    ax.set_xlabel("V_D [V]")
    ax.set_ylabel("|J_n| tied-body 1D DD [A/cm$^2$]")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(figdir / "stage20_dd.png", dpi=130)
    plt.close(fig)


# ============================================================ main
def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--workers", type=int, default=min(12, os.cpu_count() or 1))
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
    man = build_manifest(REPO, sys.argv, [])
    man["dev"] = bool(args.dev)
    (out / "manifest.json").write_text(json.dumps(man, indent=2) + "\n")
    t_start = time.perf_counter()
    specs = []
    for d in D_BODY_UM:
        for lv in LEVELS:
            specs.append(("dense", d, lv))
            specs.append(("vpt", d, lv))
        specs.append(("onstate", d))
    # longest first
    specs.sort(
        key=lambda s: (
            0 if s[0] == "dense" and s[2] == 1 else 1 if s[0] == "vpt" and s[2] == 1 else 2
        )
    )
    res = {}
    ctx = mp.get_context("fork")
    with ProcessPoolExecutor(args.workers, mp_context=ctx) as pool:
        futs = {pool.submit(run_job, s): s for s in specs}
        for f in as_completed(futs):
            r = f.result()
            key = (r["kind"], r["d_body"]) + ((r["level"],) if r["kind"] != "onstate" else ())
            res[key] = r
            print(
                f"done {key} {r.get('seconds', float('nan')):.0f}s {r.get('error', '')}", flush=True
            )
    sys.stderr.flush()
    stderr = log_path.read_text(errors="replace") if log_path.exists() else ""
    n_blas = len([ln for ln in stderr.splitlines() if "On entry to" in ln])
    rows, extra, crit = evaluate(res, n_blas)
    runtime = time.perf_counter() - t_start
    results = {
        "config": {
            "d_body_um": list(D_BODY_UM), "levels": list(LEVELS), "n_dense": N_DENSE,
            "barrier_level_VT": BARRIER_LEVEL_VT, "tolerances": TOL, "g3_range_V": list(G3_RANGE_V),
            "v_ceiling_V": V_CEIL, "dd_grid_V": [DD_GRID_V[0], DD_GRID_V[-1]],
            "ny_body": {str(d): body_design(DESIGN, d).mesh_kw["ny_body"] for d in D_BODY_UM},
        },
        "per_design": rows,
        "criteria_detail": extra,
        "on_entry_to_lines": n_blas,
        "jobs": {"|".join(map(str, k)): v for k, v in res.items()},
        "runtime_seconds": runtime,
        "peak_rss_MiB": rss_mib(),
    }  # fmt: skip
    figdir = (REPO / "figures") if out == DEFAULT_OUT.resolve() else out / "figures"
    label = f"stage 20 | commit {str(man.get('git_commit'))[:8]} | {man['timestamp_utc'][:19]}Z"
    try:
        if not any("error" in v for v in res.values()):
            plot_all(figdir, rows, res, label)
    except Exception as exc:
        results["figures_error"] = f"{type(exc).__name__}: {exc}"
    (out / "results.json").write_text(json.dumps(results, indent=1, default=float) + "\n")
    summary = {"criteria": crit, "gate_passed": bool(all(crit.values())), "dev": bool(args.dev),
               "per_design": rows, "criteria_detail": extra, "runtime_seconds": runtime}  # fmt: skip
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    print(json.dumps(crit, indent=2))
    print(json.dumps(extra["R2"], indent=1, default=float))
    print(f"runtime {runtime:.1f} s")


if __name__ == "__main__":
    main()
