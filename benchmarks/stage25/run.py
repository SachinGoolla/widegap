"""Stage-25 benchmark (docs/benchmarks/stage25_generality_hatakeyama.md): H1-H3 (gating), R1, P25a-c.

The stage-23 artifact budget re-evaluated with the Hatakeyama 2004 impact-ionization set
(``SiC4H.hatakeyama()``, the stage-12 set: a, b of the E || c electron/hole fits) instead of the
Niwa 2015 set (``SiC4H.niwa()``, the default).  No re-optimization; no library change: the
coefficient set is the ``ev.ch`` attribute that ``SiCAdjoint`` already exposes (default
``sc.niwa()``), set right after construction, before any solve.  It is used by the planar-cut BV,
the adjoint root find and ``DenseCuts``/``dense_bv`` alike.

Choices fixed BEFORE any official run (recorded in results.json["config"]):

* AB1': the stage-23 AB1 jobs (sharp stage-13 device, original off-state, mesh h, planar cut and 49 + 1
  dense lines, same tolerances) with ``ev.ch`` = Hatakeyama.
* AB3': the stage-23 AB3 dense job (0.8 um geometry via ``stage20.make_ev(0.8, 0)``, 49 lines, 12c
  model, h) with Hatakeyama; blocking = min(BV_Hat, V_pt), V_pt = the stage-23 DD onset
  (results/stage23/results.json AB3/<design>/dd_onset_V, unchanged).
* AB6: the stage-21 1.4 um geometry via ``stage20.make_ev(1.4, level)`` (``body_design``, the same
  geometry stage 21 installs; the stage-21 diagnostic windows do not enter the BV), the 49 + 10
  objective-cut lines of the stage-13b ``job_dense`` and the recorded stage-21 initial (theta = 0)
  and final theta; dense true-min BV gain on h (level 0) under Hatakeyama vs the recorded Niwa gain.
  Level 1 is computed under Niwa only (H1).
* H1: with Niwa selected the stage-23 AB1 (planar and dense, both designs) and AB3 dense values and
  the stage-21 dense BVs (levels 0 and 1, initial and final) are recomputed and compared with the
  records (relative <= 1e-6).
* H2: no job error, every dense BV resolved with failed_evals = 0 (lines that never bracket lie above
  2x the first-crossing voltage; counted and flagged), no "On entry to" line.
* H3: Hatakeyama dense minimum of the 0.8 um initial design, level 0 vs level 1, relative <= 3 %.
* Predictions P25a (AB1' artifact fraction >= 0.5), P25b (every AB3' design has V_pt < BV_Hat and
  artifact fraction >= 0.7), P25c (AB6 Hatakeyama gain > 0 and within +-30 % of the Niwa gain) are
  reported confirmed/refuted, never changed.
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

from widegap.adjoint_sic import SiCAdjoint
from widegap.adjoint_sic_corner import DenseCuts, dense_bv
from widegap.artifactbudget import Cites, blocking_voltage, entry
from widegap.manifest import build_manifest
from widegap.sic import SiC4H, TracedRegistry

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[2]
DEFAULT_OUT = REPO / "results" / "stage25"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


M23 = _load("stage23_run", REPO / "benchmarks/stage23/run.py")  # fixed names: fork pickling
M20, M13 = M23.M20, M23.M13
RES = M23.RES
DESIGNS = M23.DESIGNS
SETS = ("niwa", "hatakeyama")
TOL = {"H1_rel": 1e-6, "H3_rel": 0.03, "P25a": 0.5, "P25b": 0.7, "P25c": 0.30}
N_DENSE = M20.N_DENSE
CUTS_UM = M20.M13B.CUTS_UM
S21_THETA = ("stage21", ("final", "theta"))


def set_ch(ev, cs):
    ev.ch = ev.sc.hatakeyama() if cs == "hatakeyama" else ev.sc.niwa()
    return ev


# ============================================================ jobs (module level: fork pool)
def job_ab1(args):
    cs, name, theta = args
    t0 = time.perf_counter()
    sc = SiC4H(TracedRegistry.from_json())
    ev = set_ch(SiCAdjoint(sc, M13.DESIGN, 0), cs)
    th = np.asarray(theta, float)
    ev.bv_cache = None
    ev.X_on = None
    ev.X_eq = None
    planar = float(ev.bv(th, grad=False, warm=False)["BV"])
    x_max_um = float(sc.length_cm(ev.prm.x_max)) * 1e4
    x_bc_um = float(sc.length_cm(ev.prm.x_bc)) * 1e4
    xs = np.linspace(0.0, x_max_um, N_DENSE)
    d = dense_bv(ev, th, DenseCuts(ev, np.r_[xs, x_bc_um]), tol_v=0.02,
                 loose_tol_v=M20.LOOSE_TOL_V[0])  # fmt: skip
    bv = np.asarray(d["bv"], float)
    bd = bv[:N_DENSE]
    jm = int(np.nanargmin(bd))
    return {
        "kind": "ab1", "cs": cs, "name": name, "planar_bv": planar, "dense_bv_min": float(bd[jm]),
        "x_argmin_um": float(xs[jm]), "dense_bv_at_x_bc": float(bv[-1]),
        "n_unbracketed": int(np.isnan(bd).sum()), "resolved_all": bool(d["resolved_all"]),
        "stats": dict(ev.stats), "seconds": time.perf_counter() - t0,
    }  # fmt: skip


def job_dense(args):
    """Dense true-min BV: geometry '0.8' (49 lines, stage 23) or '1.4' (49 + 10 cut lines, stage 21)."""
    cs, geom, name, theta, level = args
    t0 = time.perf_counter()
    ev = set_ch(M20.make_ev(float(geom), level), cs)
    sc = ev.sc
    xs = np.linspace(0.0, float(sc.length_cm(ev.prm.x_max)) * 1e4, N_DENSE)
    xs_all = np.r_[xs, np.asarray(CUTS_UM)] if geom == "1.4" else xs
    d = dense_bv(ev, np.asarray(theta, float), DenseCuts(ev, xs_all), tol_v=0.02,
                 loose_tol_v=M20.LOOSE_TOL_V[level])  # fmt: skip
    bv = np.asarray(d["bv"], float)[:N_DENSE]
    jm = int(np.nanargmin(bv))
    return {
        "kind": "dense", "cs": cs, "geom": geom, "name": name, "level": level,
        "bv_min": float(bv[jm]), "x_argmin_um": float(xs[jm]),
        "n_unbracketed": int(np.isnan(bv).sum()), "resolved_all": bool(d["resolved_all"]),
        "stats": dict(ev.stats), "seconds": time.perf_counter() - t0,
    }  # fmt: skip


def run_job(spec):
    fn = {"ab1": job_ab1, "dense": job_dense}[spec[0]]
    try:
        return fn(spec[1:])
    except Exception as exc:
        import traceback

        return {"kind": spec[0], "cs": spec[1], "error": repr(exc),
                "traceback": traceback.format_exc()}  # fmt: skip


# ============================================================ analysis
def rel(a, b):
    return float(abs(a / b - 1.0))


def evaluate(res, c: Cites, n_blas, thetas, th21):
    errors = [f"{k}: {v['error']}" for k, v in res.items() if v.get("error")]
    g = lambda *k: res.get(k)  # noqa: E731
    table = []
    # ---------------- H1 (Niwa reproduction)
    h1 = {}
    s23 = RES.format("stage23")
    for n in ("initial", "13"):
        r = g("ab1", "niwa", n)
        rec_p = c.ref(s23, "T1", "recomputed", "bv_initial" if n == "initial" else "bv_final")
        rec_d = c.ref(s23, "detail", "AB1", n, "dense_bv_min")
        if r and not r.get("error"):
            h1[f"AB1|{n}|planar"] = rel(r["planar_bv"], rec_p["value"])
            h1[f"AB1|{n}|dense"] = rel(r["dense_bv_min"], rec_d["value"])
    for n in DESIGNS:
        r = g("dense", "niwa", n, "0.8", 0)
        rec = c.ref(s23, "AB3", n, "bv_recomputed")
        if r and not r.get("error"):
            h1[f"AB3|{n}|dense"] = rel(r["bv_min"], rec["value"])
    for lv in (0, 1):
        for n, tag in (("initial", "initial"), ("final", "final")):
            r = g("dense", "niwa", n, "1.4", lv)
            rec = c.ref(RES.format("stage21"), "dense", f"{tag}_L{lv}", "bv_min")
            if r and not r.get("error"):
                h1[f"stage21|{n}|L{lv}"] = rel(r["bv_min"], rec["value"])
    h1_ok = bool(h1) and len(h1) == 4 + 5 + 4 and max(h1.values()) <= TOL["H1_rel"]
    # ---------------- AB1'
    a = {cs: {n: g("ab1", cs, n) for n in ("initial", "13")} for cs in SETS}
    ab1_ok = all(a[cs][n] and not a[cs][n].get("error") for cs in SETS for n in a[cs])
    pred = {}
    if ab1_ok:
        for cs in SETS:
            gc = a[cs]["13"]["planar_bv"] - a[cs]["initial"]["planar_bv"]
            gr = a[cs]["13"]["dense_bv_min"] - a[cs]["initial"]["dense_bv_min"]
            table.append(
                entry(
                    f"AB1[{cs}]",
                    f"stage 13 planar -> dense, old off-state, sharp device, h, {cs}",
                    "BV gain initial -> final",
                    gc,
                    gr,
                    "V",
                    {
                        "claimed": "recomputed planar-cut BV",
                        "real": "recomputed dense true-min BV (49 lines)",
                    },
                    f"planar {a[cs]['initial']['planar_bv']:.2f} -> {a[cs]['13']['planar_bv']:.2f}; dense "
                    f"{a[cs]['initial']['dense_bv_min']:.2f} -> {a[cs]['13']['dense_bv_min']:.2f} V",
                    claimed_endpoints=[a[cs]["initial"]["planar_bv"], a[cs]["13"]["planar_bv"]],
                    real_endpoints=[a[cs]["initial"]["dense_bv_min"], a[cs]["13"]["dense_bv_min"]],
                    coefficient_set=cs,
                )
            )
        af = [t for t in table if t["id"] == "AB1[hatakeyama]"][0]["artifact_fraction"]
        pred["P25a"] = {"artifact_fraction": af, "threshold": TOL["P25a"],
                        "confirmed": bool(af is not None and af >= TOL["P25a"])}  # fmt: skip
    # ---------------- AB3'
    ab3, p25b_ok, p25b_rows = {}, True, {}
    for n in DESIGNS:
        rh = g("dense", "hatakeyama", n, "0.8", 0)
        if not rh or rh.get("error"):
            p25b_ok = False
            continue
        on = c.ref(s23, "AB3", n, "dd_onset_V")
        rc = c.ref(s23, "AB3", n, "dd_V_reached")
        bn = c.ref(s23, "AB3", n, "bv_recomputed")
        blk = blocking_voltage(rh["bv_min"], on["value"], rc["value"])
        ab3[n] = {"bv_hat": rh["bv_min"], "bv_niwa": bn["value"], "V_pt": on["value"], **blk}
        row = entry(
            f"AB3[{n}][hatakeyama]", f"{M23.DESIGN_LABEL[n]}, 0.8 um, 12c model, h, Hatakeyama",
            "blocking voltage min(BV, V_pt) vs avalanche BV", rh["bv_min"], blk["blocking_V"], "V",
            {"claimed": "recomputed dense true-min BV (Hatakeyama)",
             "real": [on, rc]}, f"V_pt (DD onset J_n = 1 A/cm2, stage 23) = {on['value']:.2f} V; "
            f"limited by {blk['limited_by']}", coefficient_set="hatakeyama",
        )  # fmt: skip
        table.append(row)
        p25b_rows[n] = {"V_pt": on["value"], "BV_hat": rh["bv_min"],
                        "artifact_fraction": row["artifact_fraction"],
                        "punch_through_limited": bool(on["value"] < rh["bv_min"]),
                        "ok": bool(on["value"] < rh["bv_min"] and row["artifact_fraction"] >= TOL["P25b"])}  # fmt: skip
        p25b_ok = p25b_ok and p25b_rows[n]["ok"]
    pred["P25b"] = {"by_design": p25b_rows, "threshold": TOL["P25b"],
                    "confirmed": bool(p25b_ok and len(p25b_rows) == len(DESIGNS))}  # fmt: skip
    # ---------------- AB6
    d6 = {cs: {n: g("dense", cs, n, "1.4", 0) for n in ("initial", "final")} for cs in SETS}
    gains = {}
    s21 = RES.format("stage21")
    gn = (
        c.ref(s21, "dense", "final_L0", "bv_min")["value"]
        - c.ref(s21, "dense", "initial_L0", "bv_min")["value"]
    )
    if all(d6[cs][n] and not d6[cs][n].get("error") for cs in SETS for n in d6[cs]):
        gains = {cs: d6[cs]["final"]["bv_min"] - d6[cs]["initial"]["bv_min"] for cs in SETS}
        for cs in SETS:
            table.append(
                entry(
                    f"AB6[{cs}]",
                    f"stage 21 optimum on the 1.4 um geometry, {cs}, h",
                    "dense true-min BV gain over initial",
                    gains["niwa"],
                    gains[cs],
                    "V",
                    {
                        "claimed": "Niwa gain (recomputed; equals results/stage21 dense final_L0 - initial_L0)",
                        "real": f"{cs} gain (recomputed); theta: results/stage21/results.json final/theta",
                    },
                    f"BV {d6[cs]['initial']['bv_min']:.2f} -> {d6[cs]['final']['bv_min']:.2f} V",
                    coefficient_set=cs,
                    gain_V=gains[cs],
                )
            )
        r = gains["hatakeyama"] / gains["niwa"] - 1.0
        pred["P25c"] = {"gain_hatakeyama_V": gains["hatakeyama"], "gain_niwa_V": gains["niwa"],
                        "gain_niwa_recorded_V": gn, "rel_difference": r, "tolerance": TOL["P25c"],
                        "confirmed": bool(gains["hatakeyama"] > 0 and abs(r) <= TOL["P25c"])}  # fmt: skip
    # ---------------- H3
    h30, h31 = (
        g("dense", "hatakeyama", "initial", "0.8", 0),
        g("dense", "hatakeyama", "initial", "0.8", 1),
    )
    h3 = {}
    if h30 and h31 and not h30.get("error") and not h31.get("error"):
        h3 = {
            "bv_h": h30["bv_min"],
            "bv_h2": h31["bv_min"],
            "rel_change": rel(h31["bv_min"], h30["bv_min"]),
        }
        h3["ok"] = bool(h3["rel_change"] <= TOL["H3_rel"])
    # ---------------- H2
    unbr, h2_ok = {}, not errors and n_blas == 0
    for k, v in res.items():
        if v.get("error"):
            continue
        unbr["|".join(map(str, k))] = v["n_unbracketed"]
        h2_ok = h2_ok and bool(v["resolved_all"]) and v["stats"].get("failed_evals", 1) == 0
    crit = {"H1": bool(h1_ok), "H2": bool(h2_ok), "H3": bool(h3.get("ok"))}
    extra = {"H1": {"rel_errors": h1, "tol": TOL["H1_rel"], "n_checks": len(h1)},
             "H2": {"job_errors": errors, "on_entry_to_lines": n_blas,
                    "n_unbracketed_lines_flagged_above_2x_first_crossing": unbr},
             "H3": h3, "predictions": pred, "AB3": ab3, "AB6_gains": gains}  # fmt: skip
    return table, extra, crit


# ============================================================ figure
def plot_compare(figdir: Path, table: list, extra: dict, label: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figdir.mkdir(parents=True, exist_ok=True)
    byid = {r["id"]: r for r in table}
    cols = {"niwa": ("#9ecae1", "#08519c"), "hatakeyama": ("#fdbe85", "#d94801")}
    names = {"niwa": "Niwa 2015", "hatakeyama": "Hatakeyama 2004"}
    fig, ax = plt.subplots(1, 3, figsize=(18, 5.2), gridspec_kw={"width_ratios": [1.2, 3.2, 1.4]})

    def grouped(a, groups, getter, ylabel, title):
        x = np.arange(len(groups))
        w = 0.2
        for j, (cs, kind) in enumerate([(c, k) for c in SETS for k in (0, 1)]):
            vals = [getter(gr, cs)[kind] for gr in groups]
            xx = x + (j - 1.5) * w
            a.bar(xx, vals, w, color=cols[cs][kind],
                  label=f"{names[cs]} {'claimed' if kind == 0 else 'real'}")  # fmt: skip
            for xi, v in zip(xx, vals, strict=True):
                a.text(xi, v, f"{v:.0f}", ha="center", va="bottom", fontsize=6)
        a.set_xticks(x)
        a.set_ylabel(ylabel, fontsize=8)
        a.set_title(title, fontsize=9)
        a.spines[["top", "right"]].set_visible(False)

    def tick(a, groups, labels, getter):
        a.set_xticklabels(
            [f"{lb}\nartifact N {getter(g_, 'niwa')[2] * 100:+.0f} % / H {getter(g_, 'hatakeyama')[2] * 100:+.0f} %"
             for g_, lb in zip(groups, labels, strict=True)], fontsize=7)  # fmt: skip

    def pick(i_fmt):
        def f(gr, cs):
            r = byid[i_fmt.format(gr=gr, cs=cs)]
            return r["claimed"], r["real"], r["artifact_fraction"]

        return f

    f1 = pick("AB1[{cs}]".replace("{cs}", "{cs}"))
    grouped(
        ax[0],
        ["AB1"],
        lambda gr, cs: f1(gr, cs),
        "BV gain [V]",
        "AB1': stage-13 planar vs dense gain",
    )
    tick(ax[0], ["AB1"], ["stage 13"], lambda gr, cs: f1(gr, cs))
    f3 = pick("AB3[{gr}][hatakeyama]")

    def f3b(gr, cs):
        if cs == "hatakeyama":
            return f3(gr, cs)
        n = extra["AB3"][gr]
        return (
            n["bv_niwa"],
            min(n["bv_niwa"], n["V_pt"]),
            1.0 - min(n["bv_niwa"], n["V_pt"]) / n["bv_niwa"],
        )

    grp = [n for n in DESIGNS if n in extra["AB3"]]
    grouped(
        ax[1], grp, f3b, "voltage [V]", "AB3': avalanche BV (dense, h) vs blocking min(BV, V_pt)"
    )
    tick(ax[1], grp, [M23.DESIGN_LABEL[n] for n in grp], f3b)
    ax[1].legend(fontsize=7, frameon=False, ncol=2, loc="upper right")
    gains = extra["AB6_gains"]
    if gains:
        f6 = lambda gr, cs: (gains["niwa"], gains[cs], 1.0 - gains[cs] / gains["niwa"])  # noqa: E731
        grouped(ax[2], ["AB6"], f6, "BV gain [V]", "AB6: stage-21 gain (1.4 um)")
        tick(ax[2], ["AB6"], ["stage 21"], f6)
    fig.suptitle(
        "Stage 25: artifact budget under Niwa 2015 vs Hatakeyama 2004 (claimed: light, real: dark)",
        fontsize=10,
    )
    fig.text(0.99, 0.005, label, ha="right", fontsize=6)
    fig.tight_layout(rect=(0, 0.02, 1, 0.95))
    fig.savefig(figdir / "stage25_niwa_vs_hatakeyama.png", dpi=140)
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
    cites = Cites(REPO)
    thetas = {n: M23.theta_of(cites, n) for n in DESIGNS}
    th21 = [float(t) for t in cites.get(RES.format(S21_THETA[0]), *S21_THETA[1])]
    specs = []
    for cs in SETS:
        for n in ("initial", "13"):
            specs.append(("ab1", cs, n, thetas[n]))
        for n in DESIGNS:
            specs.append(("dense", cs, "0.8", n, thetas[n], 0))
        for n, th in (("initial", [0.0] * 5), ("final", th21)):
            specs.append(("dense", cs, "1.4", n, th, 0))
    for n, th in (("initial", [0.0] * 5), ("final", th21)):
        specs.append(("dense", "niwa", "1.4", n, th, 1))  # H1 only
    specs.append(("dense", "hatakeyama", "0.8", "initial", thetas["initial"], 1))  # H3
    specs.sort(key=lambda s: 0 if s[0] == "dense" and s[5] == 1 else 1)
    res = {}
    ctx = mp.get_context("fork")
    with ProcessPoolExecutor(args.workers, mp_context=ctx) as pool:
        futs = {pool.submit(run_job, s): s for s in specs}
        for f in as_completed(futs):
            s, r = futs[f], f.result()
            key = ("ab1", s[1], s[2]) if s[0] == "ab1" else ("dense", s[1], s[3], s[2], s[5])
            res[key] = r
            print(
                f"done {key} {r.get('seconds', float('nan')):.0f}s {r.get('error', '')}", flush=True
            )
    sys.stderr.flush()
    stderr = log_path.read_text(errors="replace") if log_path.exists() else ""
    n_blas = len([ln for ln in stderr.splitlines() if "On entry to" in ln])
    table, extra, crit = evaluate(res, cites, n_blas, thetas, th21)
    runtime = time.perf_counter() - t_start
    seen = {(u["file"], u["key"]) for u in cites.used}
    results = {
        "config": {"tolerances": TOL, "coefficient_sets": {"niwa": "SiC4H.niwa() (default)",
                   "hatakeyama": "SiC4H.hatakeyama(): ii_{a,b}_{n,p}_par_c (stage 12)"},
                   "theta_stage21_final": th21, "theta": thetas, "n_dense": N_DENSE,
                   "note": "ev.ch set after construction; library unchanged"},  # fmt: skip
        "table": table, **extra,
        "citations": [{"file": f, "key": k} for f, k in sorted(seen)],
        "jobs": {"|".join(map(str, k)): v for k, v in res.items()},
        "on_entry_to_lines": n_blas, "runtime_seconds": runtime,
    }  # fmt: skip
    figdir = (REPO / "figures") if out == DEFAULT_OUT.resolve() else out / "figures"
    label = f"stage 25 | commit {str(man.get('git_commit'))[:8]} | {man['timestamp_utc'][:19]}Z"
    try:
        plot_compare(figdir, table, extra, label)
    except Exception as exc:
        results["figures_error"] = f"{type(exc).__name__}: {exc}"
    (out / "results.json").write_text(json.dumps(results, indent=1, default=float) + "\n")
    summary = {"criteria": crit, "gate_passed": bool(all(crit.values())), "dev": bool(args.dev),
               "table": table, "H1": extra["H1"], "H2": extra["H2"], "H3": extra["H3"],
               "predictions": extra["predictions"], "runtime_seconds": runtime}  # fmt: skip
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    print(json.dumps(crit, indent=2))
    for r in table:
        af = r["artifact_fraction"]
        print(
            f"{r['id']:<28} claimed {r['claimed']:.6g} real {r['real']:.6g} {r['unit']} "
            f"artifact {'n/a' if af is None else format(af, '+.4f')}"
        )
    for k, v in extra["predictions"].items():
        print(k, "CONFIRMED" if v["confirmed"] else "REFUTED")
    print(
        f"H1 max rel err {max(extra['H1']['rel_errors'].values(), default=float('nan')):.3g}; "
        f"H3 {extra['H3']}"
    )
    print(f"runtime {runtime:.1f} s")


if __name__ == "__main__":
    main()
