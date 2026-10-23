"""Stage-15f benchmark (docs/benchmarks/stage15_caveat_closure.md, section 15f): V1, V2 (gating), R4.

Saturation velocity from the Ishikawa thesis (Table 3.1, Eq. 3.3) against Khan & Cooper, in the
stage-12b / 13d DD model (``SiC4H.vsat_set``).  Devices (mesh h = level 0):

* ``12b``          : the stage-12b half-cell (abrupt junctions, theta = 0);
* ``13d_initial``  : the stage-13d rounded device (r = 0.1 um) at theta = 0;
* ``13d_final``    : the same at the stage-13d final theta (results/stage13d/results.json).

R4 per device and v_sat set: R_on,sp (V_GS = 15 V, V_DS = 0.1 V), V_T (stage-12 linear
extrapolation at maximum transconductance on the stage-12 transfer grid, V_DS = 0.1 V) and I_D at
V_DS = 15 V (V_GS = 15 V, stage-12 output grid; the spec does not give V_GS, 15 V is the stage-12
on-state gate bias).  The sweeps are the stage-12 ``sweep`` routine.

Variants (``vsat``; stage-15 Amendment 1):
* ``khan_cooper``     : the existing default;
* ``ishikawa_beta1.0``: Ishikawa isotropic v_sat = 1.7e7 cm/s, beta = 1.0 (the model's lower bound);
* ``ishikawa_beta1.2``: the same v_sat, beta = 1.2 (Khan & Cooper value).
The betas are MODEL CHOICES, not registry values: the thesis alpha at N_D <= 1e16 is < 1, which the
stage-9a Caughey-Thomas implementation (ct_factor, beta >= 1) does not support.  V2: every solve of
both Ishikawa variants converges.
V1 is evaluated against the real registry (src/widegap/materials/data/sic4h.json).  The Ishikawa
entries are in sic4h_pending_15f.json until an independent verifier signs them; the check refuses
them (``verified_by`` null).  ``--dev-unverified-registry`` builds, for DEV runs only, a merged
registry in which the pending entries carry a stub verification; it is recorded in the outputs and
refused when ``--out`` lies under results/.

Writes <out>/{results,summary,manifest}.json (default results/stage15f).
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

from widegap.adjoint_sic_oxide import OxideAdjoint
from widegap.manifest import build_manifest
from widegap.materials.registry import ProvenanceError, Registry, check_entry
from widegap.sic import REGISTRY_PATH, SiC4H, SiCDMOS, TracedRegistry

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[2]
DEFAULT_OUT = REPO / "results" / "stage15f"
R13D = REPO / "results" / "stage13d"
PENDING_PATH = REGISTRY_PATH.with_name("sic4h_pending_15f.json")
V_GS_ID, V_DS_ID = 15.0, 15.0
DEV_STUB = "DEV-STUB (not a verification)"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


S12 = _load("stage12_run", REPO / "benchmarks/stage12/run.py")  # fixed names: fork pickling
M13D = _load("stage13d_run", REPO / "benchmarks/stage13d/run.py")
DESIGN, TH0, R_UM = M13D.DESIGN, M13D.TH0, M13D.M.R_UM
CASES = ("12b", "13d_initial", "13d_final")
VARIANTS = ("khan_cooper", "ishikawa_beta1.0", "ishikawa_beta1.2")
BETA = {
    "ishikawa_beta1.0": 1.0,
    "ishikawa_beta1.2": 1.2,
}  # MODEL CHOICES (Amendment 1), not registry values
NEW_NAMES = [
    "vsat_n_par_c",
    "vsat_n_perp_c_ishikawa",
    *[
        f"ct_beta_n_{o}_c_ishikawa_{n}"
        for o in ("par", "perp")
        for n in ("2e15", "5e15", "1e16", "3e16")
    ],
]


# ============================================================ registry
def registry_json(dev_unverified: bool):
    """Registry JSON path/entries for the runs (merged with stubbed pending entries if dev)."""
    main = json.loads(REGISTRY_PATH.read_text())
    if not dev_unverified:
        return main["parameters"]
    ents = list(main["parameters"])
    have = {e["name"] for e in ents}
    for e in json.loads(PENDING_PATH.read_text())["parameters"]:
        if e["name"] not in have:
            ents.append({**e, "verified_by": DEV_STUB, "verified_value_match": True})
    return ents


def v1_check():
    """V1: every new entry exists in sic4h.json and passes check_entry (incl. independent verification)."""
    main = {e["name"]: e for e in json.loads(REGISTRY_PATH.read_text())["parameters"]}
    pend = {e["name"]: e for e in json.loads(PENDING_PATH.read_text())["parameters"]}
    rows = {}
    for n in NEW_NAMES:
        row = {"in_sic4h_json": n in main, "in_pending_file": n in pend}
        entry = main.get(n) or pend.get(n)
        try:
            check_entry(entry)
            row["check_entry"] = "ok"
        except ProvenanceError as exc:
            row["check_entry"] = f"REFUSED: {exc}"
        row["usable_in_registry"] = bool(n in main and row["check_entry"] == "ok")
        rows[n] = row
    kc = {n: main[n] for n in ("vsat_n_perp_c", "ct_beta_n_perp_c")}
    return {"entries": rows, "all_usable": all(r["usable_in_registry"] for r in rows.values()),
            "khan_cooper_entries_unchanged": {n: {"value": e["value"], "verified_by": e["verified_by"]}
                                              for n, e in kc.items()}}  # fmt: skip


# ============================================================ device adapters and jobs
class Adapter:
    """What the stage-12 ``sweep`` needs from a device (model, bc, currents, scaling)."""

    def __init__(self, sc, case, level=0):
        self.sc, self.case = sc, case
        if case == "12b":
            self.d = SiCDMOS(sc, DESIGN, level)
            self.model, self.base, self.prm = self.d.model, self.d, self.d.prm
            self.ev = None
        else:
            self.ev = OxideAdjoint(sc, DESIGN, level, r_um=R_UM, cuts_um=(), tau_v=M13D.M.TAU_V)
            theta = (
                TH0
                if case == "13d_initial"
                else np.asarray(
                    json.loads((R13D / "results.json").read_text())["final"]["theta"], float
                )
            )
            self.ev.set_theta(theta)  # fmt: skip
            self.model, self.base, self.prm = self.ev.mod, self.ev.base, self.ev.prm

    @property
    def mob(self):
        return self.model.mob

    def bc(self, vg, vs, vd):
        return self.base.bc(vg, vs, vd)

    def equilibrium(self):
        return self.d.equilibrium() if self.ev is None else self.ev._eq()

    def currents(self, sol):
        n = (self.d.dev if self.ev is None else self.ev.dev).nodes
        mod = self.model
        return {"I_source": mod.contact_current(sol, n["source"]),
                "I_drain": mod.contact_current(sol, n["drain"]),
                "I_body": mod.contact_current(sol, n["body"])}  # fmt: skip

    def current_a_per_cm(self, i_scaled):
        return self.base.current_a_per_cm(i_scaled)

    @property
    def x_cm(self):
        return float(self.sc.length_cm(self.prm.x_max))


def job_case(args):
    """One device with one v_sat variant: R_on,sp, V_T, I_D(V_DS = 15 V); never raises."""
    case, variant, entries = args
    t0 = time.perf_counter()
    rec = {"case": case, "vsat": variant, "error": None, "converged_all": False}
    try:
        reg = TracedRegistry(Registry(entries))
        if variant == "khan_cooper":
            sc = SiC4H(reg)
        else:
            sc = SiC4H(reg, vsat_set="ishikawa", vsat_beta=BETA[variant])
        dev = Adapter(sc, case)
        rec["mobility"] = {"vsat_cm_s": sc.saturation()[0], "beta": dev.mob.beta_n,
                           "vsat_scaled": dev.mob.vsat_n}  # fmt: skip
        eq = dev.equilibrium()
        rec["equilibrium_converged"] = bool(eq.report.converged)
        pts, _ = S12.sweep(dev, sc, eq.X, S12.TRANSFER_VG, [S12.V_LIN], "transfer")
        vg = [p["V_GS"] for p in pts]
        idr = [abs(p["I_D_scaled"]) for p in pts]
        vt, gm = S12.vt_extract(vg, idr)
        i15 = next(p for p in pts if p["V_GS"] == S12.V_ON)
        rec["transfer"] = {"V_T_V": vt, "gm_max_scaled": gm, "points": pts,
                           "converged_all": all(p["converged"] for p in pts)}  # fmt: skip
        rec["R_on_sp_ohm_cm2"] = float(S12.V_LIN * dev.x_cm / i15["I_D_A_per_cm"])
        rec["R_on_point_converged"] = bool(i15["converged"])
        opts, _ = S12.sweep(dev, sc, eq.X, [V_GS_ID], S12.OUTPUT_VDS, "output")
        p15 = next(p for p in opts if abs(p["V_DS"] - V_DS_ID) < 1e-12)
        rec["output"] = {"points": opts, "converged_all": all(p["converged"] for p in opts)}
        rec["I_D_VDS15_A_per_cm"] = p15["I_D_A_per_cm"]
        rec["I_D_point_converged"] = bool(p15["converged"])
        rec["converged_all"] = bool(rec["equilibrium_converged"] and rec["transfer"]["converged_all"]
                                    and rec["output"]["converged_all"])  # fmt: skip
        if dev.ev is not None and variant == "khan_cooper":  # cross-check with the 13d R_on routine
            rec["R_on_sp_13d_routine_ohm_cm2"] = dev.ev.ron(dev.ev.theta, grad=False)["R"]
    except Exception as exc:  # recorded, never hidden
        rec["error"] = f"{type(exc).__name__}: {exc}"
    rec["seconds"] = time.perf_counter() - t0
    return rec


# ============================================================ main
def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--workers", type=int, default=min(9, os.cpu_count() or 1))
    ap.add_argument("--dev-unverified-registry", action="store_true",
                    help="DEV ONLY: merge the pending 15f entries with a stub verification")  # fmt: skip
    args = ap.parse_args(argv)
    out = args.out.resolve()
    if args.dev_unverified_registry and (REPO / "results") in (out, *out.parents):
        raise SystemExit("refusing --dev-unverified-registry with --out under results/")
    out.mkdir(parents=True, exist_ok=True)
    man = build_manifest(REPO, sys.argv, [])
    man["dev_unverified_registry"] = bool(args.dev_unverified_registry)
    (out / "manifest.json").write_text(json.dumps(man, indent=2) + "\n")
    t_start = time.perf_counter()
    v1 = v1_check()
    entries = registry_json(args.dev_unverified_registry)
    usable = Registry(entries).usable("vsat_n_perp_c_ishikawa")
    jobs = [(c, v, entries) for c in CASES for v in VARIANTS if usable or v == "khan_cooper"]
    results: dict = {}
    ctx = mp.get_context("fork")
    with ProcessPoolExecutor(args.workers, mp_context=ctx) as pool:
        futs = {pool.submit(job_case, a): (a[0], a[1]) for a in jobs}
        for f in as_completed(futs):
            r = f.result()
            results[futs[f]] = r
            print(f"done {futs[f]} {r['seconds']:.0f}s err={r['error']}", flush=True)
    table = []
    for c in CASES:
        row = {"case": c}
        for v in VARIANTS:
            r = results.get((c, v))
            row[v] = None if r is None else {
                "R_on_sp_ohm_cm2": r.get("R_on_sp_ohm_cm2"),
                "V_T_V": (r.get("transfer") or {}).get("V_T_V"),
                "I_D_VDS15_A_per_cm": r.get("I_D_VDS15_A_per_cm"),
                "converged_all": r["converged_all"], "error": r["error"]}  # fmt: skip
        kc = row["khan_cooper"]
        for v in VARIANTS[1:]:
            q = row[v]
            if q and kc and q["R_on_sp_ohm_cm2"] and kc["R_on_sp_ohm_cm2"]:
                q["ratio_to_khan_cooper"] = {
                    "R_on_sp": q["R_on_sp_ohm_cm2"] / kc["R_on_sp_ohm_cm2"],
                    "V_T_shift_V": q["V_T_V"] - kc["V_T_V"],
                    "I_D_VDS15": q["I_D_VDS15_A_per_cm"] / kc["I_D_VDS15_A_per_cm"]}  # fmt: skip
        table.append(row)
    v2_runs = [results.get((c, v)) for c in CASES for v in VARIANTS[1:]]
    v2 = bool(all(r is not None and r["converged_all"] and r["error"] is None for r in v2_runs))
    s13 = json.loads((R13D / "summary.json").read_text())["values"]
    consistency = {c: {"R_on_stage13d_summary": s13["R_initial" if c == "13d_initial" else "R_final"],
                       "R_on_here_khan_cooper": (results.get((c, "khan_cooper")) or {}).get("R_on_sp_ohm_cm2")}
                   for c in ("13d_initial", "13d_final")}  # fmt: skip
    crit = {"V1": bool(v1["all_usable"]), "V2": v2}
    res = {"config": {"V_GS_ID": V_GS_ID, "V_DS_ID": V_DS_ID, "transfer_vg": S12.TRANSFER_VG,
                      "output_vds": S12.OUTPUT_VDS, "V_LIN": S12.V_LIN, "V_ON": S12.V_ON, "level": 0,
                      "ishikawa_beta_model_choices": BETA,
                      "dev_unverified_registry": bool(args.dev_unverified_registry)},
           "V1": v1, "runs": {f"{c}|{v}": r for (c, v), r in sorted(results.items())},
           "R4_table": table, "consistency_13d": consistency,
           "runtime_seconds": time.perf_counter() - t_start}  # fmt: skip
    summary = {"criteria": crit, "gate_passed": all(crit.values()),
               "dev_unverified_registry": bool(args.dev_unverified_registry),
               "values": {"V1_all_usable": v1["all_usable"], "V2_errors": {f"{r['case']}|{r['vsat']}": r["error"] for r in v2_runs if r},
                "beta_note": "beta values of the Ishikawa variants are model choices (Amendment 1), not registry values",
                          "R4": table, "consistency_13d": consistency},
               "runtime_seconds": res["runtime_seconds"]}  # fmt: skip
    (out / "results.json").write_text(json.dumps(res, indent=2, default=float) + "\n")
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    print(json.dumps(crit, indent=2))
    print(f"runtime {res['runtime_seconds']:.1f} s")


if __name__ == "__main__":
    main()
