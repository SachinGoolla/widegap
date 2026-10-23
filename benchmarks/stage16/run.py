"""Stage-16 benchmark (docs/benchmarks/stage16_channel_mobility.md): H1-H6 (gating), R1-R3.

Sourced Noguchi inversion-layer mobility (N17/N20) in the stage-12b / 13d DMOSFET, FORWARD
re-evaluation only (no re-optimization).  Devices: ``12b``, ``13d_initial``, ``13d_final`` (the
stage-15f adapter).  Bias: R_on,sp at V_GS = 15 V, V_DS = 0.1 V (stage 12); V_T by the stage-12
linear extrapolation at maximum transconductance on the stage-12 transfer grid (V_DS = 0.1 V).

Jobs (each in a fork-pool worker; never raises, every failure is recorded):

* ``transfer`` (level 0, three devices): the stage-12 bulk transfer sweep (warm starts, bulk
  V_T / R_on,sp), then one Picard solve per V_GS with the Noguchi mobility (channel V_T, R_on,sp),
  H2 (bulk-mobility code path = ``mode='bulk'`` Picard vs the bulk sweep), channel share, the
  mu_ch(x) profile and the clamp counts at V_GS = 15 V.
* ``ron`` (levels 0-2, d_ch list): R_on,sp at V_GS = V_ON by the direct bias march, bulk and
  channel, for H5 / H6 / R1 (h/4) / R2 (d_ch = 10 nm).
* ``h3`` (N_A = 1e16, 1e17, 3e17): long-channel lateral n-MOSFET (widegap.channel_mobility.
  LateralNMOS, L = 10 um), V_GS found by Brent so that N_S(mid) = 3e12 cm^-2, then
  mu_eff = I_D L / (q W N_S V_DS) against mu_ch(E_eff, N_S) at the mid-channel column.

H1 is evaluated against the real registry: the Noguchi entries live in sic4h_pending_16.json until
an independent verifier signs them; ``check_entry`` refuses them (``verified_by`` null).
``--dev-unverified-registry`` builds, for DEV runs only, a merged registry in which the pending
entries carry a stub verification; it is recorded in the outputs and refused when ``--out`` lies
under results/.

Writes <out>/{results,summary,manifest}.json (default results/stage16) and figures stage16_*.png
(repo figures/ for the default output, else <out>/figures).
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

from widegap.channel_mobility import (
    NOGUCHI_NAMES,
    ChannelMap,
    LateralNMOS,
    NoguchiParams,
    depletion_sheet,
    picard_solve,
)
from widegap.fv2d.mos import continue_bias
from widegap.manifest import build_manifest
from widegap.materials.ionization import Q
from widegap.materials.registry import ProvenanceError, Registry, check_entry
from widegap.sic import REGISTRY_PATH, SiC4H, TracedRegistry

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[2]
DEFAULT_OUT = REPO / "results" / "stage16"
PENDING_PATH = REGISTRY_PATH.with_name("sic4h_pending_16.json")
DEV_STUB = "DEV-STUB (not a verification)"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


S15F = _load("stage15f_run", REPO / "benchmarks/stage15f/run.py")  # fixed names: fork pickling
S12 = S15F.S12
CASES = S15F.CASES  # ("12b", "13d_initial", "13d_final")

# ---- study parameters (spec) ----------------------------------------------------------------
V_ON, V_LIN = S12.V_ON, S12.V_LIN  # 15 V gate; V_DS = 0.1 V
D_CH_UM = 0.020  # 20 nm (spec section 2); 0.010 / 0.040 for R2 / H6
H2_TOL, H3_TOL, H5_TOL, H6_TOL = 1e-12, 0.02, 0.03, 0.01
H3_NA = (1e16, 1e17, 3e17)
H3_NS, H3_VDS, H3_L_UM = 3e12, 0.05, 10.0
GAMMA_TABLE_NA = (1e16, 3e16, 1e17, 4e17)
GAMMA_N20 = (0.20, 0.34, 0.46, 0.54)  # N20 ms p. 6 (measured)
GAMMA_SPEC = (0.229, 0.333, 0.447, 0.578)  # spec section 1 (log10)
GAMMA_SPEC_LN = (0.092, 0.331, 0.593, 0.895)


def rss_mib():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def _lst(a):
    return np.asarray(a, dtype=float).tolist()


# ============================================================ registry
def registry_json(dev_unverified: bool):
    main = json.loads(REGISTRY_PATH.read_text())
    if not dev_unverified:
        return main["parameters"]
    ents = list(main["parameters"])
    have = {e["name"] for e in ents}
    for e in json.loads(PENDING_PATH.read_text())["parameters"]:
        if e["name"] not in have:
            ents.append({**e, "verified_by": DEV_STUB, "verified_value_match": True})
    return ents


def h1_check():
    """H1: every Noguchi entry exists in sic4h.json and passes check_entry (incl. verification)."""
    main = {e["name"]: e for e in json.loads(REGISTRY_PATH.read_text())["parameters"]}
    pend = {e["name"]: e for e in json.loads(PENDING_PATH.read_text())["parameters"]}
    rows = {}
    for n in NOGUCHI_NAMES:
        row = {"in_sic4h_json": n in main, "in_pending_file": n in pend}
        entry = main.get(n) or pend.get(n)
        try:
            check_entry(entry)
            row["check_entry"] = "ok"
        except (ProvenanceError, TypeError) as exc:
            row["check_entry"] = f"REFUSED: {exc}"
        row["usable_in_registry"] = bool(n in main and row["check_entry"] == "ok")
        rows[n] = row
    return {"entries": rows, "all_usable": all(r["usable_in_registry"] for r in rows.values())}


# ============================================================ device helpers
class Dev(S15F.Adapter):
    """Stage-15f adapter plus the arrays the channel map needs."""

    def __init__(self, sc, case, level):
        super().__init__(sc, case, level)
        if self.ev is None:
            self.nd, self.na = self.d.nd, self.d.na
        else:
            self.nd, self.na = self.ev.nd, self.ev.na
        self.level = level

    @property
    def nodes(self):
        return (self.d.dev if self.ev is None else self.ev.dev).nodes


def _setup(case, level, entries):
    reg = TracedRegistry(Registry(entries))
    sc = SiC4H(reg)
    P = NoguchiParams.from_registry(reg)
    dev = Dev(sc, case, level)
    return sc, P, dev


def _id_a_per_cm(dev, sol):
    return float(dev.current_a_per_cm(abs(dev.model.contact_current(sol, dev.nodes["drain"]))))


def _bias_to(dev, sc, X_eq, vg, vds):
    """Direct bias march (V_DS, then V_GS) from equilibrium; returns the final solution."""
    out = continue_bias(
        dev.model, X_eq, lambda p: dev.bc(p[0], 0.0, p[1]), [0.0, 0.0],
        [[0.0, sc.volt_s(vds)], [sc.volt_s(vg), sc.volt_s(vds)]], max_step=20.0,
    )  # fmt: skip
    return out[-1][1]


def _picard_record(r):
    return {"converged": r["converged"], "picard_converged": r["picard_converged"],
            "newton_all_converged": r["newton_all_converged"], "iterations": r["iterations"],
            "final_rel": r["final_rel"], "newton_its": r["newton_its"],
            "rel_history": r["rel_history"]}  # fmt: skip


def _profile(cmap, r, sc):
    """mu_ch(x) and its ingredients on the gated channel columns."""
    cols = list(cmap.cols)
    keep = [k for k, i in enumerate(cols) if cmap.col_gated[i]]
    c = r["columns"]
    x_um = [float(sc.length_cm(cmap.mesh.x[cols[k]])) * 1e4 for k in keep]
    sel = lambda a: _lst(np.asarray(a)[keep])  # noqa: E731
    return {"x_um": x_um, "N_S_cm2": sel(r["ns_all"][cols]), "E_eff_MV_cm": sel(c["e_eff"]),
            "mu_ch": sel(c["mu"]), "mu_ph": sel(c["mu_ph"]), "mu_sr": sel(c["mu_sr"]),
            "mu_c": sel(c["mu_c"]), "N_A_cm3": sel(cmap.col_na_cm3[cols])}  # fmt: skip


# ============================================================ jobs
def job_transfer(args):
    """Level-0 (or given level) transfer sweeps: bulk (stage-12 path) and channel (Picard)."""
    case, level, entries = args
    t0 = time.perf_counter()
    rec = {"kind": "transfer", "case": case, "level": level, "error": None}
    try:
        sc, P, dev = _setup(case, level, entries)
        cmap = ChannelMap(sc, dev.model, dev.nd, dev.na, D_CH_UM,
                          gate_x_max=dev.prm.gate_x1)  # fmt: skip
        eq = dev.equilibrium()
        rec["equilibrium_converged"] = bool(eq.report.converged)
        vgs = list(S12.TRANSFER_VG)
        pts, out = S12.sweep(dev, sc, eq.X, vgs, [V_LIN], "transfer")
        bulk_ok = all(p["converged"] for p in pts)
        idb = [abs(p["I_D_scaled"]) for p in pts]
        vt_b, gm_b = S12.vt_extract([p["V_GS"] for p in pts], idb)
        i15 = next(p for p in pts if p["V_GS"] == V_ON)
        ron_b = float(V_LIN * dev.x_cm / i15["I_D_A_per_cm"])
        rec["bulk"] = {"V_T_V": vt_b, "gm_max_scaled": gm_b, "R_on_sp_ohm_cm2": ron_b,
                       "converged_all": bool(bulk_ok), "I_D_A_per_cm_15V": i15["I_D_A_per_cm"]}  # fmt: skip
        vds_s = float(sc.volt_s(V_LIN))
        sols = {vg: out[k + 1][1] for k, vg in enumerate(vgs)}
        # ---- channel Picard at every V_GS
        ch_pts, picard = [], []
        keep15 = None
        for vg in vgs:
            bc = dev.bc(float(sc.volt_s(vg)), 0.0, vds_s)
            r = picard_solve(cmap, P, sols[vg].X, bc)
            picard.append({"V_GS": vg, **_picard_record(r)})
            ok = r["converged"]
            i_a = _id_a_per_cm(dev, r["sol"]) if r["sol"] is not None and ok else float("nan")
            ch_pts.append({"V_GS": vg, "I_D_A_per_cm": i_a,
                           "I_D_scaled": abs(dev.model.contact_current(r["sol"], dev.nodes["drain"])) if ok else float("nan")})  # fmt: skip
            if vg == V_ON:
                keep15 = r
        all_ok = all(p["converged"] for p in picard)
        rec["picard"] = picard
        if all_ok:
            vt_c, gm_c = S12.vt_extract(vgs, [p["I_D_scaled"] for p in ch_pts])
            ron_c = float(
                V_LIN * dev.x_cm / next(p for p in ch_pts if p["V_GS"] == V_ON)["I_D_A_per_cm"]
            )
            rec["channel"] = {"V_T_V": vt_c, "gm_max_scaled": gm_c, "R_on_sp_ohm_cm2": ron_c,
                              "points": ch_pts}  # fmt: skip
        # ---- channel share / counts / profile at V_GS = 15 V
        if keep15 is not None and keep15["converged"]:
            rec["share_channel_15V"] = cmap.drop_share(keep15["sol"].phi_n, vds_s)
            rec["counts_15V"] = keep15["counts"]
            rec["profile_15V"] = _profile(cmap, keep15, sc)
            c = keep15["columns"]
            rec["gamma_extrapolation_15V"] = {
                "N_S_max_cm2": float(np.nanmax(keep15["ns_all"])),
                "E_eff_range_MV_cm": [float(c["e_eff"].min()), float(c["e_eff"].max())]}  # fmt: skip
            # bulk reference share (same quantity on the bulk solution)
            cmap.model.mu_edge = cmap.mu_bulk.copy()
            rec["share_bulk_15V"] = cmap.drop_share(sols[V_ON].phi_n, vds_s)
        # ---- H2: bulk-mobility code path of the same driver
        cmap.model.mu_edge = cmap.mu_bulk.copy()
        rb = picard_solve(cmap, P, sols[V_ON].X, dev.bc(float(sc.volt_s(V_ON)), 0.0, vds_s),
                          mode="bulk")  # fmt: skip
        ron_path = (float(V_LIN * dev.x_cm / _id_a_per_cm(dev, rb["sol"]))
                    if rb["converged"] else float("nan"))  # fmt: skip
        rec["h2"] = {"R_on_sp_code_path": ron_path, "R_on_sp_stage12_sweep": ron_b,
                     "rel_diff": abs(ron_path / ron_b - 1.0), "converged": bool(rb["converged"]),
                     "picard": _picard_record(rb)}  # fmt: skip
        cmap.model.mu_edge = cmap.mu_bulk.copy()
        rec["converged_all"] = bool(bulk_ok and all_ok and eq.report.converged)
    except Exception as exc:  # recorded, never hidden
        rec["error"] = f"{type(exc).__name__}: {exc}"
    rec["seconds"] = time.perf_counter() - t0
    rec["rss_mib"] = rss_mib()
    return rec


def job_ron(args):
    """R_on,sp at V_GS = V_ON by the direct march: bulk and channel for each d_ch (nm)."""
    case, level, d_ch_nm, entries = args
    t0 = time.perf_counter()
    rec = {"kind": "ron", "case": case, "level": level, "error": None, "d_ch_nm": list(d_ch_nm)}
    try:
        sc, P, dev = _setup(case, level, entries)
        eq = dev.equilibrium()
        vds_s = float(sc.volt_s(V_LIN))
        sol = _bias_to(dev, sc, eq.X, V_ON, V_LIN)
        rec["bulk_converged"] = bool(eq.report.converged and sol.report.converged)
        rec["nodes"] = int(dev.model.mesh.n_nodes)
        rec["bulk"] = {"R_on_sp_ohm_cm2": float(V_LIN * dev.x_cm / _id_a_per_cm(dev, sol))}
        rec["channel"] = {}
        for dnm in d_ch_nm:
            cmap = ChannelMap(sc, dev.model, dev.nd, dev.na, dnm * 1e-3,
                              gate_x_max=dev.prm.gate_x1)  # fmt: skip
            dev.model.mu_edge = cmap.mu_bulk.copy()
            r = picard_solve(cmap, P, sol.X, dev.bc(float(sc.volt_s(V_ON)), 0.0, vds_s))
            ent = {"picard": _picard_record(r), "n_edges": int(len(cmap.edges)),
                   "nodes_within_d_ch_top_column": int(
                       len(next(iter(cmap.col_nodes.values())))
                   )}  # fmt: skip
            if r["converged"]:
                ent["R_on_sp_ohm_cm2"] = float(V_LIN * dev.x_cm / _id_a_per_cm(dev, r["sol"]))
                ent["share_channel"] = cmap.drop_share(r["sol"].phi_n, vds_s)
                ent["counts"] = r["counts"]
                ent["profile"] = _profile(cmap, r, sc)
            rec["channel"][str(dnm)] = ent
            dev.model.mu_edge = cmap.mu_bulk.copy()
        rec["share_bulk"] = cmap.drop_share(sol.phi_n, vds_s)
    except Exception as exc:
        rec["error"] = f"{type(exc).__name__}: {exc}"
    rec["seconds"] = time.perf_counter() - t0
    rec["rss_mib"] = rss_mib()
    return rec


def job_h3(args):
    """H3: long-channel lateral n-MOSFET at N_S(mid) = 3e12 cm^-2."""
    na, entries = args
    t0 = time.perf_counter()
    rec = {"kind": "h3", "N_A_cm3": na, "error": None}
    try:
        reg = TracedRegistry(Registry(entries))
        sc = SiC4H(reg)
        P = NoguchiParams.from_registry(reg)
        dev = LateralNMOS(sc, na, L_um=H3_L_UM)
        cmap = ChannelMap(sc, dev.model, dev.nd, dev.na, D_CH_UM,
                          gate_x_max=dev.gate_x[1])  # fmt: skip
        eq = dev.equilibrium()
        imid = int(np.argmin(np.abs(dev.mesh.x - dev.x_mid_s)))
        kmid = list(cmap.cols).index(imid)
        vds_s = float(sc.volt_s(H3_VDS))
        cache = {}

        def solve(vg):
            out = continue_bias(
                dev.model, eq.X, lambda p: dev.bc(p[0], 0.0, p[1]), [0.0, 0.0],
                [[0.0, vds_s], [float(sc.volt_s(vg)), vds_s]], max_step=20.0,
            )  # fmt: skip
            dev.model.mu_edge = cmap.mu_bulk.copy()
            r = picard_solve(cmap, P, out[-1][1].X, dev.bc(float(sc.volt_s(vg)), 0.0, vds_s))
            cache[vg] = r
            return r

        def f(vg):
            r = solve(vg)
            if not r["converged"]:
                raise RuntimeError(f"Picard/Newton not converged at V_GS = {vg}")
            return np.log(r["ns_all"][imid] / H3_NS)

        vg_star = brentq(f, 2.0, 30.0, xtol=1e-9, rtol=1e-12)
        r = cache.get(vg_star) or solve(vg_star)
        ns = float(r["ns_all"][imid])
        c = r["columns"]
        mu_ch = float(c["mu"][kmid])
        i_a = dev.drain_current_a_per_cm(r["sol"])
        mu_eff = i_a * (H3_L_UM * 1e-4) / (Q * ns * H3_VDS)
        # local diagnostic from the surface quasi-Fermi gradient next to the mid column
        pn = r["sol"].phi_n * sc.vt
        s0, s1 = cmap.surf_node[imid], cmap.surf_node[imid + 1]
        dx = float(sc.length_cm(dev.mesh.x[imid + 1] - dev.mesh.x[imid]))
        ns_loc = 0.5 * (r["ns_all"][imid] + r["ns_all"][imid + 1])
        mu_loc = i_a / (Q * ns_loc * abs(pn[s1] - pn[s0]) / dx)
        # channel-integrated drop between the two n+ edges (diagnostic)
        gcols = [i for i in cmap.cols if dev.x_s0_s <= dev.mesh.x[i] <= dev.x_d0_s]
        v_ch = float(abs(pn[cmap.surf_node[gcols[0]]] - pn[cmap.surf_node[gcols[-1]]]))
        rec.update({
            "V_GS_V": float(vg_star), "N_S_mid_cm2": ns, "mu_ch_mid": mu_ch, "mu_eff": float(mu_eff),
            "rel_diff": float(abs(mu_eff / mu_ch - 1.0)), "mu_local_gradient": float(mu_loc),
            "rel_diff_local": float(abs(mu_loc / mu_ch - 1.0)),
            "E_eff_mid_MV_cm": float(c["e_eff"][kmid]), "N_dep_cm2": float(c["n_dep"][kmid]),
            "gamma_mid": float(c["gamma"][kmid]), "I_D_A_per_cm": float(i_a),
            "V_channel_drop_V": v_ch, "V_DS_V": H3_VDS, "picard": _picard_record(r),
            "mesh": [int(dev.mesh.nx), int(dev.mesh.ny)], "x_mid_um": float(sc.length_cm(dev.mesh.x[imid])) * 1e4,
            "L_um": H3_L_UM, "converged_all": bool(r["converged"] and eq.report.converged),
        })  # fmt: skip
    except Exception as exc:
        rec["error"] = f"{type(exc).__name__}: {exc}"
    rec["seconds"] = time.perf_counter() - t0
    rec["rss_mib"] = rss_mib()
    return rec


# ============================================================ R3 (no solves)
def r3_table(entries):
    reg = TracedRegistry(Registry(entries))
    sc = SiC4H(reg)
    P = NoguchiParams.from_registry(reg)
    na = np.array(GAMMA_TABLE_NA)
    nd = depletion_sheet(sc, na)
    g10 = P.C * np.log10(nd / P.ndep0) + P.D
    gln = P.C * np.log(nd / P.ndep0) + P.D
    e = np.array([0.02, 0.05, 0.1, 0.2, 0.5, 1.0])
    mu_ph = P.ph_coeff * e**P.ph_exp
    fig6 = float(P.ph_coeff * 0.1**P.ph_exp)
    mu_v_per_cm = float(P.ph_coeff * (0.1e6) ** P.ph_exp)
    # Fig. 4 anchors: mu_Coulomb,0 and gamma at N_dep = 3.5e12
    r = 3.5
    return {
        "gamma_vs_NA": {"N_A_cm3": _lst(na), "N_dep_cm2": _lst(nd), "gamma_log10": _lst(g10),
                         "gamma_ln_reading": _lst(gln), "gamma_N20_measured": list(GAMMA_N20),
                         "gamma_spec_log10": list(GAMMA_SPEC), "gamma_spec_ln": list(GAMMA_SPEC_LN),
                         "max_abs_diff_log10_vs_spec": float(np.max(np.abs(g10 - GAMMA_SPEC))),
                         "max_abs_diff_log10_vs_N20_measured": float(np.max(np.abs(g10 - GAMMA_N20))),
                         "max_abs_diff_ln_vs_N20_measured": float(np.max(np.abs(gln - GAMMA_N20)))},
        "mu_ph": {"E_eff_MV_cm": _lst(e), "mu_ph_cm2Vs": _lst(mu_ph),
                  "N17_Fig6_at_0p1_MV_cm": fig6, "if_E_in_V_per_cm_at_0p1_MV_cm": mu_v_per_cm},
        "N20_Fig4b_anchor_Ndep_3p5e12": {
            "gamma_log10": float(P.C * np.log10(r) + P.D), "gamma_ln": float(P.C * np.log(r) + P.D),
            "mu_coulomb0": float(P.A * r**P.B), "figure_reads": "gamma ~ 0.55, mu_C0 ~ 19"},
        "constants_used": {"eps_r": sc.eps_r, "n_i_cm3": sc.n_i, "kT_q_V": sc.vt},
    }  # fmt: skip


# ============================================================ figures
def figures(res, figdir: Path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figdir.mkdir(parents=True, exist_ok=True)
    fig, axs = plt.subplots(1, 3, figsize=(13, 3.8))
    for ax, case in zip(axs, CASES, strict=True):
        t = res["runs"].get(f"transfer|{case}|0")
        if not t or "profile_15V" not in t:
            ax.set_title(f"{case} (no data)")
            continue
        p = t["profile_15V"]
        x = np.array(p["x_um"])
        ax.plot(x, p["mu_ch"], "k-o", ms=3, label=r"$\mu_{ch}$")
        ax.plot(x, p["mu_ph"], "--", label=r"$\mu_{ph}$")
        ax.plot(x, p["mu_sr"], "--", label=r"$\mu_{SR}$")
        ax.plot(x, p["mu_c"], "--", label=r"$\mu_{C}$")
        ax.set_yscale("log")
        ax.set_xlabel("x [um] (gated p-type surface columns)")
        ax.set_ylabel(r"mobility [cm$^2$/Vs]")
        ax.set_title(f"{case}, V_GS = 15 V, level 0")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7, loc="lower right")
    fig.tight_layout()
    fig.savefig(figdir / "stage16_mu_ch_profile.png", dpi=130)
    plt.close(fig)
    rows = res["R1_table"]
    fig, ax = plt.subplots(figsize=(7.5, 4))
    w = 0.38
    xs = np.arange(len(rows))
    rb = [r["R_on_bulk_ohm_cm2"] * 1e4 if r["R_on_bulk_ohm_cm2"] else np.nan for r in rows]
    rc = [r["R_on_channel_ohm_cm2"] * 1e4 if r["R_on_channel_ohm_cm2"] else np.nan for r in rows]
    ax.bar(xs - w / 2, rb, w, label="bulk Hall mobility")
    ax.bar(xs + w / 2, rc, w, label="Noguchi channel mobility")
    ax.set_xticks(xs, [f"{r['case']}\nlevel {r['level']}" for r in rows], fontsize=8)
    ax.set_ylabel(r"$R_{on,sp}$ [$10^{-4}\ \Omega\,$cm$^2$]")
    ax.set_title("R_on,sp (V_GS = 15 V, V_DS = 0.1 V): bulk vs channel mobility")
    ax.legend()
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(figdir / "stage16_ron_comparison.png", dpi=130)
    plt.close(fig)


# ============================================================ main
def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--workers", type=int, default=min(4, os.cpu_count() or 1))
    ap.add_argument("--dev-unverified-registry", action="store_true",
                    help="DEV ONLY: merge the pending 16 entries with a stub verification")  # fmt: skip
    ap.add_argument("--quick", action="store_true", help="DEV ONLY: level 0 / 1 only (skip h/4)")
    args = ap.parse_args(argv)
    out = args.out.resolve()
    if (args.dev_unverified_registry or args.quick) and (REPO / "results") in (out, *out.parents):
        raise SystemExit("refusing --dev-unverified-registry/--quick with --out under results/")
    out.mkdir(parents=True, exist_ok=True)
    man = build_manifest(REPO, sys.argv, [])
    man["dev_unverified_registry"] = bool(args.dev_unverified_registry)
    (out / "manifest.json").write_text(json.dumps(man, indent=2) + "\n")
    t_start = time.perf_counter()
    h1 = h1_check()
    entries = registry_json(args.dev_unverified_registry)
    usable = Registry(entries).usable("noguchi_coul_A")
    if not usable:
        crit = {"H1": False}
        res = {"H1": h1, "note": "Noguchi entries unusable (unverified): no channel run possible"}
        (out / "results.json").write_text(json.dumps(res, indent=2) + "\n")
        (out / "summary.json").write_text(
            json.dumps({"criteria": crit, "gate_passed": False}, indent=2) + "\n"
        )
        print(json.dumps(crit))
        return
    jobs = []
    for c in CASES:
        jobs.append((job_transfer, ("transfer", c, 0), (c, 0, entries)))
    for lv in (0, 1) if args.quick else (0, 1, 2):
        for c in CASES:
            if c != "12b" and lv == 1:
                continue
            dch = (10, 20, 40) if (c == "12b" and lv == 1) else (20,)
            jobs.append((job_ron, ("ron", c, lv), (c, lv, dch, entries)))
    for na in H3_NA:
        jobs.append((job_h3, ("h3", na, 0), (na, entries)))
    # heaviest first
    order = {2: 0, 1: 1, 0: 2}
    jobs.sort(key=lambda j: 3 if j[1][0] == "h3" else order[j[1][2]])
    runs: dict = {}
    ctx = mp.get_context("fork")
    with ProcessPoolExecutor(args.workers, mp_context=ctx) as pool:
        futs = {pool.submit(fn, a): key for fn, key, a in jobs}
        for f in as_completed(futs):
            r = f.result()
            runs["|".join(str(k) for k in futs[f])] = r
            print(
                f"done {futs[f]} {r['seconds']:.0f}s rss={r['rss_mib']:.0f}MiB err={r['error']}",
                flush=True,
            )
    r3 = r3_table(entries)
    # ---- criteria
    crit, vals = {"H1": bool(h1["all_usable"])}, {}
    t12 = runs.get("transfer|12b|0", {})
    h2 = t12.get("h2")
    vals["H2"] = h2
    crit["H2"] = bool(h2 and h2["converged"] and h2["rel_diff"] <= H2_TOL)
    h3rows = [runs[f"h3|{na}|0"] for na in H3_NA]
    vals["H3"] = [{k: r.get(k) for k in ("N_A_cm3", "V_GS_V", "N_S_mid_cm2", "mu_ch_mid", "mu_eff",
                   "rel_diff", "mu_local_gradient", "rel_diff_local", "E_eff_mid_MV_cm", "error")}
                  for r in h3rows]  # fmt: skip
    crit["H3"] = bool(all(r["error"] is None and r.get("converged_all") and r["rel_diff"] <= H3_TOL
                          and abs(r["N_S_mid_cm2"] / H3_NS - 1.0) <= 0.01 for r in h3rows))  # fmt: skip
    # H4: every Picard evaluation converged
    picards = []
    for k, r in runs.items():
        if r.get("error"):
            picards.append((k, False))
        if r["kind"] == "transfer":
            picards += [(f"{k}@{p['V_GS']}", p["converged"]) for p in r.get("picard", [])]
            if r.get("h2"):
                picards.append((f"{k}@H2", r["h2"]["converged"]))
        elif r["kind"] == "ron":
            picards += [
                (f"{k}@{d}", e["picard"]["converged"]) for d, e in r.get("channel", {}).items()
            ]
            picards.append((f"{k}@bulk", r.get("bulk_converged", False)))
        else:
            if r.get("picard"):
                picards.append((k, r["picard"]["converged"]))
    crit["H4"] = bool(picards and all(ok for _, ok in picards))
    vals["H4"] = {"n_evaluations": len(picards), "failed": [k for k, ok in picards if not ok]}

    def ron(case, lv, d=None):
        r = runs.get(f"ron|{case}|{lv}")
        if not r or r.get("error"):
            return None, None
        ch = r["channel"].get(str(d or 20), {})
        return r["bulk"]["R_on_sp_ohm_cm2"], ch.get("R_on_sp_ohm_cm2")

    r0, r1 = ron("12b", 0)[1], ron("12b", 1)[1]
    vals["H5"] = {"R_on_channel_h": r0, "R_on_channel_h2": r1,
                  "rel_diff": None if not (r0 and r1) else abs(r1 / r0 - 1.0)}  # fmt: skip
    crit["H5"] = bool(r0 and r1 and abs(r1 / r0 - 1.0) <= H5_TOL)
    r20, r40, r10 = ron("12b", 1, 20)[1], ron("12b", 1, 40)[1], ron("12b", 1, 10)[1]
    vals["H6"] = {"R_on_20nm": r20, "R_on_40nm": r40, "R_on_10nm": r10,
                  "rel_diff_40_vs_20": None if not (r20 and r40) else abs(r40 / r20 - 1.0),
                  "rel_diff_10_vs_20": None if not (r20 and r10) else abs(r10 / r20 - 1.0)}  # fmt: skip
    crit["H6"] = bool(r20 and r40 and abs(r40 / r20 - 1.0) <= H6_TOL)
    # ---- R1 table
    r1_rows = []
    for c in CASES:
        t = runs.get(f"transfer|{c}|0", {})
        for lv in (0, 1, 2):
            rr = runs.get(f"ron|{c}|{lv}")
            if rr is None or rr.get("error"):
                continue
            ch = rr["channel"].get("20", {})
            row = {"case": c, "level": lv, "R_on_bulk_ohm_cm2": rr["bulk"]["R_on_sp_ohm_cm2"],
                   "R_on_channel_ohm_cm2": ch.get("R_on_sp_ohm_cm2"),
                   "share_channel_drop": ch.get("share_channel"), "share_bulk_drop": rr.get("share_bulk")}  # fmt: skip
            if row["R_on_channel_ohm_cm2"]:
                row["ratio"] = row["R_on_channel_ohm_cm2"] / row["R_on_bulk_ohm_cm2"]
                row["share_delta_1_minus_Rb_over_Rc"] = 1.0 - 1.0 / row["ratio"]
            if lv == 0 and t.get("channel"):
                row["V_T_bulk_V"], row["V_T_channel_V"] = t["bulk"]["V_T_V"], t["channel"]["V_T_V"]
            r1_rows.append(row)
    r2 = {c: {"counts_15V": runs.get(f"transfer|{c}|0", {}).get("counts_15V"),
              "gamma_extrapolation_15V": runs.get(f"transfer|{c}|0", {}).get("gamma_extrapolation_15V"),
              "profile_15V": runs.get(f"transfer|{c}|0", {}).get("profile_15V")} for c in CASES}  # fmt: skip
    r2["d_ch_10nm_12b_level1"] = {"R_on_channel": r10, "R_on_20nm": r20}
    runtime = time.perf_counter() - t_start
    res = {"config": {"V_ON": V_ON, "V_LIN": V_LIN, "d_ch_nm": 20, "tolerances": {
        "H2": H2_TOL, "H3": H3_TOL, "H5": H5_TOL, "H6": H6_TOL, "picard_tol": 1e-10, "picard_max": 60},
        "dev_unverified_registry": bool(args.dev_unverified_registry), "quick": bool(args.quick),
        "levels": "0 (=h), 1 (=h/2), 2 (=h/4)"},
        "H1": h1, "runs": runs, "R1_table": r1_rows, "R2": r2, "R3": r3, "runtime_seconds": runtime}  # fmt: skip
    summary = {"criteria": crit, "gate_passed": bool(all(crit.values())),
               "dev_unverified_registry": bool(args.dev_unverified_registry),
               "values": {**vals, "R1": r1_rows, "R3_gamma": r3["gamma_vs_NA"]},
               "runtime_seconds": runtime}  # fmt: skip
    (out / "results.json").write_text(json.dumps(res, indent=2, default=float) + "\n")
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    figdir = REPO / "figures" if out == DEFAULT_OUT.resolve() else out / "figures"
    try:
        figures(res, figdir)
    except Exception as exc:  # plotting must not hide results
        print("figure error:", exc)
    print(json.dumps(crit, indent=2))
    print(f"runtime {runtime:.1f} s")


if __name__ == "__main__":
    main()
