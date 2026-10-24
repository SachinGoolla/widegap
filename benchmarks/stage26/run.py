"""Stage-26 benchmark (docs/benchmarks/stage26_saha_comparison.md): C1-C3 (gating), R1.

Consistency check of the stage-16 channel-mobility DMOSFET model against ONE fabricated 4H-SiC
DMOSFET (Saha and Cooper, IEEE TED 54, 2786 (2007), "S07"): comparisons (a) cross-simulator
(constant mu = 20 cm^2/Vs vs S07 Table I), (b) closed-form mobility vs the S07 Fig. 12 fit, (c)
measured R_on,sp (6.95 mOhm cm^2) with a one-at-a-time and corner bracket.  Not a validation.

Model: ``widegap.sic.DMOSDesign`` with the stage-26a extensions (CSL layer, retrograde body,
gate overhang, explicit source-contact start), rounded junctions (r = 0.1 um as in every earlier
stage), V_GS = 20 V, V_DS = 0.1 V, channel/mobility arms by the stage-16 Picard driver.  Series
terms (``widegap.series``) are added after the solve.

Geometry mapping (recorded in results.json["interpretations"]); everything in um, half cell
x in [0, x_max] with x = 0 the JFET centre:

* x_max = pitch / 2; x_jfet = W_J / 2 (W_J full width).
* source-contact start x_cs and gate end x_g from the S07 cell (Fig. 8: gate length 3 um,
  gate-source separation 0.5 um, source contact 3 um, pitch 10 um): x_g = x_cs - 0.5,
  x_cs = x_max - L_contact (central: 2.0, 1.5, 5.0 -> 3.0 um contact; comparison (a): pitch 8,
  LS = 2 -> x_cs = 2.0, x_g = 1.5).  The gate therefore extends 1.0 um past the JFET edge
  (gate length 3 um at W_J = 1 um) and overhangs the n+ source by L_ov = x_g - x_src0.
* x_src0 = x_jfet + ch_len, with ch_len calibrated so that the NET-p surface span of the
  rounded doping (r = 0.1 um) equals L_ch ("electrical" channel length).  The literal nominal
  mapping (ch_len = L_ch between the 50 % points) leaves NO net-p surface column at N_A,s =
  2.5e17 (n+ tail of the 1e20 source, 5 corners of the pre-run failed with an empty channel) and
  is kept as the scan item ``nominal_l_ch``.  The n+ source extends to x_max - 0.04
  (p-body contact: bc_w 0.04, pad_w 0.02, so the modelled source contact length is
  L_contact - 0.04 um) and the source contact starts at x_cs (src_gap = x_cs - x_src0).
* W_J trend (1.5, 2.0 um): the S07 cell is widened keeping every other length (gate 1 um past
  the JFET edge, contact 1.5 um past it, 3 um contact): pitch = 9 + W_J.
* "W_J effective" (0.8 guess): x_jfet = 0.4 with the drawn layout (x_g, x_cs, pitch) fixed.

``--dev`` marks the run DEV (refused when --out lies under results/).
"""

# ruff: noqa: E402
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import itertools
import json
import multiprocessing as mp
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from widegap.adjoint_channel import D_CH_UM
from widegap.adjoint_sic import EvalError
from widegap.adjoint_sic_corner import DenseCuts, dense_bv, smooth_step
from widegap.bodyredesign import OffMarch, column_index_matrix, source_columns
from widegap.channel_mobility import (
    ChannelMap,
    NoguchiParams,
    depletion_sheet,
    eps_sic_f_per_cm,
    noguchi_mobility,
    picard_solve,
)
from widegap.manifest import build_manifest
from widegap.materials.ionization import Q
from widegap.offstate import CorrectedCornerAdjoint
from widegap.series import contact_resistance_sp, source_contact_length_cm
from widegap.sic import EPS0, DMOSDesign, SiC4H, TracedRegistry

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[2]
DEFAULT_OUT = REPO / "results" / "stage26"

# ---- pre-registered constants (docs/benchmarks/stage26_saha_comparison.md) ----------------
R_MEAS = 6.95e-3  # S07 measured R_on,sp [ohm cm^2] (1 um JFET)
R_MEAS_WJ = {1.0: 6.95e-3, 1.5: 7.5e-3, 2.0: 7.9e-3}
TABLE_I = {  # (rho_c [ohm cm^2], mobility) -> R_on,sp [ohm cm^2]; pitch 8 um design inferred
    (5e-4, "20"): 3.89e-3, (1e-3, "20"): 5.38e-3, (5e-4, "fit"): 4.53e-3, (1e-3, "fit"): 6.02e-3,
}  # fmt: skip
BV_S06_SIM = 1230.0
# S07 Fig. 12 fitted curve, read from the figure (APPROXIMATE): gate field [MV/cm] -> mu [cm^2/Vs]
FIG12_X = np.array([0.5, 1.0, 1.5, 2.0, 3.0, 4.0])
FIG12_Y = np.array([11.0, 13.6, 13.4, 11.5, 9.0, 7.0])
EPS_OX = 3.9  # owner-approved SiO2 (registry eps_r_SiO2_thermal is checked equal)
MU_CONST = 20.0
LEVEL_MAIN, LEVEL_FINE = 1, 2  # h/2 (reported / gated) and h/4
C2_TOL = 0.03
V_GS, V_DS = 20.0, 0.1
R_UM = 0.1  # junction rounding length (every earlier stage)
T_OX_UM = 0.052
T_OX_SIM_UM = 0.050
GATE_SEP_UM = 0.5
BC_W_UM, PAD_W_UM = 0.04, 0.02
N_S_FIT_MAX = 1.0e13
RANGES = {  # (low, high); central marked in CENTRAL
    "na_s": (2.5e17, 5e17), "na_b": (5e17, 1e18), "w_j_eff": (0.8, 1.0), "l_ch": (0.3, 0.4),
    "r_sub": (0.0, 0.7e-3), "rho_c": (5e-4, 1e-3),
}  # fmt: skip
CENTRAL = {"na_s": 5e17, "na_b": 1e18, "w_j_eff": 1.0, "l_ch": 0.3, "r_sub": 0.35e-3,
           "rho_c": 1e-3}  # fmt: skip
PDE_KEYS = ("na_s", "na_b", "w_j_eff", "l_ch")
SERIES_KEYS = ("rho_c", "r_sub")


def mdesign(cfg: dict):
    """(DMOSDesign, metadata) for a geometry/doping configuration ``cfg`` (um, cm^-3)."""
    x_max = cfg["x_max"]
    x_jfet = 0.5 * cfg["w_j_eff"]
    ch_len = electrical_ch_len(cfg) if "l_ch_el" in cfg else cfg["ch_len"]
    x_src0 = x_jfet + ch_len
    x_cs = x_max - cfg["l_contact"]
    x_g = x_cs - GATE_SEP_UM
    x_src1 = x_max - 0.5 * BC_W_UM - PAD_W_UM
    if not (x_src0 <= x_g < x_cs < x_src1):
        raise ValueError(f"inconsistent layout: x_src0 {x_src0}, x_g {x_g}, x_cs {x_cs}")
    d = DMOSDesign(
        x_jfet=x_jfet, ch_len=ch_len, src_w=x_src1 - x_src0, bc_w=BC_W_UM, pad_w=PAD_W_UM,
        d_src=0.25, d_body=0.6, w_drift=6.0, t_drain=0.5, t_ox=cfg["t_ox"],
        n_drift=1e16, n_jfet=1e17, n_body=cfg["na_s"], n_src=1e20, n_drain=1e19,
        mesh_kw=dict(ny_drift=16, ny_body=8), t_csl=cfg.get("t_csl", 0.3), n_csl=1e17,
        n_body_b=cfg.get("na_b"), l_ov=x_g - x_src0, src_gap=x_cs - x_src0,
    )  # fmt: skip
    meta = {"x_max": x_max, "x_jfet": x_jfet, "x_src0": x_src0, "x_gate_end": x_g,
            "x_contact_start": x_cs, "x_src1": x_src1, "l_ov": x_g - x_src0,
            "l_contact_model": x_src1 - x_cs, "ch_len_used": ch_len,
            "l_ch_electrical_target": cfg.get("l_ch_el")}  # fmt: skip
    return d, meta


def central_cfg(**kw) -> dict:
    c = dict(name="central", x_max=5.0, w_j_eff=1.0, l_ch_el=0.3, l_contact=3.0, na_s=5e17,
             na_b=1e18, t_ox=T_OX_UM, t_csl=0.3, r_um=R_UM)  # fmt: skip
    c.update(kw)
    return c


def sim_cfg(na, **kw) -> dict:
    """Comparison (a): S07 simulation conditions (8 um pitch, LS = 2 um, L_ch 0.5 um, 50 nm)."""
    c = dict(name=f"a_NA{na:.1e}", x_max=4.0, w_j_eff=1.0, l_ch_el=0.5, l_contact=2.0, na_s=na,
             na_b=None, t_ox=T_OX_SIM_UM, t_csl=0.3, r_um=R_UM)  # fmt: skip
    c.update(kw)
    return c


def electrical_ch_len(cfg: dict) -> float:
    """ch_len (50 % point separation) for which the NET-p surface span (y = 0, theta = 0) of the
    rounded doping equals the electrical channel length ``cfg["l_ch_el"]``.  Translation
    invariant, so only the dopings and the rounding length enter."""
    target = float(cfg["l_ch_el"])
    r = cfg.get("r_um", R_UM)
    a_s, n_j, n_s = cfg["na_s"], 1e17, 1e20
    x1 = 50.0  # far end of the source box (the right edge does not matter)

    def span(ch):
        f = lambda x: (  # noqa: E731
            a_s * smooth_step(x, r)
            - n_j * smooth_step(-x, r)
            - n_s * smooth_step(x - ch, r) * smooth_step(x1 - x, r)
        )
        g = np.linspace(-3 * r, ch + r, 4001)
        v = f(g)
        if v.max() <= 0.0:
            return 0.0  # no net-p surface at all
        k = int(np.argmax(v))
        xa = brentq(f, g[0], g[k], xtol=1e-12)
        xb = brentq(f, g[k], g[-1], xtol=1e-12)
        return xb - xa

    return float(brentq(lambda ch: span(ch) - target, target, target + 1.0, xtol=1e-9))


# ============================================================ mobility arms
class ArmMap(ChannelMap):
    """:class:`ChannelMap` whose column mobility follows the arm (the Picard driver is unchanged).

    ``noguchi``: stage 16.  ``const``: mu = 20 cm^2/Vs on the channel edges.  ``s07fit``: the
    piecewise-linear S07 Fig. 12 fit versus the gate field E_gate = q (N_dep + N_S) / eps_ox
    (equivalent to the spec conversion E_s = (eps_ox/eps_SiC) E_gate, N_S = eps_SiC E_s / q -
    N_dep); constant outside [0.5, 4] MV/cm.  The Noguchi diagnostics are kept in the result."""

    def __init__(self, *a, arm="noguchi", eps_ox_f_per_cm=None, **kw):
        super().__init__(*a, **kw)
        self.arm, self.eps_ox = arm, eps_ox_f_per_cm

    def column_mobility(self, p, n_scaled, clamp=True):
        res, ns_all = super().column_mobility(p, n_scaled, clamp=clamp)
        if self.arm == "noguchi":
            return res, ns_all
        res = dict(res)
        res["mu_noguchi"] = res["mu"]
        if self.arm == "const":
            res["mu"] = np.full_like(res["mu"], MU_CONST)
        elif self.arm == "s07fit":
            e_gate = Q * (res["n_dep"] + res["ns_used"]) / self.eps_ox * 1e-6  # MV/cm
            res["e_gate"] = e_gate
            res["mu"] = np.interp(e_gate, FIG12_X, FIG12_Y)
        else:
            raise ValueError(self.arm)
        return res, ns_all


def _lst(a):
    return np.asarray(a, dtype=float).tolist()


def registry_facts(sc):
    return {"eps_r_SiC": sc.eps_r, "eps_r_SiO2_registry": sc.reg.value("eps_r_SiO2_thermal"),
            "eps_ox_used": EPS_OX, "n_i_cm3": sc.n_i}  # fmt: skip


# ============================================================ jobs
def job_pde(args):
    """One channel-model on-state solve (R_pde, Picard record, share, profile) on ``level``."""
    key, cfg, level, arm = args
    t0 = time.perf_counter()
    rec = {"key": key, "kind": "pde", "level": level, "arm": arm, "cfg": cfg, "error": None}
    try:
        cfg = dict(cfg)
        design, meta = mdesign(cfg)
        rec["meta"] = meta
        sc = SiC4H(TracedRegistry.from_json())
        P = NoguchiParams.from_registry(sc.reg)
        ev = CorrectedCornerAdjoint(sc, design, level, r_um=cfg.get("r_um", R_UM), cuts_um=(),
                                    tau_v=5.0)  # fmt: skip
        ev.V_GS, ev.V_DS = V_GS, V_DS
        eps_ox = EPS_OX * EPS0 * 1e-2
        cmap = ArmMap(sc, ev.mod, ev.nd, ev.na, D_CH_UM, gate_x_max=ev.prm.gate_x1, arm=arm,
                      eps_ox_f_per_cm=eps_ox)  # fmt: skip
        cmap.mu_bulk = np.asarray(ev.mod.mu_edge, dtype=float).copy()
        sol_b, bc = ev.on_state()
        ev.mod.mu_edge = cmap.mu_bulk.copy()
        r = picard_solve(cmap, P, sol_b.X, bc)
        rec["picard"] = {k: r[k] for k in ("converged", "picard_converged", "newton_all_converged",
                                           "iterations", "final_rel", "newton_its")}  # fmt: skip
        rec["nodes"] = int(ev.mesh.n_nodes)
        rec["n_gated_columns"] = int(cmap.col_gated.sum())
        rec["n_channel_columns"] = int(len(cmap.cols))
        if not r["converged"]:
            raise EvalError("channel Picard/Newton did not converge")
        i_sc = ev.mod.contact_current(r["sol"], ev.base.dev.nodes["drain"])
        i_a = ev.base.current_a_per_cm(abs(i_sc))
        rec["R_pde"] = float(ev.V_DS * ev.x_cm / i_a)
        rec["x_max_cm"] = float(ev.x_cm)
        rec["l_c_cm"] = source_contact_length_cm(ev.base)
        vds_s = float(sc.volt_s(ev.V_DS))
        rec["share_channel"] = cmap.drop_share(r["sol"].phi_n, vds_s)
        c = r["columns"]
        gated = cmap.col_gated[cmap.cols]
        cols = cmap.cols
        rec["counts"] = r["counts"]
        rec["profile"] = {
            "x_um": [float(sc.length_cm(cmap.mesh.x[i])) * 1e4 for i in cols[gated]],
            "N_S_cm2": _lst(r["ns_all"][cols][gated]), "mu": _lst(c["mu"][gated]),
            "E_eff_MV_cm": _lst(c["e_eff"][gated]), "N_A_cm3": _lst(cmap.col_na_cm3[cols][gated]),
            "mu_noguchi": _lst(c.get("mu_noguchi", c["mu"])[gated]),
            "e_gate_MV_cm": _lst(c["e_gate"][gated]) if "e_gate" in c else None,
        }  # fmt: skip
        rec["ns_gated_max"] = float(np.max(r["ns_all"][cols][gated]))
        rec["ns_gated_min"] = float(np.min(r["ns_all"][cols][gated]))
        rec["mu_gated_mean"] = float(np.mean(c["mu"][gated]))
    except Exception as exc:
        rec["error"] = f"{type(exc).__name__}: {exc}"
    rec["seconds"] = time.perf_counter() - t0
    return rec


def job_off(args):
    """Punch-through (stage-19/20 method, V_pt,es) and dense avalanche BV of one design (level 0)."""
    key, cfg, level = args
    t0 = time.perf_counter()
    rec = {"key": key, "kind": "off", "level": level, "error": None}
    try:
        design, meta = mdesign(cfg)
        sc = SiC4H(TracedRegistry.from_json())
        ev = CorrectedCornerAdjoint(sc, design, level, r_um=cfg.get("r_um", R_UM), cuts_um=(),
                                    tau_v=5.0)  # fmt: skip
        ev.set_theta(np.zeros(5))
        cols = source_columns(ev)
        idx = column_index_matrix(ev, cols)
        v_ceil = 1500.0
        try:
            r = OffMarch(ev, np.zeros(5)).vpt_es(idx, 10.0, v_max_v=v_ceil)
            rec["V_pt_es"] = float(r["V_pt_es"]) if r["V_pt_es"] is not None else None
            rec["V_pt_ceiling"] = bool(r["V_pt_es"] is None)
            rec["worst_column_x_um"] = (
                None
                if r["V_pt_es"] is None
                else float(sc.length_cm(ev.mesh.x[cols[r["worst_column_row"]]])) * 1e4
            )
            rec["barrier_curve"] = [list(map(float, c[:2])) for c in r["curve"]]
        except Exception as exc:
            rec["V_pt_error"] = f"{type(exc).__name__}: {exc}"
        xs = np.linspace(0.0, meta["x_max"] - 0.1, 26)
        dc = DenseCuts(ev, xs)
        try:
            d = dense_bv(ev, np.zeros(5), dc, tol_v=0.5, loose_tol_v=3.0, cap=1.3)
            rec["bv_min"] = d["bv_min"]
            rec["x_argmin_um"] = d["x_argmin_um"]
            rec["bv_x_um"] = _lst(xs)
            rec["bv"] = _lst(d["bv"])
            rec["bv_resolved_all"] = bool(d["resolved_all"])
            rec["bv_n_unbracketed"] = int(d["n_unbracketed"])
        except Exception as exc:
            rec["bv_error"] = f"{type(exc).__name__}: {exc}"
        rec["stats"] = dict(ev.stats)
    except Exception as exc:
        rec["error"] = f"{type(exc).__name__}: {exc}"
    rec["seconds"] = time.perf_counter() - t0
    return rec


# ============================================================ (b) closed form
def comparison_b(sc, P, na_grid, gate_fields=(1.0, 2.0, 3.0, 4.0)):
    eps_sic = eps_sic_f_per_cm(sc)
    eps_ox = EPS_OX * EPS0 * 1e-2
    rows = []
    for na in na_grid:
        ndep = float(depletion_sheet(sc, na))
        pts = []
        for eg in gate_fields:
            e_s = eps_ox / eps_sic * eg  # MV/cm
            ns = eps_sic * e_s * 1e6 / Q - ndep
            fit = float(np.interp(eg, FIG12_X, FIG12_Y))
            pt = {"E_gate_MV_cm": eg, "E_s_MV_cm": e_s, "N_S_cm2": float(ns),
                  "mu_fit": fit, "included": bool(ns > 0.0)}  # fmt: skip
            if ns > 0.0:
                r = noguchi_mobility(sc, P, na, ns, clamp=False)
                pt.update({k: float(r[k]) for k in ("mu", "mu_ph", "mu_sr", "mu_c", "e_eff")})
                pt["ratio"] = pt["mu"] / fit
                pt["within_1p5"] = bool(1 / 1.5 <= pt["ratio"] <= 1.5)
                pt["ns_beyond_fit"] = bool(ns > N_S_FIT_MAX)
            pts.append(pt)
        inc = [p for p in pts if p["included"]]
        rows.append({"N_A": float(na), "N_dep_cm2": ndep, "points": pts,
                     "n_included": len(inc),
                     "pass": bool(inc and all(p["within_1p5"] for p in inc))})  # fmt: skip
    return rows


# ============================================================ figures
def plot_all(figdir: Path, res: dict, label: str):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figdir.mkdir(parents=True, exist_ok=True)

    def foot(fig):
        fig.text(0.01, 0.005, label, fontsize=6, color="0.4")

    ra = res["a"]
    fig, ax = plt.subplots(figsize=(7.2, 4.0), constrained_layout=True)
    nas = [r["N_A"] for r in ra["rows"]]
    for rho, c in ((5e-4, "C0"), (1e-3, "C1")):
        ys = [r[f"R_rho{rho:g}"] * 1e3 for r in ra["rows"]]
        ax.plot(range(len(nas)), ys, "o-", color=c, label=f"model, rho_c {rho:g}")
        ax.axhline(TABLE_I[(rho, "20")] * 1e3, color=c, ls="--", lw=1,
                   label=f"S07 Table I {TABLE_I[(rho, '20')] * 1e3:.2f}")  # fmt: skip
    ax.set_xticks(range(len(nas)), [f"{n:.1e}" for n in nas])
    ax.set_xlabel("uniform body N_A [cm^-3] (unknown in S07 simulation)")
    ax.set_ylabel("R_on,sp [mOhm cm^2]")
    ax.set_title("(a) mu = 20 cm^2/Vs, 8 um pitch")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)
    foot(fig)
    fig.savefig(figdir / "stage26_a_cross_simulator.png", dpi=130)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 4.0), constrained_layout=True)
    xs = np.linspace(0.5, 4.0, 50)
    ax.plot(xs, np.interp(xs, FIG12_X, FIG12_Y), "k-", lw=2, label="S07 Fig. 12 fit (figure-read)")
    for row in res["b"]["rows"]:
        if row["N_A"] not in res["b"]["na_plot"]:
            continue
        pts = [p for p in row["points"] if p["included"]]
        ax.plot([p["E_gate_MV_cm"] for p in pts], [p["mu"] for p in pts], "o-",
                label=f"Noguchi, N_A {row['N_A']:.1e}")  # fmt: skip
    ax.set_xlabel("gate (oxide) field [MV/cm]")
    ax.set_ylabel("mu_ch [cm^2/Vs]")
    ax.set_title("(b) closed-form Noguchi vs S07 Fig. 12 fit")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)
    foot(fig)
    fig.savefig(figdir / "stage26_b_mobility.png", dpi=130)
    plt.close(fig)

    rc = res["c"]
    fig, ax = plt.subplots(figsize=(7.2, 4.4), constrained_layout=True)
    oat = rc["oat"]
    names = [o["label"] for o in oat]
    lo = [o["delta_lo_mohm"] for o in oat]
    hi = [o["delta_hi_mohm"] for o in oat]
    y = np.arange(len(oat))
    ax.barh(y, lo, color="C0", alpha=0.8)
    ax.barh(y, hi, color="C3", alpha=0.8)
    ax.set_yticks(y, names, fontsize=7)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("change of R_on,sp vs central [mOhm cm^2] (Noguchi arm)")
    ax.set_title("(c) one-at-a-time sensitivities")
    ax.grid(alpha=0.3, axis="x")
    foot(fig)
    fig.savefig(figdir / "stage26_c_oat.png", dpi=130)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 4.0), constrained_layout=True)
    b = rc["bracket"]
    ax.axhspan(b["min"] * 1e3, b["max"] * 1e3, color="C0", alpha=0.2, label="Noguchi bracket")
    ax.axhline(R_MEAS * 1e3, color="k", lw=2, label="S07 measured 6.95")
    ax.axhline(rc["central"]["noguchi"]["R_total"] * 1e3, color="C0", label="Noguchi central")
    if rc["central"].get("s07fit"):
        ax.axhline(rc["central"]["s07fit"]["R_total"] * 1e3, color="C2", ls="--",
                   label="S07-fit arm central")  # fmt: skip
    ax.set_xticks([])
    ax.set_ylabel("R_on,sp [mOhm cm^2]")
    ax.set_title("(c) predicted vs measured R_on,sp")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)
    foot(fig)
    fig.savefig(figdir / "stage26_c_bracket.png", dpi=130)
    plt.close(fig)

    trend = rc.get("wj_trend", {})
    if trend:
        fig, ax = plt.subplots(figsize=(7.2, 4.0), constrained_layout=True)
        wj = sorted((k for k in trend if k != "increments_mohm"), key=float)
        for arm, c in (("noguchi", "C0"), ("s07fit", "C2")):
            ax.plot([float(w) for w in wj], [trend[w][arm]["R_total"] * 1e3 for w in wj], "o-",
                    color=c, label=arm)  # fmt: skip
        ax.plot([1.0, 1.5, 2.0], [R_MEAS_WJ[w] * 1e3 for w in (1.0, 1.5, 2.0)], "ks--",
                label="S07 measured")  # fmt: skip
        ax.set_xlabel("JFET width W_J [um]")
        ax.set_ylabel("R_on,sp [mOhm cm^2]")
        ax.set_title("(c) R_on vs W_J")
        ax.legend(fontsize=7)
        ax.grid(alpha=0.3)
        foot(fig)
        fig.savefig(figdir / "stage26_c_wj_trend.png", dpi=130)
        plt.close(fig)


# ============================================================ main
def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--workers", type=int, default=min(10, os.cpu_count() or 1))
    ap.add_argument("--dev", action="store_true", help="DEV run (refused under results/)")
    ap.add_argument(
        "--scan",
        action="store_true",
        help="also run the disclosed extra sensitivity items (not in the ranges)",
    )
    ap.add_argument("--no-off", action="store_true", help="skip the off-state job")
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


def nominal_cfg(cfg: dict) -> dict:
    """Same layout with ch_len = L_ch between the 50 % points (the literal nominal mapping)."""
    c = dict(cfg)
    c["ch_len"] = c.pop("l_ch_el")
    return c


def corner_cfg(vals: dict, name: str) -> dict:
    return central_cfg(name=name, na_s=vals["na_s"], na_b=vals["na_b"], w_j_eff=vals["w_j_eff"],
                       l_ch_el=vals["l_ch"])  # fmt: skip


def run(args, out: Path, log_path: Path) -> None:
    t_start = time.perf_counter()
    man = build_manifest(REPO, sys.argv, [])
    man["dev"] = bool(args.dev)
    (out / "manifest.json").write_text(json.dumps(man, indent=2) + "\n")
    sc = SiC4H(TracedRegistry.from_json())
    P = NoguchiParams.from_registry(sc.reg)
    res: dict = {"dev": bool(args.dev), "registry": registry_facts(sc)}

    def dump():
        (out / "results.json").write_text(json.dumps(res, indent=2, default=float) + "\n")

    # ------------------------------------------------------------------ jobs
    jobs: dict = {}  # key -> args

    def add(key, cfg, level, arm):
        jobs[key] = (key, cfg, level, arm)

    cen = dict(CENTRAL)
    # (a)
    na_a = (2.5e17, 5e17, 1e18)
    for na in na_a:
        add(f"a|{na:g}|L{LEVEL_MAIN}", sim_cfg(na), LEVEL_MAIN, "const")
    add(f"a|{5e17:g}|L{LEVEL_FINE}", sim_cfg(5e17), LEVEL_FINE, "const")
    # (c) central, both arms, h/2 and h/4
    for arm in ("noguchi", "s07fit"):
        for lv in (LEVEL_MAIN, LEVEL_FINE):
            add(f"c|central|{arm}|L{lv}", central_cfg(), lv, arm)
    # (c) PDE corners (Noguchi arm, h/2); the central is one of them
    corners = []
    for combo in itertools.product(*[RANGES[k] for k in PDE_KEYS]):
        vals = dict(zip(PDE_KEYS, combo, strict=True))
        key = "c|corner|" + "|".join(f"{k}={v:g}" for k, v in vals.items())
        corners.append((key, vals))
        if vals == {k: cen[k] for k in PDE_KEYS}:
            continue  # = the central run
        add(key, corner_cfg(vals, key), LEVEL_MAIN, "noguchi")
    # W_J trend (both arms, h/2); W_J = 1.0 is the central
    for wj in (1.5, 2.0):
        cfg = central_cfg(name=f"wj{wj}", w_j_eff=wj, x_max=0.5 * (9.0 + wj))
        for arm in ("noguchi", "s07fit"):
            add(f"c|wj|{wj:g}|{arm}", cfg, LEVEL_MAIN, arm)
    # disclosed extra scan (NOT part of the ranges)
    if args.scan:
        add("scan|nominal_l_ch", nominal_cfg(central_cfg()), LEVEL_MAIN, "noguchi")
        add("scan|r_um=0.05", central_cfg(r_um=0.05), LEVEL_MAIN, "noguchi")
        add("scan|r_um=0.15", central_cfg(r_um=0.15), LEVEL_MAIN, "noguchi")
        add("scan|t_csl=0.4", central_cfg(t_csl=0.4), LEVEL_MAIN, "noguchi")
        add("scan|na_uniform_5e17", central_cfg(na_b=None), LEVEL_MAIN, "noguchi")
    off_args = ("off|central", central_cfg(), 0)

    pool = ProcessPoolExecutor(args.workers, mp_context=mp.get_context("fork"))
    futs = {}
    # heavy jobs first
    order = sorted(jobs, key=lambda k: -(jobs[k][2]))
    for k in order:
        futs[pool.submit(job_pde, jobs[k])] = k
    if not args.no_off:
        futs[pool.submit(job_off, off_args)] = off_args[0]
    done: dict = {}
    for f in as_completed(futs):
        k = futs[f]
        try:
            done[k] = f.result()
        except Exception as exc:  # a worker crash is recorded, never hidden
            done[k] = {"key": k, "error": f"{type(exc).__name__}: {exc}", "kind": "crash"}
        d = done[k]
        msg = d.get("error") or f"R_pde {d.get('R_pde', float('nan')) * 1e3:.4f} mOhm cm2"
        print(f"{k}: {msg} ({d.get('seconds', 0):.0f} s)", flush=True)
    pool.shutdown()
    res["jobs_seconds"] = {k: d.get("seconds") for k, d in done.items()}
    res["raw"] = {k: {kk: vv for kk, vv in d.items() if kk != "barrier_curve"}
                  for k, d in done.items()}  # fmt: skip

    dump()  # raw job records survive an analysis failure

    def pde(key):
        d = done.get(key)
        if d is None or d.get("error") or "R_pde" not in d:
            raise RuntimeError(f"job {key} failed: {None if d is None else d.get('error')}")
        return d

    def series(d, rho_c, r_sub):
        dr = contact_resistance_sp(rho_c, d["x_max_cm"], d["l_c_cm"])
        return {"R_pde": d["R_pde"], "dR_c": dr, "R_sub": r_sub, "R_total": d["R_pde"] + dr + r_sub}

    # ------------------------------------------------------------------ (a)
    rows = []
    for na in na_a:
        d = pde(f"a|{na:g}|L{LEVEL_MAIN}")
        row = {
            "N_A": na,
            "R_pde": d["R_pde"],
            "l_c_cm": d["l_c_cm"],
            "x_max_cm": d["x_max_cm"],
            "share_channel": d["share_channel"],
            "picard": d["picard"],
            "nodes": d["nodes"],
        }
        for rho in (5e-4, 1e-3):
            row[f"R_rho{rho:g}"] = series(d, rho, 0.0)["R_total"]
        rows.append(row)
    d4 = pde(f"a|{5e17:g}|L{LEVEL_FINE}")
    lo = {rho: min(r[f"R_rho{rho:g}"] for r in rows) for rho in (5e-4, 1e-3)}
    hi = {rho: max(r[f"R_rho{rho:g}"] for r in rows) for rho in (5e-4, 1e-3)}
    cen_a = rows[1]
    crit_a = {}
    for rho in (5e-4, 1e-3):
        t = TABLE_I[(rho, "20")]
        crit_a[f"rho{rho:g}"] = {
            "table_I": t,
            "bracket": [lo[rho], hi[rho]],
            "central": cen_a[f"R_rho{rho:g}"],
            "in_bracket": bool(lo[rho] <= t <= hi[rho]),
            "central_rel_dev": cen_a[f"R_rho{rho:g}"] / t - 1.0,
            "central_within_20pct": bool(abs(cen_a[f"R_rho{rho:g}"] / t - 1.0) <= 0.20),
        }
    implied = (TABLE_I[(1e-3, "20")] - TABLE_I[(5e-4, "20")]) / 5e-4  # x_max/L_c implied by Table I
    res["a"] = {"rows": rows, "criterion_A": crit_a,
                "criterion_A_pass": bool(all(v["in_bracket"] and v["central_within_20pct"]
                                             for v in crit_a.values())),
                "h4_central_R_pde": d4["R_pde"],
                "h4_vs_h2_pde_rel": d4["R_pde"] / rows[1]["R_pde"] - 1.0,
                "table_I_implied_xmax_over_Lc": implied,
                "model_xmax_over_Lc": rows[1]["x_max_cm"] / rows[1]["l_c_cm"]}  # fmt: skip

    # ------------------------------------------------------------------ (b)
    na_grid = np.linspace(2.5e17, 5e17, 11)
    rows_b = comparison_b(sc, P, na_grid)
    ends = [r for r in rows_b if r["N_A"] in (na_grid[0], na_grid[-1])]
    res["b"] = {"rows": rows_b, "na_plot": [float(na_grid[0]), float(na_grid[-1])],
                "ends": ends, "P26b_pass": bool(any(r["pass"] for r in rows_b)),
                "n_na_passing": int(sum(r["pass"] for r in rows_b)),
                "excluded_points": [{"N_A": r["N_A"], "E_gate": p["E_gate_MV_cm"]}
                                    for r in rows_b for p in r["points"] if not p["included"]]}  # fmt: skip

    # ------------------------------------------------------------------ (c)
    def total(d, vals):
        return series(d, vals["rho_c"], vals["r_sub"])

    d_cen = {arm: {lv: pde(f"c|central|{arm}|L{lv}") for lv in (LEVEL_MAIN, LEVEL_FINE)}
             for arm in ("noguchi", "s07fit")}  # fmt: skip
    central = {}
    for arm in ("noguchi", "s07fit"):
        s1 = total(d_cen[arm][LEVEL_MAIN], cen)
        s2 = total(d_cen[arm][LEVEL_FINE], cen)
        central[arm] = {
            **s1, "R_total_h4": s2["R_total"], "h4_vs_h2_total_rel": s2["R_total"] / s1["R_total"] - 1.0,
            "h4_vs_h2_pde_rel": s2["R_pde"] / s1["R_pde"] - 1.0,
            "ratio_to_measured": s1["R_total"] / R_MEAS, "share_channel": d_cen[arm][LEVEL_MAIN]["share_channel"],
            "picard": d_cen[arm][LEVEL_MAIN]["picard"], "nodes_h2": d_cen[arm][LEVEL_MAIN]["nodes"],
            "nodes_h4": d_cen[arm][LEVEL_FINE]["nodes"], "counts": d_cen[arm][LEVEL_MAIN]["counts"],
            "profile": d_cen[arm][LEVEL_MAIN]["profile"],
            "mu_gated_mean": d_cen[arm][LEVEL_MAIN]["mu_gated_mean"],
            "ns_gated_range": [d_cen[arm][LEVEL_MAIN]["ns_gated_min"], d_cen[arm][LEVEL_MAIN]["ns_gated_max"]],
            "n_gated_columns_h2": d_cen[arm][LEVEL_MAIN]["n_gated_columns"],
            "n_gated_columns_h4": d_cen[arm][LEVEL_FINE]["n_gated_columns"],
        }  # fmt: skip
    # corner table (Noguchi arm)
    corner_rows = []
    for key, vals in corners:
        d = d_cen["noguchi"][LEVEL_MAIN] if vals == {k: cen[k] for k in PDE_KEYS} else pde(key)
        for rho, rs in itertools.product(RANGES["rho_c"], RANGES["r_sub"]):
            s = series(d, rho, rs)
            corner_rows.append({**vals, "rho_c": rho, "r_sub": rs, **s})
    rt = [r["R_total"] for r in corner_rows]
    bracket = {"min": min(rt), "max": max(rt), "width_ratio": max(rt) / min(rt),
               "n_corners": len(rt), "contains_measured": bool(min(rt) <= R_MEAS <= max(rt)),
               "weak_label": bool(max(rt) / min(rt) >= 2.0 - 1e-12)}  # fmt: skip
    # one-at-a-time around the central (PDE keys from the corner runs, series analytically)
    r_c = central["noguchi"]["R_total"]
    oat = []
    for k in PDE_KEYS + SERIES_KEYS:
        ends_k = {}
        for j, end in enumerate(RANGES[k]):
            v = dict(cen)
            v[k] = end
            if k in PDE_KEYS:
                ck = "c|corner|" + "|".join(f"{kk}={v[kk]:g}" for kk in PDE_KEYS)
                d = (
                    d_cen["noguchi"][LEVEL_MAIN]
                    if all(v[x] == cen[x] for x in PDE_KEYS)
                    else pde(ck)
                )
            else:
                d = d_cen["noguchi"][LEVEL_MAIN]
            ends_k[j] = total(d, v)["R_total"] - r_c
        oat.append({"param": k, "label": f"{k} {RANGES[k][0]:g}..{RANGES[k][1]:g}",
                    "central": cen[k], "range": list(RANGES[k]),
                    "delta_lo_mohm": ends_k[0] * 1e3, "delta_hi_mohm": ends_k[1] * 1e3,
                    "span_mohm": abs(ends_k[1] - ends_k[0]) * 1e3})  # fmt: skip
    oat.sort(key=lambda o: -o["span_mohm"])
    trend = {"1": {arm: central[arm] for arm in ("noguchi", "s07fit")}}
    for wj in (1.5, 2.0):
        trend[f"{wj:g}"] = {}
        for arm in ("noguchi", "s07fit"):
            d = pde(f"c|wj|{wj:g}|{arm}")
            trend[f"{wj:g}"][arm] = {
                **series(d, cen["rho_c"], cen["r_sub"]),
                "share_channel": d["share_channel"],
            }
    # structure of the trend: increase 1 -> 1.5 -> 2 versus measured +0.55, +0.95 mOhm cm^2
    for arm in ("noguchi", "s07fit"):
        r1, r15, r2 = (trend[w][arm]["R_total"] for w in ("1", "1.5", "2"))
        trend.setdefault("increments_mohm", {})[arm] = {
            "d_1_to_1p5": (r15 - r1) * 1e3,
            "d_1_to_2": (r2 - r1) * 1e3,
        }
    trend["increments_mohm"]["measured"] = {"d_1_to_1p5": 0.55, "d_1_to_2": 0.95}
    scan = {}
    if args.scan:
        for key in ("nominal_l_ch", "r_um=0.05", "r_um=0.15", "t_csl=0.4", "na_uniform_5e17"):
            d = done.get(f"scan|{key}")
            if d is None or d.get("error"):
                scan[key] = {"error": None if d is None else d.get("error")}
                continue
            s = series(d, cen["rho_c"], cen["r_sub"])
            scan[key] = {
                **s,
                "delta_vs_central_mohm": (s["R_total"] - r_c) * 1e3,
                "n_gated_columns": d["n_gated_columns"],
                "meta": d["meta"],
            }
    p26c = {"measured": R_MEAS, "noguchi_central": r_c, "noguchi_ratio": r_c / R_MEAS,
            "central_rel_dev": r_c / R_MEAS - 1.0,
            "central_within_20pct": bool(abs(r_c / R_MEAS - 1.0) <= 0.20),
            "bracket_contains_measured": bracket["contains_measured"],
            "P26c_pass": bool(bracket["contains_measured"] and abs(r_c / R_MEAS - 1.0) <= 0.20),
            "bracket_width_ratio": bracket["width_ratio"], "weak": bracket["weak_label"]}  # fmt: skip
    off = done.get("off|central")
    res["c"] = {"central": central, "bracket": bracket, "corners": corner_rows, "oat": oat,
                "oat_ranking": [o["param"] for o in oat], "wj_trend": trend, "scan": scan,
                "P26c": p26c, "measured_wj": {str(k): v for k, v in R_MEAS_WJ.items()},
                "off_state": None if off is None else {k: v for k, v in off.items()
                                                      if k not in ("barrier_curve",)},
                "bv_s06_sim": BV_S06_SIM}  # fmt: skip

    # ------------------------------------------------------------------ criteria
    errs = [k for k, d in done.items() if d.get("error")]
    pde_keys = [k for k, d in done.items() if d.get("kind") == "pde"]
    blas = [ln for ln in log_path.read_text(errors="replace").splitlines() if "On entry to" in ln]
    off_ok = args.no_off or (off is not None and not off.get("error"))
    c1 = bool(
        not errs
        and not blas
        and off_ok
        and all(done[k]["picard"]["newton_all_converged"] for k in pde_keys)
        and (args.no_off or ("V_pt_error" not in off and "bv_error" not in off))
    )
    c2_vals = {arm: central[arm]["h4_vs_h2_total_rel"] for arm in central}
    c2 = bool(all(abs(v) <= C2_TOL for v in c2_vals.values()))
    c3 = bool(
        all(
            done[k]["picard"]["picard_converged"] and done[k]["picard"]["final_rel"] < 1e-10
            for k in pde_keys
        )
    )
    crit = {"C1": c1, "C2": c2, "C3": c3}
    vals = {
        "C1": {"job_errors": errs, "on_entry_to_lines": len(blas), "off_state_ok": off_ok},
        "C2": {
            "h4_vs_h2_total_rel": c2_vals,
            "h4_vs_h2_pde_rel": {a: central[a]["h4_vs_h2_pde_rel"] for a in central},
            "tol": C2_TOL,
        },
        "C3": {
            "max_iterations": max(done[k]["picard"]["iterations"] for k in pde_keys),
            "max_final_rel": max(done[k]["picard"]["final_rel"] for k in pde_keys),
        },
    }
    res["interpretations"] = interpretations()
    res["criteria"], res["criteria_values"] = crit, vals
    res["runtime_seconds"] = time.perf_counter() - t_start
    dump()
    figdir = (REPO / "figures") if out == DEFAULT_OUT.resolve() else out / "figures"
    label = f"stage 26 | commit {str(man.get('git_commit'))[:8]} | {man['timestamp_utc'][:19]}Z"
    try:
        plot_all(figdir, res, label)
    except Exception as exc:
        res["figures_error"] = f"{type(exc).__name__}: {exc}"
    summary = {
        "criteria": crit, "gate_passed": bool(all(crit.values())), "dev": bool(args.dev),
        "criteria_values": vals,
        "R1": {"a": {"rows": [{k: r[k] for k in ("N_A", "R_pde", "R_rho0.0005", "R_rho0.001")}
                              for r in rows], "criterion_A": crit_a,
                     "criterion_A_pass": res["a"]["criterion_A_pass"],
                     "table_I_implied_xmax_over_Lc": implied,
                     "model_xmax_over_Lc": res["a"]["model_xmax_over_Lc"]},
               "b": {"P26b_pass": res["b"]["P26b_pass"], "ends": ends,
                     "n_na_passing": res["b"]["n_na_passing"]},
               "c": {"P26c": p26c, "bracket": bracket,
                     "central": {a: {k: v for k, v in central[a].items() if k != "profile"}
                                 for a in central},
                     "oat_ranking": res["c"]["oat_ranking"], "oat": oat,
                     "wj_trend": trend, "scan": scan,
                     "off_state": None if off is None else {
                         k: off.get(k) for k in ("V_pt_es", "V_pt_ceiling", "bv_min", "x_argmin_um",
                                                 "bv_resolved_all", "V_pt_error", "bv_error")}}},
        "runtime_seconds": res["runtime_seconds"],
    }  # fmt: skip
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    print(json.dumps(crit, indent=2))
    print(f"runtime {res['runtime_seconds']:.1f} s")


def interpretations() -> dict:
    return {
        "half_cell": "x in [0, pitch/2], x=0 JFET centre; x_jfet = W_J/2 (W_J = full width)",
        "layout": "gate end x_g = x_cs - 0.5 um (gate-source separation); source contact "
        "[x_cs, x_max]; x_cs = x_max - L_contact (central 3 um, comparison (a) LS 2 um)",
        "gate_over_source": "gate overhangs the n+ source by L_ov = x_g - x_src0 (x_src0 = "
        "x_jfet + ch_len, ch_len calibrated to the electrical L_ch)",
        "channel_length": "L_ch is the electrical (net-p surface) channel length of the rounded "
        "doping (r = 0.1 um, as in all earlier stages): ch_len between the 50 % points is "
        "calibrated per configuration. The literal mapping ch_len = L_ch is degenerate (the n+ "
        "tail removes the whole net-p surface at N_A,s = 2.5e17: 5 of 16 pre-run corners had no "
        "channel column) and is reported as scan item nominal_l_ch. This choice was made after "
        "the first DEV pre-run exposed that degeneracy.",
        "p_contact": "no separate p+ region is modelled: a 0.06 um body-contact strip at the cell "
        "edge (bc_w 0.04, pad_w 0.02); modelled source contact = L_contact - 0.04 um",
        "L_c": "Delta R_c uses the modelled n+ source contact length (central 2.96 um, (a) 1.96 um)",
        "comparison_a_assumptions": "pitch 8 um with LS = 2 um read as source-contact length; gate "
        "geometry (0.5 um separation) copied from the fabricated cell; "
        "t_ox 50 nm; V_G 20 V; V_DS 0.1 V; CSL 0.3 um at 1e17; "
        "JFET 1e17; uniform body N_A; no substrate",
        "csl_thickness": "t_csl = 0.3 um below the 0.6 um body (Fig. 1 depth); JFET column 1e17 "
        "over the body depth; drift 6 um at 1e16 below the CSL",
        "ionization_levels": "acceptor ionization levels, hole mobility and V_FB at N_A,s; the "
        "retrograde N_A(y) enters only as doping (documented approximation)",
        "gate_material": "V_FB of the existing model (n+ poly-Si gate work function); S07 gate "
        "material not used",
        "bias": "V_GS = 20 V (top-curve label of Fig. 9); V_DS = 0.1 V (stage convention; S07 bias "
        "not stated); T = 300 K",
        "R_sub_central": "R_sub central = 0.35 mOhm cm^2 (range midpoint; the spec gives only "
        "the range [0, 0.7])",
        "w_j_eff": "W_J,eff changes x_jfet only; the drawn layout (x_g, x_cs, pitch) is fixed",
        "wj_trend": "pitch = 9 + W_J (all other lengths fixed)",
        "bracket": "the 16 PDE corners (N_A,s, N_A,b, W_J,eff, L_ch) x 4 series corners (rho_c, "
        "R_sub) = 64 totals; the one-at-a-time points are corners of this set",
        "S07_arm": "E_gate = q(N_dep + N_S)/eps_ox from the column state (N_dep at the column "
        "surface N_A), piecewise-linear Fig. 12 fit, constant outside 0.5..4 MV/cm; "
        "figure-read and approximate",
        "mu_const_scope": "constant 20 cm^2/Vs on every channel edge (both nodes within d_ch of "
        "the interface, p-type surface columns); accumulation/bulk unchanged",
        "off_state": "V_pt,es: smallest V_D with min n+ -> body barrier < 10 V_T (stage 20 method, "
        "level h); BV: dense avalanche cuts (26 columns) on level h, zero gate "
        "bias, no termination",
        "C2_definition": "|R_total(h/4)/R_total(h/2) - 1| <= 3 % for both arms (the PDE-only "
        "change is also reported)",
    }


if __name__ == "__main__":
    main()
