"""Stage-9b benchmark (docs/benchmarks/stage9b_avalanche.md): avalanche breakdown, 1D.

Criteria B1-B5 (gating) and R1 (reported). Writes <out>/{results,summary,manifest}.json.

Study parameters (scaled units, alpha in 1/L_D, E in V_T/L_D; NOT material values):
* DIODE_CH  = (a_n, b_n, a_p, b_p) = (6, 30, 3, 39): p+n diodes (B2-B4, R1); BV_av ~ 66.6 V_T
  on the stage-4 p+n structure, inside the verified stage-4 bias range.
* B1: low-field set with b_n = 64, a_n solved so the uniform-field critical field is
  E* = 4 V_T/L_D for an i-layer of W = 30 (a_p = a_n/2, b_p = 1.3 b_n), contacts of 2 N_0,
  delta = 1e-3. Why a separate set and light contacts: the contact depletion layers and the
  hole/electron accumulation overshoot at the i-layer edges ionize too; their relative
  weight scales like 1/(E* W) (a design scan over contact doping, W and b/E* preceded this
  choice; mismatch fell from ~1e-2 to ~5e-4), so a steep alpha
  (b/E* = 16, realistic for SiC) and light contacts are chosen to push the edge-tail mismatch
  below the 1e-3 gate while keeping BV in the stage-4 range. The values are recorded in
  results.json.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from widegap.fv.avalanche import (
    AvalancheQFDD1D,
    Chynoweth,
    avalanche_bv,
    injected_reverse_solutions,
    ionization_integrals,
    uniform_field_critical,
    uniform_field_integral,
)
from widegap.fv.drift_diffusion import DDParams
from widegap.fv.equilibrium import cv_average_piecewise
from widegap.fv.qf_drift_diffusion import pack
from widegap.fv.reverse_bias import Diode1D, breakdown_voltage, sweep
from widegap.manifest import build_manifest

REPO = Path(__file__).resolve().parents[2]
DELTA, ECRIT = 1e-3, 15.0
DIODE_CH = Chynoweth(6.0, 30.0, 3.0, 39.0)
# stage-4 structures (benchmarks/stage4/run.py CASES), copied to stay independent of it
S4_CASES = {
    "p+n": dict(lo=-5.0, hi=60.0, breaks=[0.0], vals=[-100.0, 1.0]),
    "p+nn+": dict(lo=-5.0, hi=13.0, breaks=[0.0, 8.0], vals=[-100.0, 1.0, 100.0]),
}
B1_CFG = dict(W=30.0, contact=2.0, pad=10.0, e_target=4.0, b_n=64.0, meshes=[0.025, 0.0125])
DEFAULT_OUT = REPO / "results" / "stage9b"


def s4_diode(case, h):
    c = S4_CASES[case]
    x = np.linspace(c["lo"], c["hi"], int(round((c["hi"] - c["lo"]) / h)) + 1)
    return Diode1D(x, cv_average_piecewise(x, c["breaks"], c["vals"]), DELTA)


def rb_summary(reports):
    """Stage-4 solver bookkeeping. The stage-4 'converged' flag additionally demands a Newton
    step <= 1e-8, which is round-off limited for |psi| >~ 70 on long meshes; the residual
    criterion (final residual <= its own round-off tolerance) is reported separately."""
    return {
        "n_solves": len(reports),
        "flag_converged": int(sum(r.converged for r in reports)),
        "residual_converged": int(sum(r.residual_history[-1] <= r.tolerance for r in reports)),
        "max_iterations": int(max(r.iterations for r in reports)),
    }


def b1_chynoweth():
    w, e_t, b = B1_CFG["W"], B1_CFG["e_target"], B1_CFG["b_n"]

    def f(a):
        ch = Chynoweth(a, b, a / 2.0, 1.3 * b)
        return uniform_field_critical(ch, w, e_lo=0.2) - e_t

    a = float(brentq(f, 1.5 / w, 1e12, xtol=1e-14, rtol=1e-14))
    return Chynoweth(a, b, a / 2.0, 1.3 * b)


def b1_mesh(h):
    w, pad = B1_CFG["W"], B1_CFG["pad"]
    x = np.linspace(-pad, w + pad, int(round((w + 2 * pad) / h)) + 1)
    c = cv_average_piecewise(x, [0.0, w], [-B1_CFG["contact"], 0.0, B1_CFG["contact"]])
    return Diode1D(x, c, DELTA)


def run_b1(h, ch):
    w = B1_CFG["W"]
    d = b1_mesh(h)
    e_star = uniform_field_critical(ch, w, e_lo=0.2)
    xm = 0.5 * (d.x[1:] + d.x[:-1])
    margin = 20.0 / e_star  # hole/electron Boltzmann tails ~ exp(-E x): e^-20 below
    win = (xm > margin) & (xm < w - margin)
    a_star = ch.alpha_n(np.array([e_star]))[0]
    b_star = ch.alpha_p(np.array([e_star]))[0]

    def e_i(psi):
        return float(np.median(np.abs(d.element_field(psi))[win]))

    # reference voltage: i-layer field of the discrete electrostatics equals E*
    psi_lo, v_lo = d.solve(0.0)[0], 0.0
    while True:
        v_hi = v_lo + 2.0
        psi, _ = d.solve(v_hi, psi0=psi_lo)
        if e_i(psi) > e_star:
            break
        psi_lo, v_lo = psi, v_hi
    v_ref = float(
        brentq(lambda v: e_i(d.solve(v, psi0=psi_lo)[0]) - e_star, v_lo, v_hi, xtol=1e-11)
    )
    av = avalanche_bv(d, ch, "n", v_step=2.0)
    psi_bv = av.psi
    e_el = np.abs(d.element_field(psi_bv))
    e_in = e_el[win]
    uniformity_bv = float(np.max(np.abs(e_in / np.median(e_in) - 1.0)))
    psi_ref = d.solve(v_ref, psi0=psi_lo)[0]
    e_ref = np.abs(d.element_field(psi_ref))[win]
    uniformity_ref = float(np.max(np.abs(e_ref / np.median(e_ref) - 1.0)))
    # quadrature check: window-only field, exact closed form over the window width
    e_win = np.where(win, e_el, 0.0)
    i_win, _, _ = ionization_integrals(d.x, e_win, ch)
    w_win = float(np.sum(np.diff(d.x)[win]))
    i_closed = uniform_field_integral(
        ch.alpha_n(np.array([np.median(e_in)]))[0],
        ch.alpha_p(np.array([np.median(e_in)]))[0],
        w_win,
    )
    rel = av.bv / v_ref - 1.0
    return {
        "h": h,
        "n_nodes": len(d.x),
        "e_star": e_star,
        "alpha_n_star": float(a_star),
        "alpha_p_star": float(b_star),
        "bv_reference": v_ref,
        "bv_numerical": av.bv,
        "rel_mismatch": rel,
        "uniformity_window": [float(xm[win][0]), float(xm[win][-1])],
        "uniformity_rel_at_bv": uniformity_bv,
        "uniformity_rel_at_vref": uniformity_ref,
        "i_layer_field_at_bv": float(np.median(e_in)),
        "quadrature_window_vs_closed_form": abs(i_win - i_closed) / i_closed,
        "tail_ionization_excess_at_bv": 1.0 - i_win,  # I(full) = 1 at BV; I(window) = i_win
        "newton": rb_summary(av.reports),
    }


def discrete_identity(d, ch, psi):
    """Relative violation of 1 - I_n = exp(-A_W) (1 - I_p) on a computed field."""
    i_n, i_p, a_w = ionization_integrals(d.x, d.element_field(psi), ch)
    lhs, rhs = 1.0 - i_n, np.exp(-a_w) * (1.0 - i_p)
    return float(abs(lhs - rhs) / max(abs(lhs), abs(rhs), 1e-300)), i_n, i_p, a_w


def run_b2(ch_b1):
    out = {}
    cases = {
        "p+n_h0.05": (s4_diode("p+n", 0.05), DIODE_CH),
        "p+nn+_h0.05": (s4_diode("p+nn+", 0.05), DIODE_CH),
        "pin_h0.025": (b1_mesh(0.025), ch_b1),
    }
    for name, (d, ch) in cases.items():
        rn = avalanche_bv(d, ch, "n")
        rp = avalanche_bv(d, ch, "p")
        # identity at a bias where 1 - I_n = 0.5 (at the root it is a 0/0 comparison)
        ident, i_n, i_p, a_w = discrete_identity(d, ch, avalanche_bv(d, ch, "n", level=0.5).psi)
        out[name] = {
            "bv_n": rn.bv,
            "bv_p": rp.bv,
            "rel_diff": abs(rn.bv - rp.bv) / rn.bv,
            "integral_n_at_bv_n": rn.integral,
            "integral_p_at_bv_p": rp.integral,
            "identity_violation_at_I_n_0.5": ident,
            "I_n_I_p_A_W_at_I_n_0.5": [i_n, i_p, a_w],
        }
    return out


def run_b3(meshes):
    bvs, reps = [], []
    for h in meshes:
        r = avalanche_bv(s4_diode("p+n", h), DIODE_CH)
        bvs.append(r.bv)
        reps += r.reports
    d12, d23 = bvs[0] - bvs[1], bvs[1] - bvs[2]
    # round-off floor of a BV difference: 1e3 eps * BV
    floor = 1e3 * np.finfo(float).eps * bvs[-1]
    order = float(np.log2(abs(d12) / abs(d23))) if min(abs(d12), abs(d23)) > floor else float("nan")
    bv_rich = bvs[2] + (bvs[2] - bvs[1]) / (2.0**order - 1.0)
    return {
        "meshes": meshes,
        "bv": bvs,
        "differences": [d12, d23],
        "roundoff_floor": floor,
        "observed_order": order,
        "bv_richardson": float(bv_rich),
        "finest_vs_richardson": float(abs(bvs[2] - bv_rich) / bv_rich),
        "newton": rb_summary(reps),
    }


def run_b4(h, m_targets):
    d = s4_diode("p+n", h)
    rows = []
    for m_t in m_targets:
        rb = avalanche_bv(d, DIODE_CH, "n", level=1.0 - 1.0 / m_t)
        m0, s0, m1, s1 = injected_reverse_solutions(d, DIODE_CH, rb.bv)
        i_off = -(s0.jn + s0.jp)
        i_on = -(s1.jn + s1.jp)
        i_n, i_p, _ = m1.multiplication_integral(s1)
        m_dd = float(i_on[-1] / i_off[-1])
        m_form = float(1.0 / (1.0 - i_n))
        rows.append(
            {
                "m_target_on_zero_current_field": m_t,
                "v_reverse": rb.bv,
                "M_dd": m_dd,
                "M_formula_on_computed_field": m_form,
                "I_n_computed_field": i_n,
                "I_p_computed_field": i_p,
                "rel_error": abs(m_dd - m_form) / m_form,
                "total_flux_spread_G_on": float(i_on.max() / i_on.min() - 1.0),
                "injected_current": float(i_off[-1]),
                "hole_saturation_fraction": float(-s0.jp[-1] / i_off[-1]),
                "sign_violation": m1.sign_violation(s1.X),
                "newton_iterations": [s0.report.iterations, s1.report.iterations],
            }
        )
    gated = [r for r in rows if r["M_dd"] <= 10.0]
    return {
        "h": h,
        "rows": rows,
        "max_rel_error_M_le_10": max(r["rel_error"] for r in gated) if gated else float("nan"),
    }


def fd_check(dd, X, bc, rng, gen_only=False, n_vec=4, e=1e-6):
    """Worst per-row |J v - FD| / (|J||v|), for the full residual or the generation part alone
    (``generation_residual`` is evaluated directly, so no cancellation against the base terms)."""
    import scipy.sparse as sp

    if gen_only:
        jac = dd.generation_jacobian(X)

        def fun(y):
            return dd.generation_residual(y)
    else:
        jac = dd.jacobian(X)

        def fun(y):
            return dd.residual(y, bc)

    worst = 0.0
    for _ in range(n_vec):
        v = rng.standard_normal(len(X))
        v[[0, 1, 2, -3, -2, -1]] = 0.0
        fd = (fun(X + e * v) - fun(X - e * v)) / (2 * e)
        jv = jac @ v
        den = sp.csr_matrix(abs(jac)) @ np.abs(v)
        # Rows whose generation derivative is < 1e-12 of the largest generation row (alpha
        # ~ exp(-b/E) underflowing in the neutral bulk, den ~ 1e-239) are below the FD validity
        # floor: the perturbation exceeds the scale on which alpha varies. The full-residual
        # check covers those rows (generation is negligible against the base terms there).
        m = den > (1e-12 * den.max() if gen_only else 0.0)
        worst = max(worst, float(np.max(np.abs(jv - fd)[m] / den[m])))
    return worst


def run_b5(ch_strong):
    rng = np.random.default_rng(5)
    out = {"random_states": [], "converged_state": {}}
    for delta, rec in ((1e-3, True), (1e-3, False), (1e-24, True), (1e-24, False)):
        x = np.linspace(-3.0, 3.0, 25)
        c = cv_average_piecewise(x, [0.0], [-2.0, 1.0])
        prm = DDParams(delta=delta, mu_n=1.3, mu_p=0.7, tau_n=2.0, tau_p=3.0, recombination=rec)
        s = {
            k: (0.1 if delta > 1e-6 else 0.0) * rng.standard_normal(len(x))
            for k in ("psi", "n", "p")
        }
        dd = AvalancheQFDD1D(x, c, prm, np.log(delta), ch_strong, s)
        X = pack(rng.uniform(-3, 3, len(x)), rng.uniform(-1, 1, len(x)), rng.uniform(-1, 1, len(x)))
        bc = [tuple(X[:3]), tuple(X[-3:])]
        out["random_states"].append(
            {
                "delta": delta,
                "recombination": rec,
                "full_row_error": fd_check(dd, X, bc, rng),
                "generation_part_row_error": fd_check(dd, X, bc, rng, gen_only=True),
            }
        )
    d = s4_diode("p+n", 0.1)
    rb = avalanche_bv(d, DIODE_CH, "n", level=0.9)
    m0, s0, m1, s1 = injected_reverse_solutions(d, DIODE_CH, rb.bv)
    bc = [tuple(s1.X[:3]), tuple(s1.X[-3:])]
    out["converged_state"] = {
        "full_row_error": fd_check(m1, s1.X, bc, rng),
        "generation_part_row_error": fd_check(m1, s1.X, bc, rng, gen_only=True),
    }
    vals = [v for r in out["random_states"] for k, v in r.items() if k.endswith("row_error")]
    vals += list(out["converged_state"].values())
    out["max_row_error"] = float(max(vals))
    return out


def run_r1(h):
    out = {}
    for case in S4_CASES:
        d = s4_diode(case, h)
        sw = sweep(d, list(range(0, 121)))
        bv_proxy = breakdown_voltage(d, sw, ECRIT)[0]
        rec = {"h": h, "bv_proxy_ecrit15": float(bv_proxy)}
        try:
            r = avalanche_bv(d, DIODE_CH, "n", v_max=200.0)
            rec.update(bv_av=r.bv, bv_av_over_proxy=r.bv / bv_proxy)
        except (ValueError, RuntimeError) as exc:
            rec.update(bv_av=None, note=f"BV_av not obtained: {exc}")
        out[case] = rec
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    (out / "manifest.json").write_text(
        json.dumps(build_manifest(REPO, sys.argv, [5]), indent=2) + "\n"
    )

    ch1 = b1_chynoweth()
    res = {
        "parameters": {
            "diode_chynoweth": vars(DIODE_CH),
            "b1_chynoweth": vars(ch1),
            "b1_config": B1_CFG,
            "delta": DELTA,
            "ecrit_proxy": ECRIT,
        }
    }
    res["B1"] = [run_b1(h, ch1) for h in B1_CFG["meshes"]]
    res["B2"] = run_b2(ch1)
    res["B3"] = run_b3([0.1, 0.05, 0.025])
    res["B4"] = [run_b4(h, [2.0, 5.0, 10.0]) for h in (0.05, 0.025)]
    res["B5"] = run_b5(Chynoweth(3.0, 4.0, 1.5, 5.0))
    res["R1"] = run_r1(0.025)

    b1_ok = all(
        abs(r["rel_mismatch"]) <= 1e-3 and r["uniformity_rel_at_vref"] <= 1e-6 for r in res["B1"]
    )
    b2_ok = all(v["rel_diff"] <= 1e-8 for v in res["B2"].values())
    b3 = res["B3"]
    b3_ok = bool(0.8 <= b3["observed_order"] <= 2.5 and b3["finest_vs_richardson"] <= 0.01)
    b4_ok = all(r["max_rel_error_M_le_10"] <= 1e-2 for r in res["B4"])
    b5_ok = res["B5"]["max_row_error"] <= 1e-6
    crit = {"B1": b1_ok, "B2": b2_ok, "B3": b3_ok, "B4": b4_ok, "B5": b5_ok}
    crit = {k: bool(v) for k, v in crit.items()}
    summary = {
        "criteria": crit,
        "gate_passed": all(crit.values()),
        "values": {
            "B1_rel_mismatch": [r["rel_mismatch"] for r in res["B1"]],
            "B1_uniformity_at_vref": [r["uniformity_rel_at_vref"] for r in res["B1"]],
            "B2_rel_diff": {k: v["rel_diff"] for k, v in res["B2"].items()},
            "B3_order": b3["observed_order"],
            "B3_finest_vs_richardson": b3["finest_vs_richardson"],
            "B4_max_rel_error": [r["max_rel_error_M_le_10"] for r in res["B4"]],
            "B5_max_row_error": res["B5"]["max_row_error"],
        },
        "R1_reported": res["R1"],
    }
    (out / "results.json").write_text(json.dumps(res, indent=2, default=float) + "\n")
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    print(json.dumps(summary, indent=2, default=float))


if __name__ == "__main__":
    main()
