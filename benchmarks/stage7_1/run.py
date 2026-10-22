"""Stage-7.1 benchmark (docs/benchmarks/stage7_1_quasi_fermi_dd.md).

Cases Q1-Q4, criteria T1-T8 (gating) and R1-R3 (reported).
Writes <out>/{results,summary,manifest}.json (default results/stage7_1).
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np

from widegap.fv.drift_diffusion import DDParams, DriftDiffusion1D
from widegap.fv.equilibrium import cv_average_piecewise, neutral_potential
from widegap.fv.equilibrium import solve as solve_equilibrium
from widegap.fv.qf_drift_diffusion import QFDriftDiffusion1D, exact_sg_identity_currents, pack
from widegap.manifest import build_manifest
from widegap.materials.registry import Registry
from widegap.problems.drift_diffusion import D1_PARAMS, d1_fields
from widegap.verification import observed_orders

REPO = Path(__file__).resolve().parents[2]
REGISTRY = REPO / "src" / "widegap" / "materials" / "data" / "sic4h.json"
DELTA_Q1 = 1e-3
DELTA_Q3 = 1.0e-24
LOGD_Q3 = float(np.log(DELTA_Q3))


def _newton(sol):
    r = sol.report
    return {
        "converged": bool(r.converged),
        "iterations": r.iterations,
        "limited_steps": r.limited_steps,
        "final_max_residual_over_scale": r.residual_history[-1],
        "last_dX": r.last_dX,
    }


def _jac_check(dd, X, bc, seed=0):
    """Central FD (eps=1e-7), 5 seeded directions, per-row relative to sum_j |J_ij||v_j|."""
    rng = np.random.default_rng(seed)
    J = dd.jacobian(X)
    worst = 0.0
    for _ in range(5):
        v = rng.standard_normal(len(X))
        v[[0, 1, 2, -3, -2, -1]] = 0.0
        e = 1e-7
        fd = (dd.residual(X + e * v, bc) - dd.residual(X - e * v, bc)) / (2 * e)
        jv = J @ v
        den = abs(J) @ np.abs(v)
        m = den > 0
        m[[0, 1, 2, -3, -2, -1]] = False
        worst = max(worst, float(np.max(np.abs(jv - fd)[m] / den[m])))
    return worst


def _diode(delta, a, b, lo, hi, n_el, params):
    x = np.linspace(lo, hi, n_el + 1)
    c = cv_average_piecewise(x, [0.0], [-a, b])
    pl, pr = float(neutral_potential(-a, delta)), float(neutral_potential(b, delta))
    eq = solve_equilibrium(x, c, delta, pl, pr, max_iter=400)
    dd = QFDriftDiffusion1D(x, c, params, float(np.log(delta)))
    return x, c, pl, pr, eq, dd


# --------------------------------------------------------------------- Q1
def run_q1():
    delta = DELTA_Q1
    prm = DDParams(delta=delta, recombination=False)
    x, c, pl, pr, eq, dd = _diode(delta, 1.0, 1.0, -40.0, 40.0, 1024, prm)
    dd3 = DriftDiffusion1D(x, c, prm)
    z = np.zeros_like(x)
    Xq = pack(eq.psi, z, z)
    X3 = pack(eq.psi, eq.n, eq.p)
    nl, pl_ = delta * np.exp(pl), delta * np.exp(-pl)
    nr, pr_ = delta * np.exp(pr), delta * np.exp(-pr)
    out, jac = {}, {}
    for v in range(0, 11):
        bcq = [(pl + v, float(v), float(v)), (pr, 0.0, 0.0)]
        bc3 = [(pl + v, nl, pl_), (pr, nr, pr_)]
        sq = dd.solve(Xq, bcq)
        s3 = dd3.solve(X3, bc3, tol_mode="jacobian")
        if v == 0:
            jac["Q1_initial"] = _jac_check(dd, Xq, bcq)
        Xq, X3 = sq.X, pack(s3.psi, s3.n, s3.p)
        jq, j3 = float(np.mean(sq.j)), float(np.mean(s3.j))
        out[v] = {
            **_newton(sq),
            "stage3b_converged": bool(s3.report.converged),
            "psi_diff_max": float(np.max(np.abs(sq.psi - s3.psi))),
            "j_qf": jq,
            "j_stage3b": j3,
            "j_rel_diff": abs(jq - j3) / abs(j3) if j3 != 0 else None,
        }
        if v == 4:
            jac["Q1_solution"] = _jac_check(dd, Xq, bcq)
    return out, jac


# --------------------------------------------------------------------- Q2
def run_q2():
    levels = [32, 64, 128, 256, 512]
    errs = {"psi": [], "n": [], "p": []}
    newton = {}
    ld = float(np.log(D1_PARAMS["delta"]))
    for N in levels:
        x = np.linspace(0.0, 1.0, N + 1)
        psi, n, p, s_psi, s_n, s_p = d1_fields(x)
        dd = QFDriftDiffusion1D(
            x,
            np.zeros_like(x),
            DDParams(**D1_PARAMS),
            ld,
            {"psi": s_psi, "n": s_n, "p": s_p},
        )

        def to_q(psi_, n_, p_):
            return pack(psi_, ld + psi_ - np.log(n_), np.log(p_) - ld + psi_)

        Xe = to_q(psi, n, p)
        bc = [tuple(Xe[:3]), tuple(Xe[-3:])]
        psi0 = psi[0] + (psi[-1] - psi[0]) * x
        X0 = to_q(psi0, np.ones_like(x), 0.7 * np.ones_like(x))
        sol = dd.solve(X0, bc)
        for k, a, b in (("psi", sol.psi, psi), ("n", sol.n, n), ("p", sol.p, p)):
            errs[k].append(float(np.max(np.abs(a - b))))
        newton[N] = _newton(sol)
    orders = {k: observed_orders(v) for k, v in errs.items()}
    return {"levels": levels, "e_inf": errs, "orders": orders, "newton": newton}


# --------------------------------------------------------------- Q3 / Q4
def run_realistic(recombination, biases, with_jac=False, keep_states=None):
    prm = DDParams(delta=DELTA_Q3, tau_n=1e3, tau_p=1e3, recombination=recombination)
    x, c, pl, pr, eq, dd = _diode(DELTA_Q3, 10.0, 1.0, -30.0, 60.0, 3600, prm)
    z = np.zeros_like(x)
    X = pack(eq.psi, z, z)
    out, jac = {}, {}
    states = {}
    for v in biases:
        bc = [(pl + v, float(v), float(v)), (pr, 0.0, 0.0)]
        if with_jac and v == 0:
            jac["Q3_initial"] = _jac_check(dd, X, bc)
        t0 = time.perf_counter()
        sol = dd.solve(X, bc)
        el = time.perf_counter() - t0
        X = sol.X
        jm = float(np.mean(sol.j))
        rec = {
            **_newton(sol),
            "seconds": el,
            "j_mean": jm,
            "j_rel_spread": float(np.max(np.abs(sol.j - jm)) / abs(jm)) if jm else None,
            "jn_first_elem": float(sol.jn[0]),
            "jp_last_elem": float(sol.jp[-1]),
            "positive": bool(np.all(sol.n > 0) and np.all(sol.p > 0)),
        }
        if v == 0:
            rec["psi_vs_stage2_max"] = float(np.max(np.abs(sol.psi - eq.psi)))
            rec["phi_max_abs"] = float(max(np.max(np.abs(sol.phi_n)), np.max(np.abs(sol.phi_p))))
        if not recombination and v > 0:
            jn_e, jp_e = exact_sg_identity_currents(
                x, sol.psi, v, 0.0, v, 0.0, LOGD_Q3, prm.mu_n, prm.mu_p
            )
            rec["jn_identity"], rec["jp_identity"] = jn_e, jp_e
            rec["jn_rel_err"] = abs(sol.jn[0] - jn_e) / abs(jn_e)
            rec["jp_rel_err"] = abs(sol.jp[-1] - jp_e) / abs(jp_e)
        if with_jac and v == 60:
            jac["Q3_solution"] = _jac_check(dd, X, bc)
        if keep_states and v in keep_states:
            states[v] = sol
        out[v] = rec
    return out, jac, (x, c, pl, pr, eq)


def ideality(res, v0, v1):
    return {
        f"{v}-{v + 1}": float(np.log(res[v + 1]["j_mean"] / res[v]["j_mean"]))
        for v in range(v0, v1)
        if v in res and v + 1 in res and res[v]["j_mean"] > 0 and res[v + 1]["j_mean"] > 0
    }


# --------------------------------------------------------------------- R1
def run_r1(q3, geom, biases=(0, 5, 20, 60, 100)):
    x, c, pl, pr, eq = geom
    delta = DELTA_Q3
    prm = DDParams(delta=delta, recombination=False)
    dd3 = DriftDiffusion1D(x, c, prm)
    nl, pl_ = delta * np.exp(pl), delta * np.exp(-pl)
    nr, pr_ = delta * np.exp(pr), delta * np.exp(-pr)
    X_eq = pack(eq.psi, eq.n, eq.p)
    out = {}
    for mode in ("terms", "jacobian"):
        for v in biases:
            key = f"{mode}_V{v}"
            bc = [(pl + v, nl, pl_), (pr, nr, pr_)]
            try:
                with warnings.catch_warnings(), np.errstate(all="ignore"):
                    warnings.simplefilter("ignore")
                    s = dd3.solve(X_eq, bc, tol_mode=mode)
                    jm = float(np.mean(s.j))
                    q = q3[v]
                    out[key] = {
                        "converged": bool(s.report.converged),
                        "iterations": s.report.iterations,
                        "limited_steps": s.report.limited_steps,
                        "final_residuals": list(s.report.residual_history[-1]),
                        "tolerance": list(s.report.tolerance),
                        "j_mean": jm,
                        "jn_first_elem": float(s.jn[0]),
                        "jp_last_elem": float(s.jp[-1]),
                        "j_rel_spread": float(np.max(np.abs(s.j - jm)) / abs(jm)) if jm else None,
                        "jn_first_vs_identity_rel": float(
                            abs(s.jn[0] - q["jn_identity"]) / abs(q["jn_identity"])
                        )
                        if v > 0
                        else None,
                        "jp_last_vs_identity_rel": float(
                            abs(s.jp[-1] - q["jp_identity"]) / abs(q["jp_identity"])
                        )
                        if v > 0
                        else None,
                        "j_vs_qf_rel": float(abs(jm - q["j_mean"]) / abs(q["j_mean"])),
                        "min_n": float(np.min(s.n)),
                        "min_p": float(np.min(s.p)),
                    }
            except Exception as exc:  # documented outcome, must not crash
                out[key] = {"error": f"{type(exc).__name__}: {exc}"}
    return out


# --------------------------------------------------------------------- R3
def r3_units(q3):
    reg = Registry.from_json(REGISTRY)
    mx, mn = reg.value("mu_n_max_par_c"), reg.value("mu_n_min_par_c")
    nref, al = reg.value("N_n_ref_par_c"), reg.value("alpha_n_par_c")
    eps_r = reg.value("eps_r_parallel_c")
    n0, nd = 1e16, 1e16  # cm^-3; N_D = b N0 with b = 1
    mu0 = mn + (mx - mn) / (1.0 + (nd / nref) ** al)  # cm^2/(V s)
    q, kb, t, eps0 = 1.602176634e-19, 1.380649e-23, 300.0, 8.8541878128e-14  # eps0 in F/cm
    vt = kb * t / q
    ld = float(np.sqrt(eps_r * eps0 * vt / (q * n0)))  # cm
    j_unit = q * mu0 * vt * n0 / ld  # A/cm^2
    sel = [v for v in (20, 40, 60, 80, 90, 100, 105, 110) if v in q3]
    return {
        "label": "ILLUSTRATION ONLY: mu_n=mu_p=mu0 (Caughey-Thomas at N_D), complete ionization",
        "mu0_cm2_per_Vs": mu0,
        "eps_r": eps_r,
        "V_T_volt": vt,
        "L_D_cm": ld,
        "j_unit_A_per_cm2": j_unit,
        "j_A_per_cm2": {str(v): q3[v]["j_mean"] * j_unit for v in sel},
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=REPO / "results" / "stage7_1")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "manifest.json").write_text(
        json.dumps(build_manifest(REPO, sys.argv, [0]), indent=2) + "\n"
    )

    q1, jac1 = run_q1()
    q2 = run_q2()
    biases = list(range(0, 111))
    q3, jac3, geom = run_realistic(False, biases, with_jac=True)
    q4, _, _ = run_realistic(True, biases)
    r1 = run_r1(q3, geom)
    jac = {**jac1, **jac3}
    id3 = ideality(q3, 95, 102)
    id3_all = ideality(q3, 0, 110)
    id4_all = ideality(q4, 0, 110)
    id4 = ideality(q4, 95, 102)

    t7_v = [v for v in biases if abs(q3[v]["j_mean"]) >= 1e-6]
    crit = {
        "T1": all(q1[v]["stage3b_converged"] and q1[v]["psi_diff_max"] <= 1e-9 for v in q1)
        and all(q1[v]["j_rel_diff"] <= 1e-6 for v in q1 if v >= 3),
        "T2": all(all(1.8 <= o <= 2.2 for o in ords[-3:]) for ords in q2["orders"].values()),
        "T3": all(val <= 1e-6 for val in jac.values()),
        "T4": all(q3[v]["converged"] and q3[v]["iterations"] <= 60 for v in biases),
        "T5": q3[0]["psi_vs_stage2_max"] <= 1e-9 and q3[0]["phi_max_abs"] <= 1e-12,
        "T6": all(q3[v]["jn_rel_err"] <= 1e-8 and q3[v]["jp_rel_err"] <= 1e-8 for v in biases[1:]),
        "T7": all(q3[v]["j_rel_spread"] <= 1e-4 for v in t7_v),
        "T8": len(id3) == 7 and all(0.97 <= q <= 1.03 for q in id3.values()),
    }
    crit = {k: bool(v) for k, v in crit.items()}

    t1_vals = {
        "max_psi_diff": max(q1[v]["psi_diff_max"] for v in q1),
        "max_j_rel_diff_V_ge_3": max(q1[v]["j_rel_diff"] for v in q1 if v >= 3),
    }
    t6_vals = {
        "max_jn_rel_err": max(q3[v]["jn_rel_err"] for v in biases[1:]),
        "max_jp_rel_err": max(q3[v]["jp_rel_err"] for v in biases[1:]),
        "worst_V_jn": max(biases[1:], key=lambda v: q3[v]["jn_rel_err"]),
        "worst_V_jp": max(biases[1:], key=lambda v: q3[v]["jp_rel_err"]),
    }
    t7_vals = {
        "n_biases_gated": len(t7_v),
        "first_gated_V": t7_v[0] if t7_v else None,
        "max_rel_spread": max((q3[v]["j_rel_spread"] for v in t7_v), default=None),
    }
    reported = {
        "T1_values": t1_vals,
        "T2_orders": q2["orders"],
        "T3_jacobian_rel_err": jac,
        "T4_max_iterations": max(q3[v]["iterations"] for v in biases),
        "T4_all_converged": all(q3[v]["converged"] for v in biases),
        "T4_failed_V": [v for v in biases if not q3[v]["converged"]],
        "T5": {"psi_vs_stage2": q3[0]["psi_vs_stage2_max"], "phi_max": q3[0]["phi_max_abs"]},
        "T6_values": t6_vals,
        "T7_values": t7_vals,
        "T8_ideality": id3,
        "R1": r1,
        "R2": {
            "all_converged": all(q4[v]["converged"] for v in biases),
            "failed_V": [v for v in biases if not q4[v]["converged"]],
            "max_iterations": max(q4[v]["iterations"] for v in biases),
            "ideality_95_102": id4,
        },
        "R3": r3_units(q3),
    }
    res = {
        "Q1": {str(k): v for k, v in q1.items()},
        "Q2": q2,
        "Q3": {str(k): v for k, v in q3.items()},
        "Q3_ideality_all": id3_all,
        "Q4": {str(k): v for k, v in q4.items()},
        "Q4_ideality_all": id4_all,
        "R1": r1,
        "jacobian_rel_err": jac,
        "R3": reported["R3"],
    }
    summary = {"criteria": crit, "gate_passed": all(crit.values()), "reported": reported}
    (args.out / "results.json").write_text(json.dumps(res, indent=2, default=float) + "\n")
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    print(json.dumps(summary, indent=2, default=float))


if __name__ == "__main__":
    main()
