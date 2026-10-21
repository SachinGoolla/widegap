"""Stage-4 benchmark (docs/benchmarks/stage4_breakdown_index.md and stage4b).

--criteria v1 : the original pre-registered criteria G1-G9 (meshes 0.05/0.025/0.0125).
--criteria 4b : re-registered criteria (docs/benchmarks/stage4b_first_integral.md).
Writes <out>/{results,summary,manifest}.json.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from widegap.fv.equilibrium import cv_average_piecewise, neutral_potential
from widegap.fv.reverse_bias import Diode1D, breakdown_voltage, jacobian_banded, residual, sweep
from widegap.manifest import build_manifest

REPO = Path(__file__).resolve().parents[2]
DELTA, ECRIT, VMAX = 1e-3, 15.0, 120
CASES = {
    "B1": dict(lo=-5.0, hi=60.0, breaks=[0.0], vals=[-100.0, 1.0]),
    "B2": dict(lo=-5.0, hi=13.0, breaks=[0.0, 8.0], vals=[-100.0, 1.0, 100.0]),
}


def make(case, h):
    c = CASES[case]
    x = np.linspace(c["lo"], c["hi"], int(round((c["hi"] - c["lo"]) / h)) + 1)
    return Diode1D(x, cv_average_piecewise(x, c["breaks"], c["vals"]), DELTA)


def jac_check(d: Diode1D, psi, v_r, seed=0):
    rng = np.random.default_rng(seed)
    args = (d.x, d.c_cv, d.delta, 0.0, -float(v_r))
    ab = jacobian_banded(psi, d.x, d.delta, 0.0, -float(v_r))
    worst = 0.0
    for _ in range(5):
        v = rng.standard_normal(len(psi))
        v[0] = v[-1] = 0.0
        jv = ab[1] * v
        jv[:-1] += ab[0, 1:] * v[1:]
        jv[1:] += ab[2, :-1] * v[:-1]
        e = 1e-7
        fd = (residual(psi + e * v, *args) - residual(psi - e * v, *args)) / (2 * e)
        worst = max(worst, float(np.max(np.abs(jv[1:-1] - fd[1:-1])) / np.max(np.abs(jv[1:-1]))))
    return worst


def first_integral_field(case, d: Diode1D, psi, v_r, j0):
    """Exact |E| at node j0 from the cathode-side first integral of -psi'' = rho(psi).

    E(x0)^2 / 2 = int_{psi(x0)}^{psi_cathode} rho(psi) dpsi, with
    rho = delta(e^{phi_p - psi} - e^{psi}) + C. On each element C is the true
    (piecewise-constant) doping, so the integral is exact given nodal psi; the
    Boltzmann terms telescope exactly. Reference accuracy is O(h^2) via psi.
    """
    phi_p = -float(v_r)
    cfg = CASES[case]
    xm = 0.5 * (d.x[1:] + d.x[:-1])
    c_true = np.asarray(cfg["vals"])[np.searchsorted(cfg["breaks"], xm)]
    a, b = psi[j0:-1], psi[j0 + 1 :]
    total = np.sum(c_true[j0:] * (b - a))
    total += DELTA * (np.exp(phi_p - psi[j0]) - np.exp(phi_p - psi[-1]))
    total -= DELTA * (np.exp(psi[-1]) - np.exp(psi[j0]))
    return float(np.sqrt(2.0 * max(float(total), 0.0)))


def depletion_bv(case):
    if case == "B1":
        vbi = float(neutral_potential(1.0, DELTA) - neutral_potential(-100.0, DELTA))
        return ECRIT**2 * 101.0 / 200.0 - vbi + 2.0
    vbi = float(neutral_potential(100.0, DELTA) - neutral_potential(-100.0, DELTA))
    w = 8.0
    return (ECRIT - w / 2.0) * w - vbi


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--criteria", choices=["v1", "4b"], default="v1")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()
    out = args.out or REPO / "results" / ("stage4" if args.criteria == "v1" else "stage4b")
    out.mkdir(parents=True, exist_ok=True)
    (out / "manifest.json").write_text(
        json.dumps(build_manifest(REPO, sys.argv, [0]), indent=2) + "\n"
    )

    meshes = (
        [0.05, 0.025, 0.0125] if args.criteria == "v1" else [0.05, 0.025, 0.0125, 0.00625, 0.003125]
    )
    vs = list(range(0, VMAX + 1))
    res, crit = {}, {}
    all_conv, g1 = True, True
    for case in CASES:
        r = {"BV": {}, "B_minus_1": {}, "max_newton_iterations": {}}
        for h in meshes:
            d = make(case, h)
            sw = sweep(d, vs)
            bv, psi_bv, reps = breakdown_voltage(d, sw, ECRIT)
            allreps = sw.reports + reps
            all_conv &= all(q.converged and q.iterations <= 60 for q in allreps)
            r["BV"][h] = bv
            r["B_minus_1"][h] = d.max_field(psi_bv) / ECRIT - 1.0
            r["max_newton_iterations"][h] = max(q.iterations for q in allreps)
            g1 &= abs(r["B_minus_1"][h]) <= 1e-9
            if h == 0.025:
                V, E = np.array(sw.v), np.array(sw.e_max)
                r["e_max_h0.025"] = E.tolist()
                r["jacobian"] = {
                    "V50": jac_check(d, sw.psi[50], 50),
                    "V0_initial": jac_check(d, neutral_potential(d.c_cv, DELTA), 0),
                }
                if case == "B1":
                    m = (V >= 20) & (V <= 80)
                    r["slope_E2"] = float(np.polyfit(V[m], E[m] ** 2, 1)[0])
                    r["G9_cathode_field_ratio"] = float(
                        abs(d.element_field(psi_bv)[-1]) / d.max_field(psi_bv)
                    )
                    psi = sw.psi[50]
                    e = np.abs(d.element_field(psi))
                    xm = 0.5 * (d.x[1:] + d.x[:-1])
                    emax = e.max()
                    xn = float(xm[(xm > 0) & (e < 0.01 * emax)][0])
                    w = (xm >= 2.0) & (xm <= xn - 3.0)
                    r["G6_slope_v1_window"] = float(np.polyfit(xm[w], e[w], 1)[0])
                    n = DELTA * np.exp(psi)
                    nm = 0.5 * (n[1:] + n[:-1])
                    w2 = (xm >= 1.0) & (nm < 0.01)
                    r["G6_slope_4b_window"] = float(np.polyfit(xm[w2], e[w2], 1)[0])
                else:
                    m = (V >= 40) & (V <= 100)
                    r["slope_E"] = float(np.polyfit(V[m], E[m], 1)[0])
            if args.criteria == "4b":
                j0 = int(np.argmin(np.abs(d.x)))
                fi = {}
                for v in (20, 50, 80):
                    psi = sw.psi[v]
                    fi[v] = {
                        "numeric_e_junction_element": float(abs(d.element_field(psi)[j0])),
                        "first_integral": first_integral_field(case, d, psi, v, j0),
                    }
                r.setdefault("first_integral", {})[h] = fi
        bvs = [r["BV"][h] for h in meshes]
        r["BV_orders"] = [
            float(np.log(abs(bvs[i] - bvs[i + 1]) / abs(bvs[i + 1] - bvs[i + 2])) / np.log(2))
            for i in range(len(bvs) - 2)
        ]
        r["BV_depletion"] = depletion_bv(case)
        r["BV_over_depletion_h0.025"] = r["BV"][0.025] / r["BV_depletion"]
        res[case] = r

    b1, b2 = res["B1"], res["B2"]
    jac_ok = all(v <= 1e-6 for c in (b1, b2) for v in c["jacobian"].values())
    if args.criteria == "v1":
        crit = {
            "G1": g1,
            "G2": all(
                abs(c["BV"][0.025] - c["BV"][0.0125]) / c["BV"][0.0125] <= 0.01 for c in (b1, b2)
            ),
            "G3": all(abs(c["BV"][0.025] / c["BV_depletion"] - 1) <= 0.10 for c in (b1, b2)),
            "G4": abs(b1["slope_E2"] / (2 * 100 / 101) - 1) <= 0.03,
            "G5": abs(b2["slope_E"] / (1 / 8.0) - 1) <= 0.10,
            "G6": abs(b1["G6_slope_v1_window"] + 1.0) <= 0.01,
            "G7": jac_ok,
            "G8": all_conv,
            "G9": b1["G9_cathode_field_ratio"] <= 1e-3,
        }
    else:
        for c in (b1, b2):
            bvs = [c["BV"][h] for h in meshes]
            extrap = 2 * bvs[-1] - bvs[-2]  # first-order Richardson
            c["BV_richardson"] = extrap
            c["BV_finest_vs_extrap"] = abs(bvs[-1] - extrap) / extrap
            fi_err = []
            for v in (20, 50, 80):
                e1 = c["first_integral"][meshes[-2]][v]["numeric_e_junction_element"]
                e2 = c["first_integral"][meshes[-1]][v]["numeric_e_junction_element"]
                e_ex = 2 * e2 - e1
                ref = c["first_integral"][meshes[-1]][v]["first_integral"]
                fi_err.append(abs(e_ex - ref) / ref)
            c["first_integral_rel_err_extrapolated"] = fi_err
        crit = {
            "G1": g1,
            "G2'": all(
                all(0.8 <= q <= 1.2 for q in c["BV_orders"][-2:])
                and c["BV_finest_vs_extrap"] <= 0.01
                for c in (b1, b2)
            ),
            "G3'": all(max(c["first_integral_rel_err_extrapolated"]) <= 2e-3 for c in (b1, b2)),
            "G6'": abs(b1["G6_slope_4b_window"] + 1.0) <= 0.01,
            "G7": jac_ok,
            "G8": all_conv,
            "G9": b1["G9_cathode_field_ratio"] <= 1e-3,
        }
    crit = {k: bool(v) for k, v in crit.items()}
    summary = {"criteria_set": args.criteria, "criteria": crit, "gate_passed": all(crit.values())}
    (out / "results.json").write_text(json.dumps(res, indent=2, default=float) + "\n")
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
