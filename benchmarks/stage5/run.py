"""Stage-5 benchmark (docs/benchmarks/stage5_adjoint.md): adjoint vs finite differences.

Writes results/stage5/{results,summary,manifest}.json.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import spsolve

from widegap import adjoint as A
from widegap.fv import equilibrium as eq
from widegap.fv import reverse_bias as rb
from widegap.fv.drift_diffusion import DDParams, DriftDiffusion1D
from widegap.fv.equilibrium import control_volumes, cv_average_piecewise, neutral_potential
from widegap.manifest import build_manifest

REPO = Path(__file__).resolve().parents[2]
SEEDS = [0, 1, 2, 3, 4]


def direction(x, c, seed):
    xi = (x - x[0]) / (x[-1] - x[0])
    a = np.random.default_rng(seed).standard_normal(4)
    m = sum(a[k] * np.sin((k + 1) * np.pi * xi) for k in range(4))
    m = m / np.max(np.abs(m))
    v = np.maximum(np.abs(c), 1.0) * m
    v[0] = v[-1] = 0.0
    return v


def single_nodes(n_nodes, seed=11, k=5):
    return sorted(
        np.random.default_rng(seed).choice(np.arange(1, n_nodes - 1), size=k, replace=False)
    )


def unit(n, i, c):
    v = np.zeros(n)
    v[i] = max(abs(c[i]), 1.0)
    return v


def check(J_of, grad, v, eps):
    jp, jm = J_of(eps * v), J_of(-eps * v)
    fd = (jp[0] - jm[0]) / (2 * eps)
    ad = float(grad @ v)
    return {
        "adjoint": ad,
        "fd": fd,
        "abs_err": abs(ad - fd),
        "rel_err": abs(ad - fd) / abs(fd) if fd else None,
        "extra": [jp[1], jm[1]],
    }


# ----- A1 -------------------------------------------------------------------
def case_a1():
    d = 1e-8
    x = np.linspace(-40.0, 40.0, 1025)
    c = cv_average_piecewise(x, [0.0], [-1.0, 1.0])
    pl, pr = float(neutral_potential(-1.0, d)), float(neutral_potential(1.0, d))
    psi = A.polish_equilibrium(eq.solve(x, c, d, pl, pr).psi, x, c, d)
    J, g = A.equilibrium_integral_psi(psi, x, d)
    vol = control_volumes(x)

    def J_of(dc):
        cc = c + dc
        s = eq.solve(x, cc, d, pl, pr, psi0=psi)
        p2 = A.polish_equilibrium(s.psi, x, cc, d)
        return float(np.sum(vol * p2)), bool(s.report.converged)

    dirs = [check(J_of, g, direction(x, c, s), 1e-4) for s in SEEDS]
    nodes = [check(J_of, g, unit(len(x), i, c), 1e-4) for i in single_nodes(len(x))]
    return {"J": J, "directional": dirs, "single_node": nodes, "gamma_C": g.tolist()}


# ----- A2 -------------------------------------------------------------------
def _dd_setup(params, v_bias=6, n_el=1024, L=40.0, delta=1e-3):
    x = np.linspace(-L, L, n_el + 1)
    c = cv_average_piecewise(x, [0.0], [-1.0, 1.0])
    pl, pr = float(neutral_potential(-1.0, delta)), float(neutral_potential(1.0, delta))
    eqs = eq.solve(x, c, delta, pl, pr)

    def bc(v):
        return [
            (pl + v, delta * np.exp(pl), delta * np.exp(-pl)),
            (pr, delta * np.exp(pr), delta * np.exp(-pr)),
        ]

    dd = DriftDiffusion1D(x, c, params)
    X = np.empty(3 * len(x))
    X[0::3], X[1::3], X[2::3] = eqs.psi, eqs.n, eqs.p
    for v in range(0, v_bias + 1):  # bias continuation in 1 V_T steps
        X = _pack(dd.solve(X, bc(v), tol_mode="jacobian"))
    return x, c, params, bc(v_bias), A.polish_dd(dd, X, bc(v_bias))


def _pack(sol):
    X = np.empty(3 * len(sol.psi))
    X[0::3], X[1::3], X[2::3] = sol.psi, sol.n, sol.p
    return X


def case_a2(params, timing=False):
    x, c, prm, bcv, X = _dd_setup(params)
    dd = DriftDiffusion1D(x, c, prm)
    J, g = A.dd_mean_current(dd, X)

    def J_of(dc):
        d2 = DriftDiffusion1D(x, c + dc, prm)
        sol = d2.solve(X, bcv, tol_mode="jacobian")
        X2 = A.polish_dd(d2, _pack(sol), bcv)
        p, n_, pp = d2.split(X2)
        jn, jp = d2.fluxes(p, n_, pp)
        return float(np.mean(jn + jp)), bool(sol.report.converged)

    out = {
        "J": J,
        "directional": [check(J_of, g, direction(x, c, s), 1e-3) for s in SEEDS],
        "single_node": [check(J_of, g, unit(len(x), i, c), 1e-3) for i in single_nodes(len(x))],
        "gamma_C": g.tolist(),
    }
    if timing:
        t_adj, t_newton = [], []
        F = dd.residual(X, bcv)
        for _ in range(20):
            t0 = time.perf_counter()
            A.dd_mean_current(dd, X)
            t_adj.append(time.perf_counter() - t0)
            t0 = time.perf_counter()
            spsolve(dd.jacobian(X), -F)
            t_newton.append(time.perf_counter() - t0)
        t0 = time.perf_counter()
        J_of(np.zeros_like(c))
        t_fwd = time.perf_counter() - t0
        n_design = len(x) - 2
        out["timing"] = {
            "adjoint_median_s": float(np.median(t_adj)),
            "newton_iteration_median_s": float(np.median(t_newton)),
            "ratio": float(np.median(t_adj) / np.median(t_newton)),
            "forward_solve_s": t_fwd,
            "fd_gradient_estimate_s": 2 * n_design * t_fwd,
            "fd_over_adjoint": 2 * n_design * t_fwd / float(np.median(t_adj)),
        }
    return out


# ----- A3 / A4 --------------------------------------------------------------
def _b1(h=0.025):
    x = np.linspace(-5.0, 60.0, int(round(65.0 / h)) + 1)
    return x, cv_average_piecewise(x, [0.0], [-100.0, 1.0])


def case_a3_a4(e_crit=15.0, v_r=50):
    x, c = _b1()
    d = rb.Diode1D(x, c, 1e-3)
    sw = rb.sweep(d, range(0, 121))
    psi50 = A.polish_reverse(d, sw.psi[v_r], v_r)
    B, gB, dBdV, k0 = A.breakdown_index(d, psi50, v_r, e_crit)

    def B_of(dc):
        d2 = rb.Diode1D(x, c + dc, 1e-3)
        p, rep = d2.solve(v_r, psi0=psi50)
        p = A.polish_reverse(d2, p, v_r)
        e = d2.element_field(p)
        k = int(np.argmax(np.abs(e)))
        return float(abs(e[k]) / e_crit), k

    a3 = {
        "B": B,
        "dB_dV": dBdV,
        "e_star": k0,
        "directional": [check(B_of, gB, direction(x, c, s), 1e-4) for s in SEEDS],
        "single_node": [check(B_of, gB, unit(len(x), i, c), 1e-4) for i in single_nodes(len(x))],
    }
    # dB/dV check by FD in V_R
    eps_v = 1e-4
    bp_ = B_of_v(d, psi50, v_r + eps_v, e_crit)
    bm_ = B_of_v(d, psi50, v_r - eps_v, e_crit)
    a3["dB_dV_fd"] = (bp_ - bm_) / (2 * eps_v)

    bv, psi_bv, _ = rb.breakdown_voltage(d, sw, e_crit, xtol=1e-12)
    psi_bv = A.polish_reverse(d, psi_bv, bv)
    gBV, kbv = A.breakdown_voltage_sensitivity(d, psi_bv, bv, e_crit)

    def BV_of(dc):
        d2 = rb.Diode1D(x, c + dc, 1e-3)
        sw2 = rb.sweep(d2, range(0, 121))
        bv2, p2, _ = rb.breakdown_voltage(d2, sw2, e_crit, xtol=1e-12)
        k = int(np.argmax(np.abs(d2.element_field(p2))))
        return float(bv2), k

    a4 = {
        "BV": bv,
        "e_star": kbv,
        "directional": [check(BV_of, gBV, direction(x, c, s), 1e-3) for s in SEEDS],
        "single_node": [
            check(BV_of, gBV, unit(len(x), i, c), 1e-3) for i in single_nodes(len(x))[:3]
        ],
    }
    return a3, a4


def B_of_v(d, psi0, v, e_crit):
    p, _ = d.solve(v, psi0=psi0)
    p = A.polish_reverse(d, p, v)
    return d.max_field(p) / e_crit


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=REPO / "results" / "stage5")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "manifest.json").write_text(
        json.dumps(build_manifest(REPO, sys.argv, SEEDS), indent=2) + "\n"
    )

    res = {"A1": case_a1()}
    res["A2"] = case_a2(DDParams(delta=1e-3, recombination=False), timing=True)
    res["A2s"] = case_a2(DDParams(delta=1e-3, tau_n=1e3, tau_p=1e3))
    res["A3"], res["A4"] = case_a3_a4()

    def worst(case):
        # an undefined relative error (FD exactly zero) counts as a failure
        return max(
            np.inf if r["rel_err"] is None else r["rel_err"] for r in res[case]["directional"]
        )

    estar3 = all(r["extra"] == [res["A3"]["e_star"]] * 2 for r in res["A3"]["directional"])
    estar4 = all(r["extra"] == [res["A4"]["e_star"]] * 2 for r in res["A4"]["directional"])
    crit = {
        "S5a": worst("A1") <= 1e-6,
        "S5b": worst("A2") <= 1e-3,
        "S5c": worst("A3") <= 1e-6 and estar3,
        "S5d": worst("A4") <= 1e-4 and estar4,
        "S5e": res["A2"]["timing"]["ratio"] <= 3.0,
    }
    crit = {k: bool(v) for k, v in crit.items()}
    summary = {
        "criteria": crit,
        "gate_passed": all(crit.values()),
        "worst_directional_rel_err": {c: worst(c) for c in ("A1", "A2", "A2s", "A3", "A4")},
        "worst_single_node_rel_err": {
            c: max(
                (r["rel_err"] for r in res[c]["single_node"] if r["rel_err"] is not None),
                default=None,
            )
            for c in ("A1", "A2", "A2s", "A3", "A4")
        },
        "single_node_zero_fd_adjoint_and_abs_err": {
            c: [(r["adjoint"], r["abs_err"]) for r in res[c]["single_node"] if r["rel_err"] is None]
            for c in ("A1", "A2", "A2s", "A3", "A4")
        },
        "timing_A2": res["A2"]["timing"],
        "dB_dV": {"adjoint": res["A3"]["dB_dV"], "fd": res["A3"]["dB_dV_fd"]},
    }
    (args.out / "results.json").write_text(json.dumps(res, indent=2, default=float) + "\n")
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
