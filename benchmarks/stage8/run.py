"""Stage-8 benchmark (docs/benchmarks/stage8_2d_fv.md): P1-P9 (gating), R1-R2 (reported).

Writes <out>/{results,summary,manifest}.json (default results/stage8).
"""

from __future__ import annotations

import argparse
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

from widegap.fv.drift_diffusion import DDParams
from widegap.fv.equilibrium import cv_average_piecewise, neutral_potential
from widegap.fv.equilibrium import solve as solve_eq1d
from widegap.fv.qf_drift_diffusion import QFDriftDiffusion1D
from widegap.fv.qf_drift_diffusion import pack as pack1d
from widegap.fv2d.equilibrium import solve_equilibrium
from widegap.fv2d.mesh import TensorMesh
from widegap.fv2d.poisson import boundary_fluxes, solve_poisson
from widegap.fv2d.qf_dd import QFDriftDiffusion, pack
from widegap.manifest import build_manifest
from widegap.verification import observed_orders

REPO = Path(__file__).resolve().parents[2]
LEVELS = [16, 32, 64, 128, 256]
DELTA_Q3 = 1.0e-24
LOGD_Q3 = float(np.log(DELTA_Q3))


def rss_mib():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def last3_ok(orders, lo, hi):
    return bool(all(lo <= o <= hi for o in orders[-3:]))


# ------------------------------------------------------------ P1 / P3 / P4
def p1_fields(x, y):
    u = np.sin(np.pi * x) * np.sin(np.pi * y) + x * y
    ux = np.pi * np.cos(np.pi * x) * np.sin(np.pi * y) + y
    uy = np.pi * np.sin(np.pi * x) * np.cos(np.pi * y) + x
    lap = -2 * np.pi**2 * np.sin(np.pi * x) * np.sin(np.pi * y)
    eps = 1 + x**2 + y**2
    f = -(eps * lap + 2 * x * ux + 2 * y * uy)
    return u, f, eps


def eps_fun(x, y):
    return 1 + x**2 + y**2


def graded(n0, ratio, level):
    """Geometric spacing (ratio) with n0 intervals on [0,1], each interval halved `level` times."""
    h = ratio ** np.arange(n0)
    x = np.r_[0.0, np.cumsum(h)] / h.sum()
    for _ in range(level):
        x = np.sort(np.r_[x, 0.5 * (x[:-1] + x[1:])])
    x[-1] = 1.0
    return x


def solve_p1(mesh):
    xy = mesh.node_xy
    u, f, _ = p1_fields(xy[:, 0], xy[:, 1])
    bn = mesh.boundary()
    uh = solve_poisson(mesh, mesh.cell_coeff(eps_fun), f, bn, u[bn])
    return uh, u


def run_p1():
    errs = []
    for N in LEVELS:
        g = np.linspace(0, 1, N + 1)
        uh, u = solve_p1(TensorMesh(g, g))
        errs.append(float(np.max(np.abs(uh - u))))
    return {"N": LEVELS, "e_inf": errs, "orders": observed_orders(errs)}


def run_p3():
    errs = []
    for lev in range(len(LEVELS)):
        mesh = TensorMesh(graded(16, 1.05, lev), graded(16, 0.95, lev))
        uh, u = solve_p1(mesh)
        errs.append(float(np.max(np.abs(uh - u))))
    return {"N": LEVELS, "e_inf": errs, "orders": observed_orders(errs)}


def run_p4(N=128):
    g = np.linspace(0, 1, N + 1)
    mesh = TensorMesh(g, g)
    xy = mesh.node_xy
    u, f, _ = p1_fields(xy[:, 0], xy[:, 1])
    dn = np.unique(np.r_[mesh.left, mesh.right, mesh.bottom])
    eps = mesh.cell_coeff(eps_fun)
    uh = solve_poisson(mesh, eps, f, dn, u[dn])
    fl, bal, scale = boundary_fluxes(mesh, eps, uh, f, dn)
    return {
        "N": N,
        "sum_boundary_flux": float(fl.sum()),
        "sum_Vf": float(np.sum(mesh.vol * f)),
        "balance": bal,
        "scale_sum_V_abs_f": scale,
        "relative_balance": abs(bal) / scale,
    }


# --------------------------------------------------------------------- P2
def run_p2():
    e1, kap = 1.0, 2.5
    s1 = 1.0 / (0.5 + 0.5 * e1 / kap)  # slope on the eps=1 side; s2 = s1 e1/kap
    s2 = s1 * e1 / kap

    def exact(t):
        return np.where(t < 0.5, s1 * t, 0.5 * s1 + s2 * (t - 0.5))

    out = {}
    for name, (nx, ny) in {
        "x_interface_32x16cells": (32, 16),
        "y_interface_16x32cells": (16, 32),
    }.items():
        mesh = TensorMesh(np.linspace(0, 1, nx + 1), np.linspace(0, 1, ny + 1))
        ax = 0 if name.startswith("x") else 1
        eps = np.where(mesh.cell_centres[:, ax] < 0.5, e1, kap)
        lo = mesh.left if ax == 0 else mesh.bottom
        hi = mesh.right if ax == 0 else mesh.top
        dn = np.r_[lo, hi]
        dv = np.r_[np.zeros(len(lo)), np.ones(len(hi))]
        u = solve_poisson(mesh, eps, np.zeros(mesh.n_nodes), dn, dv)
        out[name] = {"max_abs_err": float(np.max(np.abs(u - exact(mesh.node_xy[:, ax]))))}
    return out


# --------------------------------------------------------------------- P5
def p5_fields(x, y, delta):
    T = np.tanh(4 * (x - 0.5))
    c = np.cos(np.pi * y / 2)
    psi = 3 * T * c
    lap = -96 * T * (1 - T**2) * c - (np.pi / 2) ** 2 * psi
    s = -lap - delta * (np.exp(-psi) - np.exp(psi))
    return psi, s


def run_p5(delta=1e-3):
    errs, newton = [], {}
    for N in LEVELS:
        g = np.linspace(0, 1, N + 1)
        mesh = TensorMesh(g, g)
        psi, s = p5_fields(mesh.node_xy[:, 0], mesh.node_xy[:, 1], delta)
        bn = mesh.boundary()
        sol = solve_equilibrium(mesh, 1.0, np.zeros(mesh.n_nodes), delta, bn, psi[bn], s_nodes=s)
        errs.append(float(np.max(np.abs(sol.psi - psi))))
        r = sol.report
        newton[N] = {
            "converged": r.converged,
            "iterations": r.iterations,
            "limited": r.limited_steps,
        }
    return {"N": LEVELS, "e_inf": errs, "orders": observed_orders(errs), "newton": newton}


# --------------------------------------------------------------------- P6
def run_p6():
    delta, a, b = 1e-8, 1.0, 1.0
    x = np.linspace(-40.0, 40.0, 1025)
    c = cv_average_piecewise(x, [0.0], [-a, b])
    pl, pr = float(neutral_potential(-a, delta)), float(neutral_potential(b, delta))
    ref = solve_eq1d(x, c, delta, pl, pr)
    ny = 8
    mesh = TensorMesh(x, np.linspace(0.0, 1.0, ny + 1))
    c2 = np.repeat(c, ny + 1)  # node k = i*(ny+1)+j
    dn = np.r_[mesh.left, mesh.right]
    dv = np.r_[np.full(ny + 1, pl), np.full(ny + 1, pr)]
    sol = solve_equilibrium(mesh, 1.0, c2, delta, dn, dv)
    P = sol.psi.reshape(len(x), ny + 1)
    diff = np.abs(P - ref.psi[:, None])
    return {
        "max_row_diff": float(diff.max()),
        "per_row_max": [float(d) for d in diff.max(axis=0)],
        "newton_converged": sol.report.converged,
        "iterations": sol.report.iterations,
        "stage2_converged": bool(ref.report.converged),
    }


# --------------------------------------------------------- Q3 set-up / P7 / P9
def q3_setup(lo=-30.0, hi=60.0, n_el=3600, a=10.0, b=1.0):
    x = np.linspace(lo, hi, n_el + 1)
    c = cv_average_piecewise(x, [0.0], [-a, b])
    pl, pr = float(neutral_potential(-a, DELTA_Q3)), float(neutral_potential(b, DELTA_Q3))
    eq = solve_eq1d(x, c, DELTA_Q3, pl, pr, max_iter=400)
    return x, c, pl, pr, eq


def run_p7(
    ny=4, ramp=(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100, 105), targets=(0, 50, 100, 105)
):
    x, c, pl, pr, eq = q3_setup()
    prm = DDParams(delta=DELTA_Q3, recombination=False)
    dd1 = QFDriftDiffusion1D(x, c, prm, LOGD_Q3)
    mesh = TensorMesh(x, np.linspace(0.0, 1.0, ny + 1))
    dd2 = QFDriftDiffusion(mesh, np.repeat(c, ny + 1), prm, LOGD_Q3)
    Ly = 1.0
    z = np.zeros_like(x)
    X1 = pack1d(eq.psi, z, z)
    X2 = pack(np.repeat(eq.psi, ny + 1), np.zeros(mesh.n_nodes), np.zeros(mesh.n_nodes))
    dn = np.r_[mesh.left, mesh.right]
    nl = len(mesh.left)
    out = {}
    for v in ramp:
        bc = [(pl + v, float(v), float(v)), (pr, 0.0, 0.0)]
        s1 = dd1.solve(X1, bc)
        dv = np.vstack([np.tile(bc[0], (nl, 1)), np.tile(bc[1], (nl, 1))])
        t0 = time.perf_counter()
        s2 = dd2.solve(X2, dn, dv)
        el = time.perf_counter() - t0
        X1, X2 = s1.X, s2.X
        if v not in targets:
            continue
        P = s2.psi.reshape(len(x), ny + 1)
        Pn = s2.phi_n.reshape(len(x), ny + 1)
        Pp = s2.phi_p.reshape(len(x), ny + 1)
        d_psi = np.abs(P - s1.psi[:, None]).max()
        d_pn = np.abs(Pn - s1.phi_n[:, None]).max()
        d_pp = np.abs(Pp - s1.phi_p[:, None]).max()
        j1 = float(np.mean(s1.j))
        cuts = np.array([dd2.cut_current(s2, 0, k) for k in range(len(x) - 1)]) / Ly
        j2 = float(np.mean(cuts))
        jc = dd2.contact_current(s2, mesh.left) / Ly
        out[v] = {
            "converged_1d": bool(s1.report.converged),
            "converged_2d": bool(s2.report.converged),
            "iterations_2d": s2.report.iterations,
            "seconds_2d": el,
            "max_abs_diff_psi": float(d_psi),
            "max_abs_diff_phi_n": float(d_pn),
            "max_abs_diff_phi_p": float(d_pp),
            "j_1d": j1,
            "j_2d_mean_cut": j2,
            "j_rel_diff_mean_cut": abs(j2 - j1) / abs(j1) if j1 else None,
            "j_2d_left_contact": jc,
            "j_1d_first_edge": float(s1.j[0]),
            "j_rel_diff_contact": abs(jc - s1.j[0]) / abs(s1.j[0]) if s1.j[0] else None,
            "gated_current": bool(abs(j1) >= 1e-6),
        }
    return out


def run_p9():
    path = REPO / "results" / "stage7_1" / "results.json"
    x, c, pl, pr, eq = q3_setup()
    prm = DDParams(delta=DELTA_Q3, tau_n=1e3, tau_p=1e3, recombination=False)
    mesh = TensorMesh(x)
    dd = QFDriftDiffusion(mesh, c, prm, LOGD_Q3)
    X = pack(eq.psi, np.zeros_like(x), np.zeros_like(x))
    if path.exists():
        ref = {int(k): v["j_mean"] for k, v in json.loads(path.read_text())["Q3"].items()}
        src = "results/stage7_1/results.json (committed Q3 j_mean)"
        ref_solver = None
    else:
        ref, ref_solver = {}, QFDriftDiffusion1D(x, c, prm, LOGD_Q3)
        src = "QFDriftDiffusion1D run in-script"
    ref1 = ref_solver
    X1 = pack1d(eq.psi, np.zeros_like(x), np.zeros_like(x))
    rel, conv = {}, True
    for v in range(0, 111):
        bc = [(pl + v, float(v), float(v)), (pr, 0.0, 0.0)]
        s = dd.solve(X, [0, len(x) - 1], np.array(bc))
        X = s.X
        conv &= s.report.converged
        j = float(np.mean(s.j))
        if ref1 is not None:
            s1 = ref1.solve(X1, bc)
            X1 = s1.X
            ref[v] = float(np.mean(s1.j))
        rel[v] = abs(j - ref[v]) / abs(ref[v]) if ref[v] else float("inf")
    big = [v for v in rel if abs(ref[v]) >= 1e-6]
    # resolved currents: first-element electron and last-element hole current are
    # well determined at every bias (the mean of jn+jp is round-off noise when the
    # true current is below ~1e-28, see stage-7.1 j_rel_spread)
    rel_res = {}
    if path.exists():
        q3 = json.loads(path.read_text())["Q3"]
        X = pack(eq.psi, np.zeros_like(x), np.zeros_like(x))
        for v in range(0, 111):
            bc = [(pl + v, float(v), float(v)), (pr, 0.0, 0.0)]
            s = dd.solve(X, [0, len(x) - 1], np.array(bc))
            X = s.X
            if v > 0:
                rel_res[v] = max(
                    abs(s.jn[0] - q3[str(v)]["jn_first_elem"]) / abs(q3[str(v)]["jn_first_elem"]),
                    abs(s.jp[-1] - q3[str(v)]["jp_last_elem"]) / abs(q3[str(v)]["jp_last_elem"]),
                )
        spread_ok = [
            v
            for v in rel
            if q3[str(v)]["j_rel_spread"] is not None and q3[str(v)]["j_rel_spread"] <= 1e-3
        ]
    return {
        "reference": src,
        "all_converged": bool(conv),
        "max_rel_diff_all": max(rel.values()),
        "worst_V_all": max(rel, key=rel.get),
        "max_rel_diff_abs_j_ge_1e-6": max(rel[v] for v in big),
        "n_biases": len(rel),
        "max_rel_diff_resolved_jn_first_jp_last_V1_110": max(rel_res.values()) if rel_res else None,
        "worst_V_resolved": max(rel_res, key=rel_res.get) if rel_res else None,
        "max_rel_diff_j_mean_where_committed_spread_le_1e-3": max(rel[v] for v in spread_ok)
        if rel_res and spread_ok
        else None,
        "first_V_with_spread_le_1e-3": min(spread_ok) if rel_res and spread_ok else None,
        "n_V_noise_dominated_j_mean_rel_gt_1e-12": sum(1 for v in rel if rel[v] > 1e-12),
        "rel_diff_at": {str(v): rel[v] for v in (0, 1, 50, 100, 105, 110)},
    }


# --------------------------------------------------------------------- P8
def jac_check(dd, X, dn, dv, rng, n_dirs=5, e=1e-7):
    J = dd.jacobian(X, dn)
    worst = 0.0
    bcrows = np.concatenate([3 * dn, 3 * dn + 1, 3 * dn + 2])
    for _ in range(n_dirs):
        v = rng.standard_normal(len(X))
        v[bcrows] = 0.0
        fd = (dd.residual(X + e * v, dn, dv) - dd.residual(X - e * v, dn, dv)) / (2 * e)
        den = abs(J) @ np.abs(v)
        m = den > 0
        m[bcrows] = False
        worst = max(worst, float(np.max(np.abs(J @ v - fd)[m] / den[m])))
    return worst


def run_p8():
    out = {}
    rng = np.random.default_rng(8)
    for delta, rec in ((1e-3, True), (1e-24, True), (1e-24, False)):
        for seed in range(3):
            ld = float(np.log(delta))
            gx = np.cumsum(np.r_[0.0, rng.uniform(0.5, 1.5, 8)])
            gy = np.cumsum(np.r_[0.0, rng.uniform(0.5, 1.5, 6)])
            mesh = TensorMesh(gx, gy)
            prm = DDParams(delta=delta, tau_n=1e3, tau_p=1e3, recombination=rec)
            c = rng.standard_normal(mesh.n_nodes) * 2
            dd = QFDriftDiffusion(mesh, c, prm, ld)
            xy = mesh.node_xy
            psi = 4 * np.sin(0.7 * xy[:, 0] + 0.3 * seed) * np.cos(
                0.5 * xy[:, 1]
            ) + 0.3 * rng.standard_normal(mesh.n_nodes)
            X = pack(
                psi, 3 * rng.standard_normal(mesh.n_nodes), 3 * rng.standard_normal(mesh.n_nodes)
            )
            dn = np.r_[mesh.left, mesh.right]
            dv = rng.standard_normal((len(dn), 3))
            st = rng.bit_generator.state
            r7 = jac_check(dd, X, dn, dv, rng, e=1e-7)
            rng.bit_generator.state = st
            r5 = jac_check(
                dd, X, dn, dv, rng, e=1e-5
            )  # same directions, larger FD step (informative)
            out[f"delta{delta:g}_srh{int(rec)}_seed{seed}"] = {"fd_1e-7": r7, "fd_1e-5": r5}
    return out


# --------------------------------------------------------------------- R1
def diode_2d(nx, ny, lx=10.0, ly=10.0, delta=1e-6, a=1.0, b=1.0, anode_y=(0.0, 10.0), srh=True):
    x = np.linspace(-lx, lx, nx + 1)
    y = np.linspace(0.0, ly, ny + 1)
    mesh = TensorMesh(x, y)
    cc = np.where(mesh.cell_centres[:, 0] < 0.0, -a, b)
    # junction on a grid line when nx is even: cell-wise doping, exact CV average
    c = mesh.cv_average(cc)
    xy = mesh.node_xy
    anode = mesh.left[(xy[mesh.left, 1] >= anode_y[0]) & (xy[mesh.left, 1] <= anode_y[1])]
    cathode = mesh.right
    pl, pr = float(neutral_potential(-a, delta)), float(neutral_potential(b, delta))
    dn = np.r_[anode, cathode]
    dv = np.r_[np.full(len(anode), pl), np.full(len(cathode), pr)]
    eq = solve_equilibrium(mesh, 1.0, c, delta, dn, dv)
    prm = DDParams(delta=delta, tau_n=1e2, tau_p=1e2, recombination=srh)
    dd = QFDriftDiffusion(mesh, c, prm, float(np.log(delta)))
    return mesh, c, eq, dd, anode, cathode, pl, pr, dn


def qf_bias(dd, mesh, X, anode, cathode, pl, pr, v):
    dn = np.r_[anode, cathode]
    dv = np.vstack(
        [np.tile((pl + v, v, v), (len(anode), 1)), np.tile((pr, 0.0, 0.0), (len(cathode), 1))]
    )
    return dd.solve(X, dn, dv)


def run_r1():
    out = {"poisson": [], "qf_dd": []}
    for N in (16, 32, 64, 128, 192, 256):
        g = np.linspace(0, 1, N + 1)
        mesh = TensorMesh(g, g)
        xy = mesh.node_xy
        u, f, _ = p1_fields(xy[:, 0], xy[:, 1])
        bn = mesh.boundary()
        eps = mesh.cell_coeff(eps_fun)
        t0 = time.perf_counter()
        solve_poisson(mesh, eps, f, bn, u[bn])
        out["poisson"].append(
            {
                "nodes": mesh.n_nodes,
                "N": N,
                "seconds": time.perf_counter() - t0,
                "peak_rss_mib": rss_mib(),
            }
        )
    for N in (16, 32, 48, 64, 96):
        mesh, c, eq, dd, an, ca, pl, pr, dn = diode_2d(N, N)
        X0 = pack(eq.psi, np.zeros(mesh.n_nodes), np.zeros(mesh.n_nodes))
        t0 = time.perf_counter()
        s = qf_bias(dd, mesh, X0, an, ca, pl, pr, 4.0)
        out["qf_dd"].append(
            {
                "nodes": mesh.n_nodes,
                "N": N,
                "V": 4.0,
                "seconds": time.perf_counter() - t0,
                "iterations": s.report.iterations,
                "converged": bool(s.report.converged),
                "peak_rss_mib": rss_mib(),
            }
        )
    return out


# --------------------------------------------------------------------- R2
def run_r2(figdir=None, nx=80, ny=40, v_end=14.0, delta=1e-3):
    mesh, c, eq, dd, an, ca, pl, pr, dn = diode_2d(nx, ny, delta=delta, anode_y=(7.5, 10.0))
    X = pack(eq.psi, np.zeros(mesh.n_nodes), np.zeros(mesh.n_nodes))
    steps = {}
    s = None
    for v in np.arange(0.0, v_end + 1e-9, 2.0):
        s = qf_bias(dd, mesh, X, an, ca, pl, pr, float(v))
        X = s.X
        steps[f"{v:g}"] = {
            "converged": bool(s.report.converged),
            "iterations": s.report.iterations,
            "I_anode": dd.contact_current(s, an),
            "I_cathode": -dd.contact_current(s, ca),
        }
    # current density profile in x-direction across the junction plane (x = 0)
    k0 = int(np.argmin(np.abs(mesh.x)))
    jx_profile = []
    for j in range(ny + 1):
        e = np.flatnonzero((mesh.edge_dir == 0) & (mesh.ea == mesh.node(k0, j)))
        jx_profile.append(float(s.j[e].sum() / mesh.w[e].sum()))
    jx = np.array(jx_profile)
    res = {
        "mesh": [nx, ny],
        "delta": delta,
        "anode": "left edge, y in [7.5, 10] of [0, 10]; cathode: full right edge",
        "steps": steps,
        "crowding_ratio_jx_max_over_min_at_x0": float(jx.max() / max(jx.min(), 1e-300)),
        "jx_at_x0_top_over_bottom": float(jx[-1] / jx[0]),
    }
    if figdir is not None:
        try:
            res["figure"] = _r2_figure(mesh, dd, s, figdir)
        except Exception as exc:  # figure is optional
            res["figure"] = f"skipped: {type(exc).__name__}: {exc}"
    return res


def _r2_figure(mesh, dd, s, figdir):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    nx, ny = mesh.nx, mesh.ny
    shape = (nx + 1, ny + 1)
    psi = s.psi.reshape(shape)
    Ex, Ey = np.gradient(-psi, mesh.x, mesh.y)
    dens = s.j / mesh.w
    jx = np.zeros(mesh.n_nodes)
    jy = np.zeros(mesh.n_nodes)
    cnt_x = np.zeros(mesh.n_nodes)
    cnt_y = np.zeros(mesh.n_nodes)
    for d, jv, cnt in ((0, jx, cnt_x), (1, jy, cnt_y)):
        sel = mesh.edge_dir == d
        for nd in (mesh.ea[sel], mesh.eb[sel]):
            np.add.at(jv, nd, dens[sel])
            np.add.at(cnt, nd, 1.0)
    jx, jy = (jx / np.maximum(cnt_x, 1)).reshape(shape), (jy / np.maximum(cnt_y, 1)).reshape(shape)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    X, Y = np.meshgrid(mesh.x, mesh.y, indexing="ij")
    m = ax[0].pcolormesh(X, Y, np.hypot(Ex, Ey), shading="gouraud", cmap="viridis")
    fig.colorbar(m, ax=ax[0], label="|E| (scaled)")
    ax[0].set_title("Field magnitude")
    jm = np.hypot(jx, jy)
    m = ax[1].pcolormesh(
        X, Y, np.log10(jm + 1e-300), shading="gouraud", cmap="magma", vmin=np.log10(jm.max()) - 4
    )
    fig.colorbar(m, ax=ax[1], label="log10 |j|")
    ax[1].streamplot(mesh.x, mesh.y, jx.T, jy.T, color="w", linewidth=0.6, density=1.2)
    ax[1].set_title("Current density and streamlines")
    for a in ax:
        a.set_xlabel("x")
        a.set_ylabel("y")
    Path(figdir).mkdir(parents=True, exist_ok=True)
    p = Path(figdir) / "stage8_current_crowding.png"
    fig.savefig(p, dpi=110)
    plt.close(fig)
    return str(p)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=REPO / "results" / "stage8")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "manifest.json").write_text(
        json.dumps(build_manifest(REPO, sys.argv, [0, 8]), indent=2) + "\n"
    )
    T = {}

    def timed(name, fn, *a, **k):
        t0 = time.perf_counter()
        r = fn(*a, **k)
        T[name] = time.perf_counter() - t0
        print(f"{name}: {T[name]:.1f}s", flush=True)
        return r

    res = {
        "P1": timed("P1", run_p1),
        "P2": timed("P2", run_p2),
        "P3": timed("P3", run_p3),
        "P4": timed("P4", run_p4),
        "P5": timed("P5", run_p5),
        "P6": timed("P6", run_p6),
        "P7": timed("P7", run_p7),
        "P8": timed("P8", run_p8),
        "P9": timed("P9", run_p9),
        "R1": timed("R1", run_r1),
        "R2": timed("R2", run_r2, REPO / "figures"),
    }
    p7 = res["P7"]
    crit = {
        "P1": last3_ok(res["P1"]["orders"], 1.9, 2.1),
        "P2": all(v["max_abs_err"] <= 1e-12 for v in res["P2"].values()),
        "P3": last3_ok(res["P3"]["orders"], 1.8, 2.2),
        "P4": res["P4"]["relative_balance"] <= 1e-12,
        "P5": last3_ok(res["P5"]["orders"], 1.9, 2.1)
        and all(v["converged"] for v in res["P5"]["newton"].values()),
        "P6": res["P6"]["max_row_diff"] <= 1e-11 and res["P6"]["newton_converged"],
        "P7": all(
            v["converged_1d"]
            and v["converged_2d"]
            and max(v["max_abs_diff_psi"], v["max_abs_diff_phi_n"], v["max_abs_diff_phi_p"])
            <= 1e-10
            and (not v["gated_current"] or v["j_rel_diff_mean_cut"] <= 1e-9)
            for v in p7.values()
        ),
        "P8": max(v["fd_1e-7"] for v in res["P8"].values()) <= 1e-6,
    }
    r9 = res["P9"]
    edge = r9["max_rel_diff_resolved_jn_first_jp_last_V1_110"]
    # P9 as amended (spec section 4): resolvable currents only
    crit["P9"] = (
        r9["all_converged"]
        and edge is not None
        and edge <= 1e-12
        and r9["max_rel_diff_abs_j_ge_1e-6"] <= 1e-12
    )
    crit = {k: bool(v) for k, v in crit.items()}
    # original wording (j_mean at every bias), reported only
    res["P9"]["literal_v1_pass"] = bool(r9["max_rel_diff_all"] <= 1e-12)
    summary = {
        "criteria": crit,
        "gate_passed": all(v for k, v in crit.items() if not k.endswith("_alt")),
        "values": {
            "P1_orders": res["P1"]["orders"],
            "P2": res["P2"],
            "P3_orders": res["P3"]["orders"],
            "P4_relative_balance": res["P4"]["relative_balance"],
            "P5_orders": res["P5"]["orders"],
            "P6_max_row_diff": res["P6"]["max_row_diff"],
            "P7": {str(k): v for k, v in p7.items()},
            "P8_worst_fd1e-7": max(v["fd_1e-7"] for v in res["P8"].values()),
            "P8_worst_fd1e-5": max(v["fd_1e-5"] for v in res["P8"].values()),
            "P9": res["P9"],
        },
        "seconds": T,
        "peak_rss_mib": rss_mib(),
    }
    (args.out / "results.json").write_text(json.dumps(res, indent=2, default=float) + "\n")
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    print(json.dumps(summary, indent=2, default=float))


if __name__ == "__main__":
    main()
