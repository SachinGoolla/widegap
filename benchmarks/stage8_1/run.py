"""Stage-8.1 benchmark (docs/benchmarks/stage8_1_mosfet.md): G1-G8 (gating), R1-R2 (reported).

Cases: C1 1D MOS capacitor (exact reference), C2 its 2D extrusion, C3 symmetric MOSFET
half-cell.  Writes <out>/{results,summary,manifest}.json (default results/stage8_1).
"""

from __future__ import annotations

import argparse
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

from widegap import mos_reference as ref
from widegap.fv.drift_diffusion import DDParams
from widegap.fv2d.mesh import TensorMesh
from widegap.fv2d.mos import Dirichlet, MOSDriftDiffusion, continue_bias
from widegap.fv2d.qf_dd import EPS_MACH
from widegap.fv2d.regions import Regions
from widegap.manifest import build_manifest

REPO = Path(__file__).resolve().parents[2]
DELTA = 1.0e-24
LOGD = float(np.log(DELTA))
EPS_R = 0.38  # eps_ox / eps_s (study parameter)
DP = DDParams(delta=DELTA, recombination=False)

# ---- C1 / C2 geometry (spec section 3)
A_SUB, T_OX, DEPTH, VFB_C1 = 10.0, 0.5, 30.0, 0.0
C1_N, C1_RATIO, C1_NOX = 40, 1.2, 4  # base semiconductor cells (geometric, finest at surface)
C1_LEVELS = (1, 2)  # bisection levels of the two Richardson meshes
PSI_F = abs(ref.bulk(A_SUB, DELTA)[0])
PSI_S_TARGETS = [-5, -3, -1, -0.3, -0.04, -0.01, 0.0, 0.01, 0.04, 0.3, 1, 2, 5, 10, 20, 40, 60]
PSI_S_TARGETS += [80, 100, 110, 2 * PSI_F, 2 * PSI_F + 2, 2 * PSI_F + 5, 2 * PSI_F + 10]

# ---- C3 geometry (all in units of the Debye length used in scaling; recorded in results)
C3 = {
    "x_breaks": [0.0, 3.0, 4.0, 8.0, 9.0, 12.0],
    "x_cells_base": [3, 2, 4, 2, 3],
    "y_semiconductor": [-8.0, 0.0],
    "well_depth": 2.0,
    "well_x": [[0.0, 4.0], [8.0, 12.0]],
    "channel_x": [4.0, 8.0],
    "channel_length": 4.0,
    "gate_x": [3.0, 9.0],
    "t_ox": T_OX,
    "oxide_over": "whole top (field oxide outside the gate), gate Dirichlet on x in [3, 9]",
    "a_body": 10.0,
    "b_wells": 100.0,
    "y_cells_base": {
        "body[-8,-2] ratio 0.7": 4,
        "well[-2,0] ratio 0.6 toward surface": 8,
        "oxide": 2,
    },
    "contacts": "source x=0 and drain x=12 on y in (-2, 0]; body: bottom y=-8; gate top y=t_ox",
    "mu_n": 1.0,
    "mu_p": 1.0,
    "W": 1.0,
    "eps_ratio": EPS_R,
}
V_DS = 2.0
V_ON = 200.0
TRANSFER_VG = [0, 60, 100, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 240]
LEVELS_G7 = (0, 1, 2)  # h, h/2, h/4


def rss_mib():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def refine(x, level):
    x = np.asarray(x, dtype=float)
    for _ in range(level):
        x = np.sort(np.r_[x, 0.5 * (x[:-1] + x[1:])])
    return x


def seg(lo, hi, n, q):
    h = q ** np.arange(n)
    h = h * (hi - lo) / h.sum()
    out = lo + np.r_[0.0, np.cumsum(h)]
    out[-1] = hi
    return out


# ============================================================ C1 / C2 / G1-G3
def c1_coords(level):
    h = C1_RATIO ** np.arange(C1_N)[::-1]
    h = h * DEPTH / h.sum()
    ys = np.r_[-DEPTH, -DEPTH + np.cumsum(h)]
    ys[-1] = 0.0
    y = np.r_[ys, np.linspace(0.0, T_OX, C1_NOX + 1)[1:]]
    return refine(y, level)


def make_cap(mesh):
    reg = Regions.from_function(mesh, lambda *c: c[-1] < 0)
    c = np.where(reg.semi_node, -A_SUB, 0.0)
    return reg, c, MOSDriftDiffusion(reg, c, DP, LOGD, eps_ratio=EPS_R)


def cap_sweep(mesh, reg, c, mod, bottom, top, pss):
    """Equilibrium solutions at gate voltages from the exact reference for each psi_s target."""
    vg = {p: ref.gate_voltage(p, A_SUB, DELTA, T_OX, EPS_R, VFB_C1) for p in pss}

    def bc_of(p):
        return Dirichlet().ohmic(c, DELTA, bottom, 0.0).gate(top, p[0], VFB_C1)

    X0 = mod.initial_guess(bc_of([vg[0.0]]))
    s0 = mod.solve_bc(X0, mod.equilibrium_bc(bc_of([vg[0.0]])))
    sols = {0.0: (s0, s0.report.iterations)}
    up = [p for p in pss if p > 0]
    down = sorted([p for p in pss if p < 0], reverse=True)
    for branch in (up, down):
        out = continue_bias(
            mod, s0.X, bc_of, [vg[0.0]], [[vg[p]] for p in branch], equilibrium=True
        )
        for p, (_, s, it, _) in zip(branch, out, strict=False):
            sols[p] = (s, it)
    return vg, sols


def gate_flux(mod, sol, top):
    """(L psi) summed over the gate nodes = eps_ox * oxide field integrated at the gate."""
    f = mod.G * (sol.psi[mod.a] - sol.psi[mod.b])
    acc = np.bincount(mod.a, f, minlength=mod.mesh.n_nodes) - np.bincount(
        mod.b, f, minlength=mod.mesh.n_nodes
    )
    return float(acc[top].sum()), acc


def run_c1():
    lv = {}
    for level in C1_LEVELS:
        y = c1_coords(level)
        m = TensorMesh(y)
        reg, c, mod = make_cap(m)
        isurf = int(np.argmin(np.abs(y)))
        vg, sols = cap_sweep(m, reg, c, mod, [0], [m.n_nodes - 1], PSI_S_TARGETS)
        lv[level] = (m, reg, c, mod, isurf, vg, sols)
    psi_b = ref.bulk(A_SUB, DELTA)[0]
    rows, worst_rel, all_ok, conv = [], 0.0, True, True
    g2_rows, g2_ok, worst_g2, worst_g2_floor = [], True, 0.0, 0.0
    for p in PSI_S_TARGETS:
        v = lv[C1_LEVELS[0]][5][p]
        exact = ref.psi_s_of_vg(v, A_SUB, DELTA, T_OX, EPS_R, VFB_C1)  # independent root find
        num = []
        for level in C1_LEVELS:
            m, reg, c, mod, isurf, vg, sols = lv[level]
            s = sols[p][0]
            conv &= bool(s.report.converged)
            num.append(float(s.psi[isurf] - psi_b))
        rich = (4.0 * num[1] - num[0]) / 3.0  # second-order Richardson (h, h/2)
        err = abs(rich - exact)
        small = abs(exact) < 0.05
        ok = err <= 1e-4 if small else err <= 2e-3 * abs(exact)
        all_ok &= ok
        if not small:
            worst_rel = max(worst_rel, err / abs(exact))
        rows.append(
            {"psi_s_target": p, "V_G": v, "exact": exact, "num_h": num[0], "num_h2": num[1],
             "richardson": rich, "abs_err": err, "ok": bool(ok)}
        )  # fmt: skip
        # G2 on the finer mesh
        m, reg, c, mod, isurf, vg, sols = lv[C1_LEVELS[-1]]
        s = sols[p][0]
        top = [m.n_nodes - 1]
        fl, acc = gate_flux(mod, s, top)
        q_num = mod.space_charge(s)
        gerr = abs(fl + q_num)
        e_top = np.flatnonzero(mod.b == top[0])
        gscale = float(
            np.sum(mod.vol * (s.p + s.n + np.abs(mod.c)))
            + np.sum(mod.G[e_top] * (abs(s.psi[mod.b[e_top]]) + abs(s.psi[mod.a[e_top]])))
        )
        floor = 1e2 * EPS_MACH * gscale
        q_exact = ref.Qs(exact, A_SUB, DELTA)
        ok2 = gerr <= 1e-10 * abs(q_num) + floor
        g2_ok &= ok2
        worst_g2 = max(worst_g2, gerr / max(abs(q_num), 1e-300))
        worst_g2_floor = max(worst_g2_floor, gerr / floor)
        g2_rows.append(
            {"psi_s_target": p, "gate_flux": fl, "Q_numeric": q_num, "Q_exact": q_exact,
             "abs_err": gerr, "rel_err": gerr / max(abs(q_num), 1e-300), "roundoff_floor": floor,
             "ok": bool(ok2)}
        )  # fmt: skip
    res1 = {
        "mesh_nodes": {str(k): v[0].n_nodes for k, v in lv.items()},
        "surface_cell": {str(k): float(np.diff(v[0].x)[np.argmin(np.abs(v[0].x[:-1]))- 1]) for k, v in lv.items()},
        "all_converged": conv,
        "newton_its_max": max(it for v in lv.values() for (_, it) in v[6].values()),
        "worst_rel_err_nonsmall": worst_rel,
        "points": rows,
        "ok": bool(all_ok and conv),
    }  # fmt: skip
    res2 = {
        "worst_rel_err": worst_g2,
        "worst_err_over_floor": worst_g2_floor,
        "points": g2_rows,
        "ok": bool(g2_ok),
    }
    return res1, res2, lv[C1_LEVELS[0]]


def run_c2(c1_level):
    """C2: 3-column 2D extrusion of the level-1 C1 mesh; every column vs the 1D solution."""
    m1, reg1, c1, mod1, isurf1, vg1, sols1 = c1_level
    y = m1.x
    m2 = TensorMesh(np.array([0.0, 1.0, 2.0, 3.0]), y)
    reg2, c2, mod2 = make_cap(m2)
    ny1 = m2.ny + 1
    bottom = [m2.node(i, 0) for i in range(m2.nx + 1)]
    top = [m2.node(i, m2.ny) for i in range(m2.nx + 1)]
    vg2, sols2 = cap_sweep(m2, reg2, c2, mod2, bottom, top, PSI_S_TARGETS)
    worst, conv = 0.0, True
    for p in PSI_S_TARGETS:
        s2 = sols2[p][0]
        conv &= bool(s2.report.converged)
        for i in range(m2.nx + 1):
            d = np.abs(s2.psi[i * ny1 : (i + 1) * ny1] - sols1[p][0].psi).max()
            worst = max(worst, float(d))
    return {
        "mesh": [m2.nx, m2.ny],
        "all_converged": conv,
        "max_abs_diff_psi": worst,
        "ok": bool(conv and worst <= 1e-10),
    }


# ============================================================ C3 MOSFET
def build_c3(level):
    xs = [0.0]
    for lo, hi, n in zip(C3["x_breaks"][:-1], C3["x_breaks"][1:], C3["x_cells_base"], strict=True):
        xs += list(seg(lo, hi, n, 1.0)[1:])
    ys = np.r_[seg(-8.0, -2.0, 4, 0.7), seg(-2.0, 0.0, 8, 0.6)[1:], np.linspace(0, T_OX, 3)[1:]]
    x, y = refine(xs, level), refine(ys, level)
    m = TensorMesh(x, y)
    reg = Regions.from_function(m, lambda xc, yc: yc < 0)
    cc = m.cell_centres
    well = ((cc[:, 0] < 4.0) | (cc[:, 0] > 8.0)) & (cc[:, 1] > -C3["well_depth"])
    c = reg.semi_average(np.where(well, C3["b_wells"], -C3["a_body"]))
    mod = MOSDriftDiffusion(reg, c, DP, LOGD, eps_ratio=EPS_R)
    xy = m.node_xy
    sd = (xy[:, 1] > -C3["well_depth"] + 1e-9) & (xy[:, 1] <= 0.0)
    nodes = {
        "source": np.flatnonzero((xy[:, 0] == 0.0) & sd),
        "drain": np.flatnonzero((xy[:, 0] == 12.0) & sd),
        "body": np.flatnonzero(xy[:, 1] == -8.0),
        "gate": np.flatnonzero((xy[:, 1] == T_OX) & (xy[:, 0] >= 3 - 1e-9) & (xy[:, 0] <= 9 + 1e-9)),
    }  # fmt: skip
    return m, reg, c, mod, nodes


class Device:
    def __init__(self, level):
        self.level = level
        self.m, self.reg, self.c, self.mod, self.nodes = build_c3(level)
        self.vfb = -ref.bulk(C3["a_body"], DELTA)[0]  # flat band at V_G = 0 (study value)
        self.cache = {}

    def bc_of(self, p):
        vg, vl, vr = p
        n = self.nodes
        return (
            Dirichlet()
            .ohmic(self.c, DELTA, n["source"], vl)
            .ohmic(self.c, DELTA, n["drain"], vr)
            .ohmic(self.c, DELTA, n["body"], 0.0)
            .gate(n["gate"], vg, self.vfb)
        )

    def equilibrium(self):
        bc = self.bc_of([0.0, 0.0, 0.0])
        return self.mod.solve_bc(self.mod.initial_guess(bc), self.mod.equilibrium_bc(bc))

    def currents(self, sol):
        n, mod = self.nodes, self.mod
        return {
            "I_source": mod.contact_current(sol, n["source"]),
            "I_drain": mod.contact_current(sol, n["drain"]),
            "I_body": mod.contact_current(sol, n["body"]),
        }

    def march(self, X, p_from, targets):
        return continue_bias(self.mod, X, self.bc_of, p_from, targets)


def bias_sweep(dev, vl, vr, vgs, X_eq):
    """V_G = 0 with the drain/source bias applied, then the V_G sweep (continuation)."""
    out = dev.march(X_eq, [0.0, 0.0, 0.0], [[0.0, vl, vr]] + [[vg, vl, vr] for vg in vgs])
    pts = []
    for tgt, s, it, ns in out[1:]:
        cur = dev.currents(s)
        pts.append(
            {"V_G": float(tgt[0]), "V_left": vl, "V_right": vr, "converged": bool(s.report.converged),
             "newton_its": it, "solves": ns, "final_row_ratio": s.report.residual_history[-1], **cur}
        )  # fmt: skip
        pts[-1]["_sol"] = s
    return out[0][1].report.converged and out[0][1], pts


def kcl(pt):
    return abs(pt["I_source"] + pt["I_drain"] + pt["I_body"])


def run_c3_transfer(dev):
    s_eq = dev.equilibrium()
    first, pts = bias_sweep(dev, 0.0, V_DS, TRANSFER_VG, s_eq.X)
    return s_eq, first, pts


def run_g6(dev, s_eq):
    rows = []
    for v in (1.0, 2.0):
        for g in (170.0, 200.0):
            _, A = bias_sweep(dev, 0.0, v, [g], s_eq.X)
            _, B = bias_sweep(dev, v, 0.0, [g], s_eq.X)
            iA, iB = A[0]["I_drain"], B[0]["I_drain"]  # drain = right contact in both
            rel = abs(iA + iB) / abs(iA)
            rows.append(
                {"V_DS": v, "V_G": g, "converged": bool(A[0]["converged"] and B[0]["converged"]),
                 "I_D(+V_DS)": iA, "I_D(-V_DS, swapped)": iB, "I_body_A": A[0]["I_body"],
                 "rel_antisymmetry": rel, "gated": abs(iA) >= 1e-6}
            )  # fmt: skip
    ok = all(r["converged"] and (r["rel_antisymmetry"] <= 1e-8 or not r["gated"]) for r in rows)
    return {"points": rows, "ok": bool(ok)}


def run_g7():
    vals, info = [], []
    for level in LEVELS_G7:
        dev = Device(level)
        s_eq = dev.equilibrium()
        _, pts = bias_sweep(dev, 0.0, V_DS, [100.0, 150.0, V_ON], s_eq.X)
        pt = pts[-1]
        vals.append(pt["I_drain"])
        info.append(
            {"level": level, "cells": [dev.m.nx, dev.m.ny], "nodes": dev.m.n_nodes,
             "converged": bool(pt["converged"]), "kcl_abs": kcl(pt), "I_drain": pt["I_drain"],
             "I_source": pt["I_source"], "I_body": pt["I_body"],
             "surface_cell_height": float(dev.m.y[dev.m.y.searchsorted(0.0)] - dev.m.y[dev.m.y.searchsorted(0.0) - 1])}
        )  # fmt: skip
    i0, i1, i2 = vals
    d01, d12 = i0 - i1, i1 - i2
    order = float(np.log2(d01 / d12)) if d01 * d12 > 0 else float("nan")
    ok_order = bool(np.isfinite(order) and 0.8 <= order <= 2.5)
    if np.isfinite(order):
        i_ex = i2 + (i2 - i1) / (2.0**order - 1.0)
        rel = abs(i2 - i_ex) / abs(i_ex)
    else:
        i_ex, rel = float("nan"), float("nan")
    return {
        "bias": {"V_GS": V_ON, "V_DS": V_DS},
        "levels": info,
        "observed_order": order,
        "richardson_I_D": i_ex,
        "finest_rel_dev_from_richardson": rel,
        "ok": bool(ok_order and rel <= 0.02 and all(r["converged"] for r in info)),
    }


def run_g8(dev, s_eq):
    mod = dev.mod
    res = {}
    for name, vg in (("subthreshold_V_G=100", 100.0), ("on_state_V_G=200", V_ON)):
        out = dev.march(s_eq.X, [0.0, 0.0, 0.0], [[0.0, 0.0, V_DS], [vg, 0.0, V_DS]])
        s = out[-1][1]
        bc = dev.bc_of([vg, 0.0, V_DS])
        rows, _ = mod._all_bc(bc)
        Jan = mod.jacobian_bc(s.X, rows).toarray()
        entry = {"converged": bool(s.report.converged)}
        for h in (1e-5, 1e-6):
            Jfd = mod.fd_jacobian(s.X, bc, h=h).toarray()
            rowmax = np.abs(Jan).max(axis=1)
            err = np.abs(Jfd - Jan).max(axis=1) / rowmax
            entry[f"worst_row_err_h{h:g}"] = float(err.max())
            entry[f"worst_row_index_h{h:g}"] = [int(err.argmax() // 3), int(err.argmax() % 3)]
        res[name] = entry
    res["worst"] = max(
        v[k] for v in res.values() for k in v if k.startswith("worst_row_err_h1e-05")
    )
    res["ok"] = bool(
        res["worst"] <= 1e-6 and all(v["converged"] for k, v in res.items() if isinstance(v, dict))
    )
    return res


def run_r1(pts):
    vg = np.array([p["V_G"] for p in pts])
    idr = np.array([p["I_drain"] for p in pts])
    gm = np.diff(idr) / np.diff(vg)
    k = int(np.argmax(gm))
    vmid = 0.5 * (vg[k] + vg[k + 1])
    imid = 0.5 * (idr[k] + idr[k + 1])
    vt = float(vmid - imid / gm[k])
    cox = EPS_R / T_OX
    pref = C3["mu_n"] * C3["W"] / C3["channel_length"] * cox * V_DS
    rows = [
        {"V_G": float(v), "I_D": float(i), "I_gca": float(pref * max(v - vt, 0.0)),
         "ratio": float(i / (pref * (v - vt))) if v > vt + 1e-9 else None}
        for v, i in zip(vg, idr, strict=True)
    ]  # fmt: skip
    return {
        "V_T_linear_extrapolation_max_gm": vt,
        "V_T_minus_VDS_over_2": vt - V_DS / 2,
        "gm_max": float(gm[k]),
        "gm_gca": float(pref),
        "gm_max_over_gca": float(gm[k] / pref),
        "C_ox": cox,
        "W_over_L": C3["W"] / C3["channel_length"],
        "note": "L = source/drain well spacing (4); gate overlaps each well by 1; no mobility degradation",
        "curve": rows,
    }


def r2_figure(dev, sol, figdir):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    m, reg = dev.m, dev.reg
    shape = (m.nx + 1, m.ny + 1)
    X, Y = np.meshgrid(m.x, m.y, indexing="ij")
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6), constrained_layout=True)
    mp = ax[0].pcolormesh(X, Y, sol.psi.reshape(shape), shading="gouraud", cmap="viridis")
    fig.colorbar(mp, ax=ax[0], label="psi (scaled, V_T)")
    ax[0].set_title(f"Potential, V_GS={V_ON:g}, V_DS={V_DS:g}")
    with np.errstate(divide="ignore"):
        ln = np.where(reg.semi_node, np.log10(np.maximum(sol.n, 1e-300)), np.nan)
    mp = ax[1].pcolormesh(X, Y, ln.reshape(shape), shading="gouraud", cmap="magma", vmin=-2)
    fig.colorbar(mp, ax=ax[1], label="log10 n (scaled)")
    ax[1].set_title("Electron density (semiconductor)")
    for a in ax:
        a.axhline(0.0, color="w", lw=0.6)
        a.set_xlabel("x")
        a.set_ylabel("y")
    Path(figdir).mkdir(parents=True, exist_ok=True)
    p = Path(figdir) / "stage8_1_mosfet_on_state.png"
    fig.savefig(p, dpi=110)
    plt.close(fig)
    return str(p)


def main() -> None:
    default_out = REPO / "results" / "stage8_1"
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=default_out)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    figdir = REPO / "figures" if args.out == default_out else args.out / "figures"
    (args.out / "manifest.json").write_text(
        json.dumps(build_manifest(REPO, sys.argv, []), indent=2) + "\n"
    )
    T = {}

    def timed(name, fn, *a):
        t0 = time.perf_counter()
        r = fn(*a)
        T[name] = time.perf_counter() - t0
        print(f"{name}: {T[name]:.1f}s", flush=True)
        return r

    g1, g2, lvl1 = timed("C1", run_c1)
    g3 = timed("C2", run_c2, lvl1)
    dev = Device(2)
    s_eq, first, pts = timed("C3_transfer", run_c3_transfer, dev)
    dev1 = Device(1)
    g6 = timed("G6", run_g6, dev1, dev1.equilibrium())
    g7 = timed("G7", run_g7)
    dev0 = Device(0)
    g8 = timed("G8", run_g8, dev0, dev0.equilibrium())
    r1 = run_r1(pts)
    on = [p for p in pts if p["V_G"] == V_ON][0]
    try:
        fig = r2_figure(dev, on["_sol"], figdir)
    except Exception as exc:  # figure is optional
        fig = f"skipped: {type(exc).__name__}: {exc}"
    transfer = [{k: v for k, v in p.items() if k != "_sol"} for p in pts]

    # G4: every C3 solve converged (transfer, G6, G7, G8)
    g4_all = (
        bool(first)
        and all(p["converged"] for p in transfer)
        and all(r["converged"] for r in g6["points"])
        and all(r["converged"] for r in g7["levels"])
        and all(v["converged"] for k, v in g8.items() if isinstance(v, dict))
    )
    # G5: KCL at every C3 point with |I_D| >= 1e-6
    gated5 = [(kcl(p), abs(p["I_drain"])) for p in transfer if abs(p["I_drain"]) >= 1e-6]
    gated5 += [(r["kcl_abs"], abs(r["I_drain"])) for r in g7["levels"] if abs(r["I_drain"]) >= 1e-6]
    worst5 = max(a / b for a, b in gated5) if gated5 else float("nan")
    g5 = {"n_points": len(gated5), "worst_rel_kcl": worst5, "ok": bool(gated5 and worst5 <= 1e-8)}
    g4 = {"all_converged": g4_all, "ok": bool(g4_all)}

    res = {
        "C3_geometry": C3,
        "G1": g1, "G2": g2, "G3": g3, "G4": g4, "G5": g5, "G6": g6, "G7": g7, "G8": g8,
        "C3_transfer_level2": {"mesh": [dev.m.nx, dev.m.ny], "V_DS": V_DS,
                               "V_FB": dev.vfb, "points": transfer},
        "R1": r1, "R2": {"figure": fig},
    }  # fmt: skip
    crit = {k: bool(res[k]["ok"]) for k in ("G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8")}
    summary = {
        "criteria": crit,
        "gate_passed": all(crit.values()),
        "values": {
            "G1_worst_rel_err": g1["worst_rel_err_nonsmall"],
            "G1_worst_abs_err": max(p["abs_err"] for p in g1["points"]),
            "G2_worst_rel_err": g2["worst_rel_err"],
            "G2_worst_err_over_roundoff_floor": g2["worst_err_over_floor"],
            "G3_max_abs_diff_psi": g3["max_abs_diff_psi"],
            "G5": g5,
            "G6_worst_rel": max(r["rel_antisymmetry"] for r in g6["points"]),
            "G7": {k: g7[k] for k in ("observed_order", "richardson_I_D",
                                      "finest_rel_dev_from_richardson")},
            "G8_worst_row_err": g8["worst"],
            "R1_V_T": r1["V_T_linear_extrapolation_max_gm"],
            "R1_gm_max_over_gca": r1["gm_max_over_gca"],
        },
        "seconds": T,
        "peak_rss_mib": rss_mib(),
    }  # fmt: skip
    (args.out / "results.json").write_text(json.dumps(res, indent=2, default=float) + "\n")
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    print(json.dumps(summary, indent=2, default=float))


if __name__ == "__main__":
    main()
