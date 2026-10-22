"""Stage-9a benchmark (docs/benchmarks/stage9a_dd_physics.md): A1-A7 (gating), R1 (reported).

Incomplete ionization, Caughey-Thomas mobility and interface charge on the stage-8.1 MOS
drift-diffusion core (``widegap.fv2d.physics``).  Writes <out>/{results,summary,manifest}.json
(default results/stage9a).
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
from widegap.fv.ionized_equilibrium import IonizedJunction
from widegap.fv2d.mesh import TensorMesh
from widegap.fv2d.mos import Dirichlet, MOSDriftDiffusion, continue_bias
from widegap.fv2d.physics import (
    InterfaceCharge,
    Ionization,
    Mobility,
    PhysicsMOSDD,
    ct_factor,
    gate_voltage_it,
    neutral_psi,
    psi_s_of_vg_it,
)
from widegap.fv2d.regions import Regions
from widegap.manifest import build_manifest
from widegap.materials.ionization import Level

REPO = Path(__file__).resolve().parents[2]
EPS_MACH = np.finfo(float).eps
DELTA = 1.0e-24
LOGD = float(np.log(DELTA))
EPS_R = 0.38
DP = DDParams(delta=DELTA, recombination=False)

# ---- study parameters (all scaled units; chosen to exercise each model)
LEV_D = (Level(1.0, -62.0),)  # MOSFET wells / A5: mildly incomplete
LEV_A = (Level(1.0, -62.0),)
ITF_STUDY = InterfaceCharge(q_f=3.0, dit=0.2, e1=-25.0, e2=15.0)  # A5
ITF_C3 = InterfaceCharge(q_f=5.0, dit=0.1, e1=-20.0, e2=20.0)  # A7 / R1
MOB_C3 = Mobility(vsat_n=0.5, vsat_p=0.5, beta_n=2.0, beta_p=1.0)
A_SUB, T_OX, DEPTH = 10.0, 0.5, 30.0

# ---- C3 geometry (stage-8.1, copied)
C3 = {
    "x_breaks": [0.0, 3.0, 4.0, 8.0, 9.0, 12.0],
    "x_cells_base": [3, 2, 4, 2, 3],
    "well_depth": 2.0,
    "t_ox": T_OX,
    "a_body": 10.0,
    "b_wells": 100.0,
    "mu_n": 1.0,
    "mu_p": 1.0,
    "eps_ratio": EPS_R,
}
V_DS = 2.0
TRANSFER_VG = [0, 60, 100, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 240]
OUTPUT_VDS = [0.5, 1, 2, 4, 8, 16]


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


def all_semi(mesh):
    return Regions.from_function(mesh, lambda *c: np.ones_like(c[0], dtype=bool))


# ============================================================ A1
def run_a1():
    """Equilibrium p-n junction, incomplete ionization: DD core (phi = 0) vs IonizedJunction."""
    lev = (Level(1.0, -58.5),)
    t = np.linspace(-1.0, 1.0, 241)
    x = 12.0 * np.sinh(3.0 * t) / np.sinh(3.0)
    mesh = TensorMesh(x)
    nd = np.where(x > 0, 10.0, 0.0)
    na = np.where(x < 0, 10.0, 0.0)
    nd[x == 0] = na[x == 0] = 5.0
    reg = all_semi(mesh)
    ion = Ionization(nd, na, lev, lev)
    dd = PhysicsMOSDD(reg, nd - na, DP, LOGD, ionization=ion)
    pl = neutral_psi(0.0, 10.0, DELTA, lev, lev)
    pr = neutral_psi(10.0, 0.0, DELTA, lev, lev)
    bc = dd.equilibrium_bc(Dirichlet().psi([0, mesh.n_nodes - 1], [pl, pr]))
    sol = dd.solve_bc(dd.initial_guess(bc), bc)
    jn = IonizedJunction(x, nd, na, DELTA, lev, lev)
    psi_ref, rep = jn.solve(sol.psi * 0.0 + np.where(x > 0, pr, pl), pl, pr)
    d = float(np.abs(sol.psi - psi_ref).max())
    complete = float(np.abs(nd - na).max())
    return {
        "nodes": mesh.n_nodes,
        "psi_range": [float(sol.psi.min()), float(sol.psi.max())],
        "max_abs_dpsi": d,
        "dd_converged": bool(sol.report.converged),
        "ref_converged": bool(rep.converged),
        "dd_iterations": sol.report.iterations,
        "ref_iterations": rep.iterations,
        "ionized_fraction_donor_bulk": float(
            ion.nd[-1] / ion.nd[-1] / (1.0 + np.exp(sol.psi[-1] + lev[0].log_gamma))
        ),
        "doping": complete,
        "ok": bool(d <= 1e-10 and sol.report.converged and rep.converged),
    }


# ============================================================ A2 (MMS)
MMS_DELTA = 0.2
MMS_LD = (Level(0.6, 0.3), Level(0.4, 1.0))
MMS_LA = (Level(1.0, -0.2),)


def _mms_sources(dim):
    import sympy as s

    xs = s.symbols("x y")[:dim]
    x = xs[0]
    y = xs[1] if dim == 2 else 0
    pi = s.pi
    psi = 0.5 + 0.8 * x - 0.4 * y + 0.3 * s.sin(pi * x) * s.cos(pi * y / 2) + 0.2 * x * y
    pn = 0.4 * x + 0.2 * s.sin(2 * pi * x) * (1 + 0.5 * y) + 0.15 * x * y
    pp = -0.3 * x + 0.25 * y + 0.1 * s.cos(pi * x) * s.sin(pi * y)
    nd = 3.0 + s.sin(pi * x) + 0.5 * y
    na = 2.0 + s.cos(pi * x) * (1 - 0.3 * y)
    dl = MMS_DELTA
    n = dl * s.exp(psi - pn)
    p = dl * s.exp(pp - psi)
    ndp = sum(lv.fraction * nd / (1 + s.exp(psi - pn + lv.log_gamma)) for lv in MMS_LD)
    nam = sum(lv.fraction * na / (1 + s.exp(pp - psi + lv.log_gamma)) for lv in MMS_LA)
    s_psi = -sum(s.diff(psi, c, 2) for c in xs) - (p - n + ndp - nam)
    s_n = sum(s.diff(s.diff(n, c) - n * s.diff(psi, c), c) for c in xs)
    s_p = sum(s.diff(-(s.diff(p, c) + p * s.diff(psi, c)), c) for c in xs)
    names = ("psi", "pn", "pp", "nd", "na", "s_psi", "s_n", "s_p")
    exprs = (psi, pn, pp, nd, na, s_psi, s_n, s_p)
    f = s.lambdify(xs, exprs, "numpy")

    def evaluate(*coords):
        vals = f(*coords)
        shape = np.broadcast(*coords).shape
        return {
            k: np.broadcast_to(np.asarray(v, dtype=float), shape)
            for k, v in zip(names, vals, strict=True)
        }

    return evaluate


def mms_solve(dim, ncell):
    ev = _mms_sources(dim)
    g = np.linspace(0.0, 1.0, ncell + 1)
    mesh = TensorMesh(g) if dim == 1 else TensorMesh(g, g)
    f = ev(*mesh.node_xy.T)
    reg = all_semi(mesh)
    dp = DDParams(delta=MMS_DELTA, recombination=False)
    ion = Ionization(f["nd"], f["na"], MMS_LD, MMS_LA)
    dd = PhysicsMOSDD(
        reg,
        f["nd"] - f["na"],
        dp,
        float(np.log(MMS_DELTA)),
        sources={"psi": f["s_psi"], "n": f["s_n"], "p": f["s_p"]},
        ionization=ion,
    )
    b = mesh.boundary()
    bc = Dirichlet().psi(b, f["psi"][b]).phi(b, f["pn"][b], f["pp"][b])
    X0 = np.zeros(3 * mesh.n_nodes)
    X0[0::3], X0[1::3], X0[2::3] = f["psi"], f["pn"], f["pp"]
    sol = dd.solve_bc(X0, bc)
    err = [
        float(np.abs(sol.psi - f["psi"]).max()),
        float(np.abs(sol.phi_n - f["pn"]).max()),
        float(np.abs(sol.phi_p - f["pp"]).max()),
    ]
    return err, bool(sol.report.converged)


def run_a2():
    out, ok = {}, True
    for dim, ns in ((1, (16, 32, 64, 128, 256)), (2, (8, 16, 32, 64, 128))):
        errs, conv = [], True
        for n in ns:
            e, c = mms_solve(dim, n)
            errs.append(e)
            conv &= c
        errs = np.array(errs)
        orders = np.log2(errs[:-1] / errs[1:])  # (len-1, 3)
        last3 = orders[-3:]
        ok_d = bool(conv and np.all((last3 >= 1.8) & (last3 <= 2.2)))
        ok &= ok_d
        out[f"{dim}D"] = {
            "cells": list(ns),
            "errors_psi_phin_phip": errs.tolist(),
            "orders": orders.tolist(),
            "last_three_orders": last3.tolist(),
            "converged": conv,
            "ok": ok_d,
        }
    out["ok"] = bool(ok)
    return out


# ============================================================ A3
def build_c3(level, ion=True, mobility=None, interface=None, base_class=False):
    xs = [0.0]
    for lo, hi, n in zip(C3["x_breaks"][:-1], C3["x_breaks"][1:], C3["x_cells_base"], strict=True):
        xs += list(seg(lo, hi, n, 1.0)[1:])
    ys = np.r_[seg(-8.0, -2.0, 4, 0.7), seg(-2.0, 0.0, 8, 0.6)[1:], np.linspace(0, T_OX, 3)[1:]]
    m = TensorMesh(refine(xs, level), refine(ys, level))
    reg = Regions.from_function(m, lambda xc, yc: yc < 0)
    cc = m.cell_centres
    well = ((cc[:, 0] < 4.0) | (cc[:, 0] > 8.0)) & (cc[:, 1] > -C3["well_depth"])
    c = reg.semi_average(np.where(well, C3["b_wells"], -C3["a_body"]))
    nd = reg.semi_average(np.where(well, C3["b_wells"], 0.0))
    na = reg.semi_average(np.where(well, 0.0, C3["a_body"]))
    if base_class:
        mod = MOSDriftDiffusion(reg, c, DP, LOGD, eps_ratio=EPS_R)
    else:
        io = Ionization(nd, na, LEV_D, LEV_A) if ion else None
        mod = PhysicsMOSDD(
            reg, c, DP, LOGD, eps_ratio=EPS_R, ionization=io, mobility=mobility, interface=interface
        )
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
    def __init__(self, level, **kw):
        self.level = level
        self.m, self.reg, self.c, self.mod, self.nodes = build_c3(level, **kw)
        if isinstance(self.mod, PhysicsMOSDD) and self.mod.ion is not None:
            X = self.mod.initial_guess()
            self.vfb = -float(X[0::3][self.nodes["body"][0]])
        else:
            self.vfb = -ref.bulk(C3["a_body"], DELTA)[0]

    def contact(self, bc, name, v):
        n, mod = self.nodes[name], self.mod
        if isinstance(mod, PhysicsMOSDD):
            return mod.ohmic(bc, n, v)
        return bc.ohmic(self.c, DELTA, n, v)

    def bc_of(self, p):
        vg, vl, vr = p
        bc = Dirichlet()
        bc = self.contact(bc, "source", vl)
        bc = self.contact(bc, "drain", vr)
        bc = self.contact(bc, "body", 0.0)
        return bc.gate(self.nodes["gate"], vg, self.vfb)

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

    def point(self, X_eq, vg, vl, vr):
        out = self.march(X_eq, [0.0, 0.0, 0.0], [[0.0, vl, vr], [vg, vl, vr]])
        return out[-1][1]


def run_a3():
    vg = 200.0
    res = {}
    cases = {
        "vs_stage8.1_class": (
            dict(base_class=True),
            dict(ion=False, mobility=Mobility(1e30, 1e30, 2.0, 2.0)),
        ),
        "all_other_models_on": (
            dict(ion=True, mobility=None, interface=ITF_C3),
            dict(ion=True, mobility=Mobility(1e30, 1e30, 2.0, 1.0), interface=ITF_C3),
        ),
    }
    ok = True
    for name, (kw_ref, kw_new) in cases.items():
        d_ref, d_new = Device(0, **kw_ref), Device(0, **kw_new)
        sr = d_ref.point(d_ref.equilibrium().X, vg, 0.0, V_DS)
        sn = d_new.point(d_new.equilibrium().X, vg, 0.0, V_DS)
        jr, jn = sr.j, sn.j
        gated = np.abs(jr) >= 1e-6
        rel = float(np.max(np.abs(jn[gated] - jr[gated]) / np.abs(jr[gated])))
        cur_r, cur_n = d_ref.currents(sr), d_new.currents(sn)
        rel_cur = max(abs(cur_n[k] - cur_r[k]) / abs(cur_r[k]) for k in ("I_drain", "I_source"))
        e = {
            "V_G": vg,
            "V_DS": V_DS,
            "converged": bool(sr.report.converged and sn.report.converged),
            "n_gated_edges": int(gated.sum()),
            "max_rel_edge_current_diff": rel,
            "max_rel_terminal_current_diff": float(rel_cur),
            "I_drain_ref": cur_r["I_drain"],
            "ok": bool(sr.report.converged and sn.report.converged and rel <= 1e-12),
        }
        ok &= e["ok"]
        res[name] = e
    res["ok"] = bool(ok)
    return res


# ============================================================ A4
BAR = dict(
    n_d=10.0, delta=1e-3, length=40.0, ncell=80, vsat_n=0.5, vsat_p=0.5, beta_n=2.0, beta_p=1.0
)
BAR_V_FACTORS = [0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]  # V / (vsat L / mu0)


def run_a4():
    B = BAR
    L, nd0, dl = B["length"], B["n_d"], B["delta"]
    mesh = TensorMesh(np.linspace(0.0, L, B["ncell"] + 1))
    reg = all_semi(mesh)
    c = np.full(mesh.n_nodes, nd0)
    dp = DDParams(delta=dl, recombination=False)
    mob = Mobility(B["vsat_n"], B["vsat_p"], B["beta_n"], B["beta_p"])
    dd = PhysicsMOSDD(reg, c, dp, float(np.log(dl)), mobility=mob)
    xc = 0.5 * (mesh.x[:-1] + mesh.x[1:])
    interior = (xc >= 5.0) & (xc <= L - 5.0)  # >= 5 L_D (L_D <= 1 in scaled units)
    rows, ok = [], True
    for fct in BAR_V_FACTORS:
        v = fct * B["vsat_n"] * L / dp.mu_n
        bc = Dirichlet().ohmic(c, dl, [0], 0.0).ohmic(c, dl, [mesh.n_nodes - 1], v)
        X0 = np.zeros(3 * mesh.n_nodes)
        psi_eq = np.arcsinh(nd0 / (2 * dl))
        X0[0::3] = psi_eq + v * mesh.x / L
        X0[1::3] = X0[2::3] = v * mesh.x / L
        sol = dd.solve_bc(X0, bc)
        F = v / L
        m, _ = ct_factor(np.array([v]), L, dp.mu_n, B["vsat_n"], B["beta_n"])
        j_closed = nd0 * dp.mu_n * float(m[0]) * F
        j_num = -sol.j[interior]  # current flows towards decreasing phi
        jn_num = -sol.jn[interior]
        rel = float(np.max(np.abs(j_num - j_closed)) / j_closed)
        rel_n = float(np.max(np.abs(jn_num - j_closed)) / j_closed)
        row = {
            "V": v,
            "V_over_vsatL": fct,
            "F": F,
            "mu_F_over_vsat": dp.mu_n * F / B["vsat_n"],
            "j_closed": j_closed,
            "j_num_mean": float(j_num.mean()),
            "max_rel_err_total_j": rel,
            "max_rel_err_jn_only": rel_n,
            "j_over_q_n_vsat": float(j_num.mean() / (nd0 * B["vsat_n"])),
            "converged": bool(sol.report.converged),
            "newton_its": sol.report.iterations,
        }
        rows.append(row)
        ok &= bool(sol.report.converged and rel <= 1e-6)
    sat = rows[-1]["j_over_q_n_vsat"]
    sat_expected = rows[-1]["mu_F_over_vsat"] * (1 + rows[-1]["mu_F_over_vsat"] ** B["beta_n"]) ** (
        -1.0 / B["beta_n"]
    )
    # saturation reproduced: numerical j / (q n vsat) equals the closed-form limit approach
    sat_ok = bool(abs(sat - sat_expected) <= 1e-6 * sat_expected and sat >= 0.99)
    return {
        "params": B,
        "points": rows,
        "j_over_qnvsat_at_max_V": sat,
        "j_over_qnvsat_closed_form": sat_expected,
        "saturation_ok": sat_ok,
        "ok": bool(ok and sat_ok),
    }


# ============================================================ A5
C1_N, C1_RATIO, C1_NOX = 40, 1.2, 4
C1_LEVELS = (1, 2)
PSI_F = abs(ref.bulk(A_SUB, DELTA)[0])
PSI_S_TARGETS = [-5, -3, -1, -0.3, -0.04, -0.01, 0.0, 0.01, 0.04, 0.3, 1, 2, 5, 10, 20, 40, 60]
PSI_S_TARGETS += [80, 100, 110, 2 * PSI_F, 2 * PSI_F + 2, 2 * PSI_F + 5, 2 * PSI_F + 10]


def c1_coords(level):
    h = C1_RATIO ** np.arange(C1_N)[::-1]
    h = h * DEPTH / h.sum()
    ys = np.r_[-DEPTH, -DEPTH + np.cumsum(h)]
    ys[-1] = 0.0
    y = np.r_[ys, np.linspace(0.0, T_OX, C1_NOX + 1)[1:]]
    return refine(y, level)


def run_a5():
    itf = ITF_STUDY
    psi_b = ref.bulk(A_SUB, DELTA)[0]
    lv = {}
    for level in C1_LEVELS:
        y = c1_coords(level)
        m = TensorMesh(y)
        reg = Regions.from_function(m, lambda *c: c[-1] < 0)
        c = np.where(reg.semi_node, -A_SUB, 0.0)
        mod = PhysicsMOSDD(reg, c, DP, LOGD, eps_ratio=EPS_R, interface=itf)
        isurf = int(np.argmin(np.abs(y)))
        bottom, top = [0], [m.n_nodes - 1]
        vg = {p: gate_voltage_it(p, A_SUB, DELTA, T_OX, EPS_R, 0.0, itf) for p in PSI_S_TARGETS}

        def bc_of(p, c=c, bottom=bottom, top=top):
            return Dirichlet().ohmic(c, DELTA, bottom, 0.0).gate(top, p[0], 0.0)

        X0 = mod.initial_guess(bc_of([vg[0.0]]))
        s0 = mod.solve_bc(X0, mod.equilibrium_bc(bc_of([vg[0.0]])))
        sols = {0.0: s0}
        up = [p for p in PSI_S_TARGETS if p > 0]
        down = sorted([p for p in PSI_S_TARGETS if p < 0], reverse=True)
        for branch in (up, down):
            out = continue_bias(
                mod, s0.X, bc_of, [vg[0.0]], [[vg[p]] for p in branch], equilibrium=True
            )
            for p, (_, s, _, _) in zip(branch, out, strict=False):
                sols[p] = s
        lv[level] = (y, isurf, vg, sols, mod)
    rows, ok, conv, worst_rel = [], True, True, 0.0
    for p in PSI_S_TARGETS:
        v = lv[C1_LEVELS[0]][2][p]
        exact = psi_s_of_vg_it(v, A_SUB, DELTA, T_OX, EPS_R, 0.0, itf)  # independent root find
        num = []
        for level in C1_LEVELS:
            _, isurf, _, sols, _ = lv[level]
            conv &= bool(sols[p].report.converged)
            num.append(float(sols[p].psi[isurf] - psi_b))
        rich = (4.0 * num[1] - num[0]) / 3.0
        err = abs(rich - exact)
        small = abs(exact) < 0.05
        okp = err <= 1e-4 if small else err <= 2e-3 * abs(exact)
        ok &= okp
        if not small:
            worst_rel = max(worst_rel, err / abs(exact))
        q_it = float(itf.charge(psi_b + exact)[0] - itf.q_f)
        rows.append(
            {"psi_s_target": p, "V_G": v, "exact": exact, "num_h": num[0], "num_h2": num[1],
             "richardson": rich, "abs_err": err, "Q_it_exact": q_it, "ok": bool(okp)}
        )  # fmt: skip
    return {
        "interface": {"q_f": itf.q_f, "dit": itf.dit, "e1": itf.e1, "e2": itf.e2},
        "mesh_nodes": {str(k): len(v[0]) for k, v in lv.items()},
        "all_converged": conv,
        "worst_rel_err_nonsmall": worst_rel,
        "worst_abs_err": max(r["abs_err"] for r in rows),
        "points": rows,
        "ok": bool(ok and conv),
    }


# ============================================================ A6
def run_a6():
    ys = np.r_[np.linspace(-3, -0.3, 4), np.linspace(-0.3, 0.0, 5)[1:], [0.25, 0.5]]
    m = TensorMesh(np.linspace(0, 1, 4), ys)
    reg = Regions.from_function(m, lambda xc, yc: yc < 0)
    rng = np.random.default_rng(20261021)
    rows, worst = [], 0.0
    for delta in (1e-24, 1e-3, 1e-12):
        nd = np.where(reg.semi_node, rng.uniform(1, 5, m.n_nodes), 0.0)
        na = np.where(reg.semi_node, rng.uniform(1, 5, m.n_nodes), 0.0)
        lg = float(np.log(delta))
        ion = Ionization(nd, na, (Level(0.6, lg + 6), Level(0.4, lg + 14)), (Level(1.0, lg + 4),))
        itf = InterfaceCharge(1.0, 0.1, -10.0, 15.0)
        for vs, beta_n, beta_p in ((0.7, 2.0, 1.0), (0.05, 1.0, 2.0), (3.0, 1.5, 1.5)):
            dd = PhysicsMOSDD(
                reg, nd - na, DDParams(delta=delta, recombination=False), lg, ionization=ion,
                mobility=Mobility(vs, 0.6 * vs, beta_n, beta_p), interface=itf,
            )  # fmt: skip
            scale = abs(lg) if delta < 1e-10 else 3.0
            for _ in range(3):
                X = np.zeros(3 * m.n_nodes)
                semi = reg.semi_node
                X[0::3] = np.where(
                    semi, rng.uniform(-scale, scale, m.n_nodes), rng.uniform(-3, 3, m.n_nodes)
                )
                X[1::3] = np.where(semi, rng.uniform(-1.5, 1.5, m.n_nodes), 0.0)
                X[2::3] = np.where(semi, rng.uniform(-1.5, 1.5, m.n_nodes), 0.0)
                bc = Dirichlet()
                r_, _ = dd._all_bc(bc)
                Ja = dd.jacobian_bc(X, r_).toarray()
                Jf = dd.fd_jacobian(X, bc, h=1e-6).toarray()
                err = np.abs(Jf - Ja).max(axis=1) / np.abs(Ja).max(axis=1)
                worst = max(worst, float(err.max()))
                rows.append(
                    {"delta": delta, "vsat": vs, "worst_row_err": float(err.max()),
                     "worst_row": [int(err.argmax() // 3), int(err.argmax() % 3)]}
                )  # fmt: skip
    return {"h": 1e-6, "n_states": len(rows), "worst_row_err": worst, "states": rows,
            "ok": bool(worst <= 1e-6)}  # fmt: skip


# ============================================================ A7 / R1
def kcl(pt):
    return abs(pt["I_source"] + pt["I_drain"] + pt["I_body"])


def transfer(dev, vl, vr, vgs, X_eq):
    out = dev.march(X_eq, [0.0, 0.0, 0.0], [[0.0, vl, vr]] + [[vg, vl, vr] for vg in vgs])
    pts = []
    for tgt, s, it, ns in out[1:]:
        pts.append(
            {"V_G": float(tgt[0]), "V_left": vl, "V_right": vr, "converged": bool(s.report.converged),
             "newton_its": it, "solves": ns, "final_row_ratio": s.report.residual_history[-1],
             **dev.currents(s)}
        )  # fmt: skip
    return bool(out[0][1].report.converged), pts


def vt_extract(vg, idr):
    vg, idr = np.asarray(vg), np.asarray(idr)
    gm = np.diff(idr) / np.diff(vg)
    k = int(np.argmax(gm))
    vmid, imid = 0.5 * (vg[k] + vg[k + 1]), 0.5 * (idr[k] + idr[k + 1])
    return float(vmid - imid / gm[k]), float(gm[k])


def run_a7_r1(level, run_output=True):
    out = {}
    pts_all = {}
    configs = {
        "none": dict(ion=True),
        "Qf": dict(ion=True, interface=InterfaceCharge(ITF_C3.q_f, 0.0)),
        "Qf+Qit": dict(ion=True, interface=ITF_C3),
        "Qf+Qit+mobility": dict(ion=True, interface=ITF_C3, mobility=MOB_C3),
    }
    devs = {}
    for name, kw in configs.items():
        t0 = time.perf_counter()
        dev = Device(level, **kw)
        eq = dev.equilibrium()
        conv0, pts = transfer(dev, 0.0, V_DS, TRANSFER_VG, eq.X)
        devs[name] = (dev, eq)
        vt, gm = vt_extract([p["V_G"] for p in pts], [p["I_drain"] for p in pts])
        pts_all[name] = pts
        out[name] = {
            "equilibrium_converged": bool(eq.report.converged),
            "first_step_converged": conv0,
            "V_T_max_gm": vt,
            "gm_max": gm,
            "I_D_at_VG_240": pts[-1]["I_drain"],
            "seconds": time.perf_counter() - t0,
            "points": pts,
        }
        print(f"  transfer[{name}]: {out[name]['seconds']:.1f}s", flush=True)
    if run_output:
        for name in ("Qf+Qit", "Qf+Qit+mobility"):
            dev, eq = devs[name]
            curve = []
            X, p_prev = None, None
            for v in OUTPUT_VDS:
                if X is None:
                    s = dev.point(eq.X, 200.0, 0.0, v)
                else:
                    s = dev.march(X, [200.0, 0.0, p_prev], [[200.0, 0.0, v]])[-1][1]
                X, p_prev = s.X.copy(), v
                cur = dev.currents(s)
                curve.append({"V_DS": v, "converged": bool(s.report.converged), **cur})
            out[name]["output_VG_200"] = curve
    # A7 on the all-models configuration
    full = pts_all["Qf+Qit+mobility"]
    gated = [(kcl(p), abs(p["I_drain"])) for p in full if abs(p["I_drain"]) >= 1e-6]
    worst = max(a / b for a, b in gated) if gated else float("nan")
    allconv = bool(
        out["Qf+Qit+mobility"]["first_step_converged"]
        and all(p["converged"] for p in full)
        and all(out[k]["equilibrium_converged"] for k in out)
    )
    a7 = {
        "level": level,
        "mesh": list(map(int, (devs["none"][0].m.nx, devs["none"][0].m.ny))),
        "models": "incomplete ionization + Caughey-Thomas + Q_f + Q_it",
        "all_converged": allconv,
        "n_points_gated": len(gated),
        "worst_rel_kcl": worst,
        "ok": bool(allconv and gated and worst <= 1e-8),
    }
    # R1 summary
    vt0 = out["none"]["V_T_max_gm"]
    r1 = {
        "level": level,
        "V_DS": V_DS,
        "interface_C3": {"q_f": ITF_C3.q_f, "dit": ITF_C3.dit, "e1": ITF_C3.e1, "e2": ITF_C3.e2},
        "mobility": {"vsat_n": MOB_C3.vsat_n, "beta_n": MOB_C3.beta_n},
        "V_T": {k: out[k]["V_T_max_gm"] for k in out},
        "dV_T_vs_none": {k: out[k]["V_T_max_gm"] - vt0 for k in out},
        "I_D_VG240": {k: out[k]["I_D_at_VG_240"] for k in out},
        "gm_max": {k: out[k]["gm_max"] for k in out},
        "I_D_VG240_ratio_mobility_over_const": out["Qf+Qit+mobility"]["I_D_at_VG_240"]
        / out["Qf+Qit"]["I_D_at_VG_240"],
        "output_VG_200": {k: out[k].get("output_VG_200") for k in ("Qf+Qit", "Qf+Qit+mobility")},
        "curves": {k: [(p["V_G"], p["I_drain"]) for p in out[k]["points"]] for k in out},
    }
    return a7, r1, out


def main() -> None:
    default_out = REPO / "results" / "stage9a"
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=default_out)
    ap.add_argument("--level", type=int, default=1, help="C3 refinement level for A7/R1")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "manifest.json").write_text(
        json.dumps(build_manifest(REPO, sys.argv, [20261021]), indent=2) + "\n"
    )
    T = {}

    def timed(name, fn, *a, **k):
        t0 = time.perf_counter()
        r = fn(*a, **k)
        T[name] = time.perf_counter() - t0
        print(f"{name}: {T[name]:.1f}s", flush=True)
        return r

    a1 = timed("A1", run_a1)
    a2 = timed("A2", run_a2)
    a3 = timed("A3", run_a3)
    a4 = timed("A4", run_a4)
    a5 = timed("A5", run_a5)
    a6 = timed("A6", run_a6)
    with np.errstate(all="ignore"):  # overflow only in rejected trial Newton iterates
        a7, r1, transfer_all = timed("A7_R1", run_a7_r1, args.level)
    res = {"A1": a1, "A2": a2, "A3": a3, "A4": a4, "A5": a5, "A6": a6, "A7": a7, "R1": r1,
           "C3_geometry": C3, "C3_transfer_full": transfer_all["Qf+Qit+mobility"]["points"]}  # fmt: skip
    crit = {k: bool(res[k]["ok"]) for k in ("A1", "A2", "A3", "A4", "A5", "A6", "A7")}
    summary = {
        "criteria": crit,
        "gate_passed": all(crit.values()),
        "values": {
            "A1_max_abs_dpsi": a1["max_abs_dpsi"],
            "A2_last_three_orders": {k: a2[k]["last_three_orders"] for k in ("1D", "2D")},
            "A3_max_rel_edge_current_diff": {
                k: a3[k]["max_rel_edge_current_diff"] for k in a3 if k != "ok"
            },
            "A4_max_rel_err": max(p["max_rel_err_total_j"] for p in a4["points"]),
            "A4_j_over_qnvsat_at_max_V": a4["j_over_qnvsat_at_max_V"],
            "A5_worst_rel_err": a5["worst_rel_err_nonsmall"],
            "A5_worst_abs_err": a5["worst_abs_err"],
            "A6_worst_row_err": a6["worst_row_err"],
            "A7_worst_rel_kcl": a7["worst_rel_kcl"],
            "R1_V_T": r1["V_T"],
            "R1_dV_T_vs_none": r1["dV_T_vs_none"],
            "R1_I_D_ratio_mobility": r1["I_D_VG240_ratio_mobility_over_const"],
        },
        "seconds": T,
        "peak_rss_mib": rss_mib(),
    }
    (args.out / "results.json").write_text(json.dumps(res, indent=2, default=float) + "\n")
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    print(json.dumps(summary, indent=2, default=float))


if __name__ == "__main__":
    main()
