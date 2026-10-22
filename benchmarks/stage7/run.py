"""Stage-7 benchmark (docs/benchmarks/stage7_sic_materials.md): M1-M6, R7a, R7b.

Reads the 4H-SiC registry src/widegap/materials/data/sic4h.json; writes
results/stage7/{results,summary,manifest}.json.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from widegap.fv.equilibrium import cv_average_piecewise
from widegap.fv.ionized_equilibrium import IonizedJunction
from widegap.manifest import build_manifest
from widegap.materials.ionization import (
    K_B,
    Level,
    Q,
    charge_antiderivative,
    ionized_donors,
    log_gamma_acceptor,
    log_gamma_donor,
    neutral_potential,
    single_donor_electrons,
)
from widegap.materials.registry import Registry
from widegap.verification import observed_orders

REPO = Path(__file__).resolve().parents[2]
REGISTRY = REPO / "src" / "widegap" / "materials" / "data" / "sic4h.json"
N0 = 1e16  # cm^-3, scaling density
REQUIRED = [
    "EGX_0",
    "passler_eps",
    "passler_theta",
    "passler_p",
    "E_exciton",
    "nu_c_band1",
    "m_c_band1",
    "m_v_band1",
    "m_v_band2",
    "Ev_band2",
    "m_v_band3",
    "Ev_band3",
    "dE_N_donor_h",
    "dE_N_donor_k",
    "N_site_fraction_h",
    "dE_Al_acceptor",
    "dE_Al_acceptor_alpha",
    "g_donor",
    "g_acceptor",
]
H_PLANCK = 6.62607015e-34  # J s, exact (SI 2019)
M0 = 9.1093837015e-31  # kg, CODATA 2018
EPS_NAMES = ["eps_r_parallel_c", "eps_r_perp_c"]


def _nc_parabolic(m_rel, t_k):
    """2 (2 pi m k T / h^2)^{3/2} in cm^-3 for one parabolic valley (spin included)."""
    return 2.0 * (2.0 * np.pi * m_rel * M0 * K_B * t_k / H_PLANCK**2) ** 1.5 * 1e-6


def derived(reg: Registry, t_k: float):
    """n_i, N_c, N_v at t_k from Persson & Lindefelt (1998) DOS masses, parabolic bands.

    N_c: nu equivalent minima of mass m_c (band 1) plus the second minimum (band 2,
    Boltzmann-weighted) if registered. N_v: three valence bands with their
    offsets, PARABOLIC approximation (the source shows the nonparabolic DOS is
    larger, so N_v is an underestimate; see stage-7 result). E_g(T) from band_gap().
    """
    kt = K_B * t_k / Q
    nc = reg.value("nu_c_band1") * _nc_parabolic(reg.value("m_c_band1"), t_k)
    if reg.usable("m_c_band2") and reg.usable("Ec_band2_offset"):
        nc += (
            reg.value("nu_c_band2")
            * _nc_parabolic(reg.value("m_c_band2"), t_k)
            * np.exp(-reg.value("Ec_band2_offset") / kt)
        )
    nv = _nc_parabolic(reg.value("m_v_band1"), t_k)
    nv += _nc_parabolic(reg.value("m_v_band2"), t_k) * np.exp(reg.value("Ev_band2") / kt)
    nv += _nc_parabolic(reg.value("m_v_band3"), t_k) * np.exp(reg.value("Ev_band3") / kt)
    eg = band_gap(reg, t_k)
    ni = np.sqrt(nc * nv) * np.exp(-eg / (2 * kt))
    return ni, nc, nv


def band_gap(reg: Registry, t_k: float) -> float:
    """E_g(T) = E_GX(T) + E_ex, simplified Passler model (Grivickas et al. 2007, Eq. 5)."""
    e0, eps, th, pp = (reg.value(n) for n in ("EGX_0", "passler_eps", "passler_theta", "passler_p"))
    egx = e0 - 0.5 * eps * th * ((1.0 + (2.0 * t_k / th) ** pp) ** (1.0 / pp) - 1.0)
    return float(egx + reg.value("E_exciton"))


def donor_levels(reg: Registry, t_k, delta, nu_c):
    """Two-site N donor (hexagonal / cubic) with the registered site fraction."""
    g = reg.value("g_donor")
    fh = reg.value("N_site_fraction_h")  # no default: must be sourced
    return (
        Level(fh, log_gamma_donor(g, reg.value("dE_N_donor_h"), t_k, delta, nu_c)),
        Level(1 - fh, log_gamma_donor(g, reg.value("dE_N_donor_k"), t_k, delta, nu_c)),
    )


def m2_m3(reg: Registry):
    out, worst = [], 0.0
    lim_ok = True
    for t in (200.0, 300.0, 500.0):
        ni, nc, nv = derived(reg, t)
        delta, nu_c = ni / N0, nc / N0
        lv = (
            Level(
                1.0,
                log_gamma_donor(reg.value("g_donor"), reg.value("dE_N_donor_h"), t, delta, nu_c),
            ),
        )
        n1 = nc / reg.value("g_donor") * np.exp(-reg.value("dE_N_donor_h") * Q / (K_B * t))
        for nd_cm in np.logspace(14, 19, 11):
            nd = nd_cm / N0
            psi = neutral_potential(nd, 0.0, delta, lv, ())
            n_num = delta * np.exp(psi) * N0
            n_ex = single_donor_electrons(nd_cm, n1)
            err = abs(n_num / n_ex - 1)
            worst = max(worst, err)
            out.append(
                {"T": t, "N_D": nd_cm, "n_numeric": n_num, "n_closed_form": n_ex, "rel_err": err}
            )
        # limits (M3)
        for nd_cm, kind in ((1e-6 * n1, "full"), (n1 / 1e-6, "freeze")):
            nd = nd_cm / N0
            psi = neutral_potential(nd, 0.0, delta, lv, ())
            frac = float(ionized_donors(psi, nd, lv) / nd)
            ok = frac > 1 - 1e-6 if kind == "full" else frac < 1e-3
            lim_ok &= ok
            out.append({"T": t, "N_D": nd_cm, "limit": kind, "ionized_fraction": frac, "ok": ok})
    return out, worst, lim_ok


def m4():
    """MMS with incomplete ionization (synthetic levels; tests the implementation)."""
    delta = 1e-6
    donors = (Level(0.5, np.log(1e-5)), Level(0.5, np.log(3e-4)))
    acceptors = (Level(1.0, np.log(2e-4)),)
    nd_f, na_f = 2.0, 0.7

    def exact(x):
        return 6.0 * np.tanh(3.0 * (x - 0.5))

    def src(x):
        t = np.tanh(3.0 * (x - 0.5))
        lap = 6.0 * 9.0 * (-2.0) * t * (1 - t**2)  # psi''
        from widegap.materials.ionization import charge

        return -lap - charge(exact(x), nd_f, na_f, delta, donors, acceptors)

    errs = []
    for n in (64, 128, 256, 512, 1024):
        x = np.linspace(0.0, 1.0, n + 1)
        j = IonizedJunction(
            x, nd_f * np.ones_like(x), na_f * np.ones_like(x), delta, donors, acceptors, src(x)
        )
        psi, rep = j.solve(np.zeros_like(x), exact(0.0), exact(1.0))
        errs.append(float(np.max(np.abs(psi - exact(x)))) if rep.converged else float("nan"))
    return {"e_inf": errs, "orders": observed_orders(errs)}


def m5(reg: Registry, h_list=(0.025, 0.0125), lo=-10.0, hi=40.0, x_bulk=(-9.0, 39.0)):
    t = 300.0
    ni, nc, nv = derived(reg, t)
    delta, nu_c, nu_v = ni / N0, nc / N0, nv / N0
    donors = donor_levels(reg, t, delta, nu_c)
    # Koizumi et al. (2009) Eq. (7): dE_A = dE_A0 - alpha * N_A^(1/3), at the junction's N_A
    de_a = reg.value("dE_Al_acceptor") - reg.value("dE_Al_acceptor_alpha") * (1e17) ** (1.0 / 3.0)
    acceptors = (Level(1.0, log_gamma_acceptor(reg.value("g_acceptor"), de_a, t, delta, nu_v)),)
    na_b, nd_b = 1e17 / N0, 1e16 / N0
    psi_p = neutral_potential(0.0, na_b, delta, donors, acceptors)
    psi_n = neutral_potential(nd_b, 0.0, delta, donors, acceptors)
    res = {"delta": delta, "psi_p": psi_p, "psi_n": psi_n, "Vbi_VT": psi_n - psi_p, "meshes": {}}
    for h in h_list:
        x = np.linspace(lo, hi, int(round((hi - lo) / h)) + 1)
        nd = cv_average_piecewise(x, [0.0], [0.0, nd_b])
        na = cv_average_piecewise(x, [0.0], [na_b, 0.0])
        j = IonizedJunction(x, nd, na, delta, donors, acceptors)
        psi0 = np.where(x < 0, psi_p, psi_n)
        psi, rep = j.solve(psi0, psi_p, psi_n)
        j0 = int(np.argmin(np.abs(x)))
        # exact first integral from the junction node to the n-side contact (true doping per element)
        a, b = psi[j0], psi[-1]
        fi = charge_antiderivative(b, nd_b, 0.0, delta, donors, acceptors) - charge_antiderivative(
            a, nd_b, 0.0, delta, donors, acceptors
        )
        e_fi = float(np.sqrt(2 * max(fi, 0.0)))
        e_el = float(abs(psi[j0 + 1] - psi[j0]) / (x[j0 + 1] - x[j0]))
        res["meshes"][str(h)] = {
            "converged": rep.converged,
            "iterations": rep.iterations,
            "limited": rep.limited_steps,
            "bulk_p_err": float(abs(psi[0] - psi_p)),
            "bulk_n_err": float(abs(psi[-1] - psi_n)),
            "psi_bulk_p": float(psi[np.argmin(np.abs(x - x_bulk[0]))]),
            "psi_bulk_n": float(psi[np.argmin(np.abs(x - x_bulk[1]))]),
            "junction_element_field": e_el,
            "first_integral_field": e_fi,
        }
    m = res["meshes"]
    e1, e2 = (
        m[str(h_list[0])]["junction_element_field"],
        m[str(h_list[1])]["junction_element_field"],
    )
    res["field_extrapolated"] = 2 * e2 - e1
    res["field_rel_err"] = (
        abs(res["field_extrapolated"] - m[str(h_list[1])]["first_integral_field"])
        / m[str(h_list[1])]["first_integral_field"]
    )
    # deep-bulk potentials vs independent neutrality roots (M5b), finest mesh
    fm = m[str(h_list[1])]
    res["bulk_vs_neutrality"] = max(abs(fm["psi_bulk_p"] - psi_p), abs(fm["psi_bulk_n"] - psi_n))
    res["Al_ionized_fraction_p_bulk"] = float(
        (delta * np.exp(-psi_p) - delta * np.exp(psi_p) + 0) / na_b
    )  # = N_A^-/N_A in neutral p bulk (p - n = N_A^-)
    return res


def m6(reg: Registry):
    names = ["mu_n_max_par_c", "mu_n_min_par_c", "N_n_ref_par_c", "alpha_n_par_c"]
    if not all(reg.usable(n) for n in names):
        return {"available": False}
    mx, mn, nref, al = (reg.value(n) for n in names)

    def mu(n):
        return mn + (mx - mn) / (1 + (n / nref) ** al)

    return {
        "available": True,
        "mu_low_N_err": abs(mu(1e-30) - mx) / mx,
        "mu_Nref_err": abs(mu(nref) - 0.5 * (mx + mn)) / mx,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=REPO / "results" / "stage7")
    ap.add_argument("--criteria", choices=["v1", "7b"], default="v1")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "manifest.json").write_text(
        json.dumps(build_manifest(REPO, sys.argv, []), indent=2) + "\n"
    )
    reg = Registry.from_json(REGISTRY)
    missing = reg.missing(REQUIRED)
    eps_ok = any(reg.usable(n) for n in EPS_NAMES)
    res = {
        "registry_names": reg.names(),
        "missing_required": missing,
        "permittivity_usable": eps_ok,
    }
    crit = {"M1": not missing and eps_ok}
    # M4 (synthetic levels) and M6 (mobility identities) do not depend on set R
    res["M4"] = m4()
    crit["M4"] = all(1.9 <= q <= 2.1 for q in res["M4"]["orders"][-3:])
    res["M6"] = m6(reg)
    if res["M6"]["available"]:
        crit["M6"] = res["M6"]["mu_low_N_err"] <= 1e-14 and res["M6"]["mu_Nref_err"] <= 1e-14
    if crit["M1"]:
        m23, worst, lim_ok = m2_m3(reg)
        res["M2_M3"] = m23
        res["M2_worst_rel_err"] = worst
        crit["M2"] = worst <= 1e-10
        crit["M3"] = lim_ok
        res["M5"] = (
            m5(reg) if args.criteria == "v1" else m5(reg, lo=-60.0, hi=80.0, x_bulk=(-30.0, 45.0))
        )
        mm = res["M5"]["meshes"].values()
        crit["M5"] = (
            all(v["converged"] for v in mm)
            and res["M5"]["bulk_vs_neutrality"] <= 1e-10
            and res["M5"]["field_rel_err"] <= 2e-3
        )
        ni3, nc3, nv3 = derived(reg, 300.0)
        res["R7a_ni_300K_cm3"] = ni3
        res["R7a_Nc_Nv_300K_cm3"] = [nc3, nv3]
        res["R7a_Eg_300K_eV"] = band_gap(reg, 300.0)
        res["R7a_Nc_vs_Grivickas_quoted_1.66e19"] = nc3 / 1.66e19
    else:
        res["blocked"] = (
            "M2, M3, M5 not evaluable: required parameters lack a readable, verified primary source"
        )
    crit = {k: bool(v) for k, v in crit.items()}
    summary = {
        "criteria": crit,
        "gate_passed": all(crit.values()) and "M5" in crit,
        "missing_required": missing,
    }
    (args.out / "results.json").write_text(json.dumps(res, indent=2, default=float) + "\n")
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
