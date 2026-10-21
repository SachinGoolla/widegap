#!/usr/bin/env python
"""Independent DEVSIM (Apache-2.0, github.com/devsim/devsim) reference for the
1D abrupt pn-diode I-V used to cross-check the widegap FV drift-diffusion solver.

Physics: steady state, isothermal, Boltzmann, complete ionization, constant
mobility, NO recombination, Scharfetter-Gummel (Bernoulli) transport.

Run with the separate venv:  $V/bin/python devsim_crosscheck.py [out.json]
Does NOT import anything from widegap.
"""

import json
import math
import os
import sys

import devsim as ds
import numpy as np
from devsim.python_packages.model_create import (
    CreateContactNodeModel,
    CreateEdgeModel,
    CreateEdgeModelDerivatives,
    CreateNodeModel,
    CreateNodeModelDerivative,
    CreateSolution,
)

# ---------------- physical inputs ----------------
T = 300.0
KB = 1.380649e-23
Q = 1.602176634e-19
EPS0 = 8.8541878128e-14  # F/cm
EPS_R = 11.7
N0 = 1e16  # cm^-3
NI = 1e13  # cm^-3
MU = 1000.0  # cm^2/Vs (mu_n = mu_p)
NNODES = 1025
HALF_WIDTH_LD = 40.0
KMAX = 10
REL_ERR = 1e-12
MAX_ITER = 30
# Requested relative tolerance is 1e-12, but the Newton update norm bottoms out at
# floating-point roundoff (~1e-9 relative on the Poisson block). The solver is therefore
# run down a ladder of relative tolerances and the tightest one that converges is recorded.
REL_LADDER = [1e-12, 1e-10, 1e-8]
ABS_LADDER = [1e-30, 1e-20, 1e-17, 1e-14, 1e-10, 1e-6, 1e-3, 1e0, 3e0, 1e1, 1e2]

VT = KB * T / Q
LD = math.sqrt(EPS_R * EPS0 * VT / (Q * N0))
J_UNIT = Q * 1000.0 * VT * N0 / LD

DEV, REG = "pn", "si"
Xmax = HALF_WIDTH_LD * LD
dx = 2 * Xmax / (NNODES - 1)


SETUP_INFO = {}


def build_mesh():
    ds.create_1d_mesh(mesh="m")
    ds.add_1d_mesh_line(mesh="m", pos=-Xmax, ps=dx, tag="anode")
    ds.add_1d_mesh_line(mesh="m", pos=0.0, ps=dx, tag="mid")
    ds.add_1d_mesh_line(mesh="m", pos=Xmax, ps=dx, tag="cathode")
    ds.add_1d_contact(mesh="m", name="anode", tag="anode", material="metal")
    ds.add_1d_contact(mesh="m", name="cathode", tag="cathode", material="metal")
    ds.add_1d_region(mesh="m", material="Si", region=REG, tag1="anode", tag2="cathode")
    ds.finalize_mesh(mesh="m")
    ds.create_device(mesh="m", device=DEV)


def neutral(net):
    """Charge-neutral (n, p) for net doping ND-NA = net, with n*p = ni^2."""
    h = 0.5 * abs(net)
    big = h + math.sqrt(h * h + NI * NI)
    small = NI * NI / big
    return (big, small) if net > 0 else (small, big)


def tight_solve():
    """DC solve with the tightest (abs, rel) tolerance pair in the ladders that converges.

    DEVSIM declares convergence only when BOTH the absolute and relative error criteria
    are met, and a failed solve restores the previous solution, so laddering is safe."""
    last = None
    tried = 0
    for ab in ABS_LADDER:
        for rel in REL_LADDER:
            tried += 1
            try:
                ds.solve(
                    type="dc", absolute_error=ab, relative_error=rel, maximum_iterations=MAX_ITER
                )
                return {
                    "converged": True,
                    "relative_error_used": rel,
                    "absolute_error_used": ab,
                    "max_iterations": MAX_ITER,
                    "attempts_before_success": tried - 1,
                }
            except ds.error as e:
                last = str(e)
    return {"converged": False, "error": last}


def setup():
    build_mesh()
    for name, val in [
        ("Permittivity", EPS_R * EPS0),
        ("ElectronCharge", Q),
        ("n_i", NI),
        ("T", T),
        ("kT", KB * T),
        ("V_t", VT),
        ("mu_n", MU),
        ("mu_p", MU),
    ]:
        ds.set_parameter(device=DEV, name=name, value=val)
    ds.set_parameter(device=DEV, name="anode_bias", value=0.0)
    ds.set_parameter(device=DEV, name="cathode_bias", value=0.0)

    x = np.array(ds.get_node_model_values(device=DEV, region=REG, name="x"))
    assert len(x) == NNODES, len(x)
    # exactly-uniform check
    assert np.allclose(np.diff(x), dx, rtol=1e-9, atol=0), "mesh not uniform"
    net = np.where(x < -0.5 * dx, -N0, np.where(x > 0.5 * dx, N0, 0.0))
    ds.node_solution(device=DEV, region=REG, name="NetDoping")
    ds.set_node_values(device=DEV, region=REG, name="NetDoping", values=list(net))

    CreateSolution(DEV, REG, "Potential")
    # ---- stage 1: equilibrium Poisson with Boltzmann n,p functions of psi ----
    en = "n_i*exp(Potential/V_t)"
    hp = "n_i*exp(-Potential/V_t)"
    CreateNodeModel(DEV, REG, "IntrinsicElectrons", en)
    CreateNodeModel(DEV, REG, "IntrinsicHoles", hp)
    CreateNodeModelDerivative(DEV, REG, "IntrinsicElectrons", en, "Potential")
    CreateNodeModelDerivative(DEV, REG, "IntrinsicHoles", hp, "Potential")
    ch = "-ElectronCharge*(IntrinsicHoles - IntrinsicElectrons + NetDoping)"
    CreateNodeModel(DEV, REG, "PotentialIntrinsicCharge", ch)
    CreateNodeModelDerivative(DEV, REG, "PotentialIntrinsicCharge", ch, "Potential")
    ef = "(Potential@n0-Potential@n1)*EdgeInverseLength"
    CreateEdgeModel(DEV, REG, "ElectricField", ef)
    CreateEdgeModelDerivatives(DEV, REG, "ElectricField", ef, "Potential")
    fl = "Permittivity*ElectricField"
    CreateEdgeModel(DEV, REG, "PotentialEdgeFlux", fl)
    CreateEdgeModelDerivatives(DEV, REG, "PotentialEdgeFlux", fl, "Potential")
    ds.equation(
        device=DEV,
        region=REG,
        name="PotentialEquation",
        variable_name="Potential",
        node_model="PotentialIntrinsicCharge",
        edge_model="PotentialEdgeFlux",
        variable_update="log_damp",
    )

    # contact potential = equilibrium neutral potential + bias (numerical constants)
    contact_net = {"anode": -N0, "cathode": +N0}
    for c, nd in contact_net.items():
        n, p = neutral(nd)
        psi0 = VT * math.log(n / NI)  # == -VT*ln(p/NI)
        assert abs(psi0 + VT * math.log(p / NI)) < 1e-12 * max(1, abs(psi0))
        expr = f"Potential - {c}_bias - ({psi0!r})"
        CreateContactNodeModel(DEV, c, c + "_potmodel", expr)
        CreateContactNodeModel(DEV, c, c + "_potmodel:Potential", "1")
        ds.contact_equation(
            device=DEV, contact=c, name="PotentialEquation", node_model=c + "_potmodel"
        )
    st1 = tight_solve()
    assert st1["converged"], st1
    SETUP_INFO["equilibrium_poisson_solve"] = st1

    # ---- stage 2: full drift-diffusion (no recombination) ----
    CreateSolution(DEV, REG, "Electrons")
    CreateSolution(DEV, REG, "Holes")
    psi = np.array(ds.get_node_model_values(device=DEV, region=REG, name="Potential"))
    ds.set_node_values(device=DEV, region=REG, name="Electrons", values=list(NI * np.exp(psi / VT)))
    ds.set_node_values(device=DEV, region=REG, name="Holes", values=list(NI * np.exp(-psi / VT)))

    pne = "-ElectronCharge*(Holes - Electrons + NetDoping)"
    CreateNodeModel(DEV, REG, "PotentialNodeCharge", pne)
    CreateNodeModelDerivative(DEV, REG, "PotentialNodeCharge", pne, "Electrons", "Holes")
    ds.equation(
        device=DEV,
        region=REG,
        name="PotentialEquation",
        variable_name="Potential",
        node_model="PotentialNodeCharge",
        edge_model="PotentialEdgeFlux",
        variable_update="log_damp",
    )

    # Bernoulli / Scharfetter-Gummel edge currents (written explicitly)
    CreateEdgeModel(DEV, REG, "vdiff", "(Potential@n0 - Potential@n1)/V_t")
    CreateEdgeModel(DEV, REG, "vdiff:Potential@n0", "V_t^(-1)")
    CreateEdgeModel(DEV, REG, "vdiff:Potential@n1", "-vdiff:Potential@n0")
    CreateEdgeModel(DEV, REG, "Bern01", "B(vdiff)")
    CreateEdgeModel(DEV, REG, "Bern01:Potential@n0", "dBdx(vdiff)*vdiff:Potential@n0")
    CreateEdgeModel(DEV, REG, "Bern01:Potential@n1", "-Bern01:Potential@n0")
    # Jn = q mu_n V_t/h * (n1 B(dv) - n0 B(-dv)),  B(-x) = B(x) + x,  dv=(psi0-psi1)/Vt
    Jn = (
        "ElectronCharge*mu_n*EdgeInverseLength*V_t*"
        "(Electrons@n1*Bern01 + Electrons@n1*vdiff - Electrons@n0*Bern01)"
    )
    Jp = (
        "-ElectronCharge*mu_p*EdgeInverseLength*V_t*"
        "(Holes@n1*Bern01 - Holes@n0*Bern01 - Holes@n0*vdiff)"
    )
    CreateEdgeModel(DEV, REG, "ElectronCurrent", Jn)
    CreateEdgeModel(DEV, REG, "HoleCurrent", Jp)
    for v in ("Electrons", "Holes", "Potential"):
        CreateEdgeModelDerivatives(DEV, REG, "ElectronCurrent", Jn, v)
        CreateEdgeModelDerivatives(DEV, REG, "HoleCurrent", Jp, v)
    CreateNodeModel(DEV, REG, "NCharge", "-ElectronCharge*Electrons")
    CreateNodeModelDerivative(DEV, REG, "NCharge", "-ElectronCharge*Electrons", "Electrons")
    CreateNodeModel(DEV, REG, "PCharge", "ElectronCharge*Holes")
    CreateNodeModelDerivative(DEV, REG, "PCharge", "ElectronCharge*Holes", "Holes")
    # R = 0: no generation/recombination node_model supplied
    ds.equation(
        device=DEV,
        region=REG,
        name="ElectronContinuityEquation",
        variable_name="Electrons",
        edge_model="ElectronCurrent",
        time_node_model="NCharge",
        variable_update="positive",
    )
    ds.equation(
        device=DEV,
        region=REG,
        name="HoleContinuityEquation",
        variable_name="Holes",
        edge_model="HoleCurrent",
        time_node_model="PCharge",
        variable_update="positive",
    )

    for c, nd in contact_net.items():
        n, p = neutral(nd)
        CreateContactNodeModel(DEV, c, c + "_nmodel", f"Electrons - ({n!r})")
        CreateContactNodeModel(DEV, c, c + "_nmodel:Electrons", "1")
        CreateContactNodeModel(DEV, c, c + "_pmodel", f"Holes - ({p!r})")
        CreateContactNodeModel(DEV, c, c + "_pmodel:Holes", "1")
        ds.contact_equation(
            device=DEV,
            contact=c,
            name="ElectronContinuityEquation",
            node_model=c + "_nmodel",
            edge_current_model="ElectronCurrent",
        )
        ds.contact_equation(
            device=DEV,
            contact=c,
            name="HoleContinuityEquation",
            node_model=c + "_pmodel",
            edge_current_model="HoleCurrent",
        )


def main(out):
    setup()
    rows = []
    for k in range(KMAX + 1):
        V = k * VT
        ds.set_parameter(device=DEV, name="anode_bias", value=V)
        info = tight_solve()
        jn = ds.get_contact_current(
            device=DEV, contact="anode", equation="ElectronContinuityEquation"
        )
        jp = ds.get_contact_current(device=DEV, contact="anode", equation="HoleContinuityEquation")
        jn_c = ds.get_contact_current(
            device=DEV, contact="cathode", equation="ElectronContinuityEquation"
        )
        jp_c = ds.get_contact_current(
            device=DEV, contact="cathode", equation="HoleContinuityEquation"
        )
        Jtot_edges = np.array(
            ds.get_edge_model_values(device=DEV, region=REG, name="ElectronCurrent")
        ) + np.array(ds.get_edge_model_values(device=DEV, region=REG, name="HoleCurrent"))
        mean = Jtot_edges.mean()
        cons = float(np.max(np.abs(Jtot_edges - mean)) / abs(mean)) if abs(mean) > 0 else None
        J_anode = jn + jp
        row = {
            "k": k,
            "V": V,
            "J_anode_signed_A_cm2": J_anode,
            "J_cathode_signed_A_cm2": jn_c + jp_c,
            "abs_J_A_cm2": abs(J_anode),
            "j_scaled": abs(J_anode) / J_UNIT,
            "abs_J_edge_mean_A_cm2": abs(mean),
            "j_scaled_edge_mean": abs(mean) / J_UNIT,
            "Jn_anode": jn,
            "Jp_anode": jp,
            "convergence": info,
            "current_conservation_max_rel_dev": cons if k >= 2 else None,
        }
        rows.append(row)
        print(
            "k=%2d V=%.6f J=%.6e A/cm2 j_scaled=%.8e conv=%s cons=%s"  # noqa: UP031
            % (
                k,
                V,
                J_anode,
                row["j_scaled"],
                info["converged"],
                row["current_conservation_max_rel_dev"],
            )
        )
    ratios = {str(k): rows[k + 1]["j_scaled"] / rows[k]["j_scaled"] for k in range(4, 9)}
    ideality = {
        str(k): 1.0 / math.log(r)
        for k, r in ((k, rows[k + 1]["j_scaled"] / rows[k]["j_scaled"]) for k in range(4, 9))
    }
    result = {
        "devsim_version": ds.__version__,
        "model": "1D abrupt pn diode, steady state, isothermal, Boltzmann, complete ionization, "
        "constant mobility, R=0, Scharfetter-Gummel/Bernoulli",
        "inputs": {
            "T_K": T,
            "k_B_J_per_K": KB,
            "q_C": Q,
            "eps0_F_per_cm": EPS0,
            "eps_r": EPS_R,
            "N0_cm3": N0,
            "n_i_cm3": NI,
            "mu_n_cm2_Vs": MU,
            "mu_p_cm2_Vs": MU,
            "mu0_cm2_Vs": 1000.0,
            "domain_half_width_in_LD": HALF_WIDTH_LD,
            "n_nodes": NNODES,
            "n_elements": NNODES - 1,
            "dx_cm": dx,
            "doping_node_at_x0": 0.0,
            "contacts": "anode (left, p-side, N_A=N0), cathode (right, n-side, N_D=N0, grounded)",
            "contact_bc": "n,p = exact charge-neutral roots (n*p=n_i^2); psi = V_T*ln(n/n_i) + V_applied",
        },
        "V_T_V": VT,
        "L_D_cm": LD,
        "j_unit_A_cm2": J_UNIT,
        "solver": {
            "type": "dc (Newton, UMFPACK direct)",
            "requested_relative_error": REL_ERR,
            "relative_ladder_tried": REL_LADDER,
            "absolute_ladder_tried": ABS_LADDER,
            "maximum_iterations_per_attempt": MAX_ITER,
            "setup_solves": SETUP_INFO,
            "initial_guess": "equilibrium (Poisson-only solve), then each bias uses previous solution",
            "bias_step_in_VT": 1,
        },
        "sweep": rows,
        "ratio_J_kp1_over_J_k_k4_to_8": ratios,
        "ideality_k4_to_8": ideality,
        "e": math.e,
    }
    with open(out, "w") as f:
        json.dump(result, f, indent=2)
    print("ratios", ratios)
    print("wrote", out)


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    out = (
        sys.argv[1]
        if len(sys.argv) > 1
        else os.path.join(here, "..", "..", "results", "stage3", "devsim.json")
    )
    main(os.path.normpath(out))
