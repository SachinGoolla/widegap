"""Nonlinear Boltzmann equilibrium on the edge graph (stage-2 physics).

    -div(eps grad psi) = delta (exp(-psi) - exp(psi)) + C + s

Newton with analytic Jacobian, stage-2 step limiter (max |dpsi| <= max_step),
per-row round-off tolerance |F_i| <= 1e2 eps S_i with
S_i = sum_j |J_ij| max(|psi_j|, 1) + V_i (delta(e^-psi + e^psi) + |C| + |s|)
and a final ||dpsi||_inf <= 1e-8.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve

from widegap.fv2d.core import apply_dirichlet_rows, stiffness
from widegap.fv2d.mesh import TensorMesh

EPS_MACH = np.finfo(float).eps


@dataclass
class EqReport:
    converged: bool = False
    iterations: int = 0
    limited_steps: int = 0
    residual_history: list[float] = field(default_factory=list)  # max_i |F_i| / S_i
    last_dX: float = float("inf")


@dataclass
class EqSolution:
    psi: np.ndarray
    report: EqReport


def solve_equilibrium(
    mesh: TensorMesh,
    eps_cell,
    c_nodes,
    delta,
    dir_nodes,
    dir_values,
    s_nodes=None,
    psi0=None,
    max_iter=60,
    max_step=2.0,
) -> EqSolution:
    L = stiffness(mesh, eps_cell)
    V = mesh.vol
    c = np.asarray(c_nodes, dtype=float)
    s = np.zeros(mesh.n_nodes) if s_nodes is None else np.asarray(s_nodes, dtype=float)
    dn = np.asarray(dir_nodes, dtype=int)
    free = np.ones(mesh.n_nodes, dtype=bool)
    free[dn] = False
    psi = np.arcsinh(c / (2.0 * delta)) if psi0 is None else np.array(psi0, dtype=float)
    psi[dn] = dir_values
    rep = EqReport()
    for _ in range(max_iter + 1):
        ep, em = np.exp(psi), np.exp(-psi)
        F = L @ psi - V * (delta * (em - ep) + c + s)
        F[dn] = 0.0
        Jd = sp.diags(V * delta * (em + ep))
        J = apply_dirichlet_rows(L + Jd, dn)
        S = abs(J) @ np.maximum(np.abs(psi), 1.0) + V * (delta * (em + ep) + np.abs(c) + np.abs(s))
        ratio = np.abs(F[free]) / S[free]
        rep.residual_history.append(float(ratio.max()))
        if np.all(ratio <= 1e2 * EPS_MACH) and rep.last_dX <= 1e-8:
            rep.converged = True
            break
        if rep.iterations == max_iter:
            break
        d = spsolve(J, -F)
        step = float(np.max(np.abs(d)))
        if step > max_step:
            d *= max_step / step
            rep.limited_steps += 1
        rep.last_dX = float(np.max(np.abs(d)))
        psi = psi + d
        rep.iterations += 1
    return EqSolution(psi, rep)
