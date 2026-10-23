"""Linear Poisson -div(eps grad u) = f; Dirichlet on node sets, zero flux elsewhere."""

from __future__ import annotations

import numpy as np

from widegap.fv2d.core import stiffness
from widegap.fv2d.mesh import TensorMesh
from widegap.linalg import safe_spsolve


def solve_poisson(mesh: TensorMesh, eps_cell, f_nodes, dir_nodes, dir_values):
    """Solve with symmetric elimination of Dirichlet nodes. Returns u (nodal)."""
    L = stiffness(mesh, eps_cell)
    n = mesh.n_nodes
    dn = np.asarray(dir_nodes, dtype=int)
    u = np.zeros(n)
    u[dn] = dir_values
    free = np.ones(n, dtype=bool)
    free[dn] = False
    fi = np.flatnonzero(free)
    rhs = mesh.vol * np.asarray(f_nodes, dtype=float)
    rhs_f = rhs[fi] - L[fi][:, dn] @ u[dn]
    u[fi] = safe_spsolve(L[fi][:, fi].tocsc(), rhs_f)
    return u


def boundary_fluxes(mesh: TensorMesh, eps_cell, u, f_nodes, dir_nodes):
    """Outward flux of eps grad u through the Dirichlet part of the boundary.

    Per Dirichlet node: (L u)_i - V_i f_i.  Discrete balance (P4):
    sum(fluxes) + sum_i V_i f_i = sum_i (L u)_i = 0 (zero flux elsewhere).
    Returns (fluxes at dir_nodes, balance residual, scale sum V|f|).
    """
    L = stiffness(mesh, eps_cell)
    r = L @ u - mesh.vol * np.asarray(f_nodes, dtype=float)
    dn = np.asarray(dir_nodes, dtype=int)
    flux = r[dn]
    sv = float(np.sum(mesh.vol * np.asarray(f_nodes, dtype=float)))
    balance = float(flux.sum() + sv)
    scale = float(np.sum(mesh.vol * np.abs(f_nodes)))
    return flux, balance, scale
