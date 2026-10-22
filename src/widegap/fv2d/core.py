"""Edge-based core: conductances and generic sparse assembly."""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from widegap.fv2d.mesh import TensorMesh


def edge_conductance(mesh: TensorMesh, k_cell) -> np.ndarray:
    """G_e = sum_c k_c w_{e,c} / h_e for a cell-wise coefficient k."""
    k = np.broadcast_to(np.asarray(k_cell, dtype=float), (mesh.n_cells,))
    return (k[mesh.cells[:, 0]] * mesh.wc[:, 0] + k[mesh.cells[:, 1]] * mesh.wc[:, 1]) / mesh.h


def stiffness(mesh: TensorMesh, k_cell) -> sp.csr_matrix:
    """Unconstrained (natural zero-flux) matrix L, (L u)_i = sum_e G_e (u_i - u_j)."""
    G = edge_conductance(mesh, k_cell)
    a, b = mesh.ea, mesh.eb
    n = mesh.n_nodes
    rows = np.r_[a, b, a, b]
    cols = np.r_[a, b, b, a]
    vals = np.r_[G, G, -G, -G]
    return sp.csr_matrix((vals, (rows, cols)), shape=(n, n))


def apply_dirichlet_rows(J: sp.spmatrix, rows) -> sp.csc_matrix:
    """Replace the given rows by identity rows."""
    rows = np.asarray(rows)
    J = sp.csr_matrix(J)
    mask = np.ones(J.shape[0])
    mask[rows] = 0.0
    J = sp.diags(mask) @ J + sp.csr_matrix((np.ones(len(rows)), (rows, rows)), shape=J.shape)
    return sp.csc_matrix(J)
