"""Semiconductor / oxide regions on a tensor mesh (stage 8.1).

Regions are cell-wise flags.  Derived objects:

* ``semi_node``  : nodes touching at least one semiconductor cell (carrier unknowns live here);
* ``oxide_node`` : nodes touching only oxide cells (identity rows for the quasi-Fermi levels);
* ``w_s``        : per-edge dual-face length counting only semiconductor cells, so every
                   edge on the semiconductor/oxide interface or inside the oxide carries
                   zero carrier flux automatically;
* ``vol_s``      : semiconductor part of each control volume (space charge and carrier
                   continuity volumes).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from widegap.fv2d.mesh import TensorMesh


@dataclass
class Regions:
    mesh: TensorMesh
    semi_cell: np.ndarray  # bool, per cell
    semi_node: np.ndarray  # bool, per node
    oxide_node: np.ndarray  # bool, per node (touches no semiconductor cell)
    w_s: np.ndarray  # per edge
    vol_s: np.ndarray  # per node

    @classmethod
    def from_cells(cls, mesh: TensorMesh, semi_cell) -> Regions:
        sc = np.asarray(semi_cell, dtype=bool)
        if sc.shape != (mesh.n_cells,):
            raise ValueError("semi_cell must have one flag per cell")
        w_s = (mesh.wc * sc[mesh.cells]).sum(axis=1)
        vol_s = np.zeros(mesh.n_nodes)
        touch = np.zeros(mesh.n_nodes, dtype=bool)
        for m in range(mesh.cell_nodes.shape[1]):
            nd = mesh.cell_nodes[:, m]
            np.add.at(vol_s, nd, np.where(sc, mesh.cell_corner_share, 0.0))
            touch[nd[sc]] = True
        return cls(mesh, sc, touch, ~touch, w_s, vol_s)

    @classmethod
    def from_function(cls, mesh: TensorMesh, is_semi) -> Regions:
        """``is_semi(*cell_centre_coords) -> bool array``."""
        return cls.from_cells(mesh, np.asarray(is_semi(*mesh.cell_centres.T), dtype=bool))

    def cell_values(self, semi_value, oxide_value) -> np.ndarray:
        return np.where(self.semi_cell, semi_value, oxide_value).astype(float)

    def semi_average(self, q_cell) -> np.ndarray:
        """Control-volume average of cell data over the semiconductor part only
        (zero at oxide-only nodes)."""
        q = np.where(self.semi_cell, np.asarray(q_cell, dtype=float), 0.0)
        acc = np.zeros(self.mesh.n_nodes)
        for m in range(self.mesh.cell_nodes.shape[1]):
            np.add.at(acc, self.mesh.cell_nodes[:, m], q * self.mesh.cell_corner_share)
        out = np.zeros_like(acc)
        ok = self.vol_s > 0
        out[ok] = acc[ok] / self.vol_s[ok]
        return out
