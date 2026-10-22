"""Tensor-product mesh (1D or 2D) with an edge list for the box method.

Node index: 1D ``k = i``; 2D ``k = i * (ny + 1) + j`` (x-major, so a thin
strip in y has a small bandwidth).  Cell index: 2D ``c = i * ny + j``.
Edge ``e = (a, b)`` has length ``h[e]`` and up to two adjacent primal cells
``cells[e, 0..1]`` (padding cell index 0 with weight 0 where absent) with
``wc[e, 0..1]`` the half-width of that cell normal to the edge (1D: weight 1,
unit cross-section).  ``w[e] = wc[e].sum()`` is the dual-face length.
"""

from __future__ import annotations

import numpy as np


class TensorMesh:
    def __init__(self, x, y=None):
        self.x = np.asarray(x, dtype=float)
        self.y = None if y is None else np.asarray(y, dtype=float)
        for v in (self.x, self.y):
            if v is not None and (v.ndim != 1 or len(v) < 2 or np.any(np.diff(v) <= 0)):
                raise ValueError("coordinates must be strictly increasing 1D arrays")
        self.dim = 1 if self.y is None else 2
        self.nx = len(self.x) - 1
        self.ny = 0 if self.y is None else len(self.y) - 1
        self._build()

    # ------------------------------------------------------------------
    def node(self, i, j=0):
        return i * (self.ny + 1) + j if self.dim == 2 else i

    def _build(self):
        x, y = self.x, self.y
        if self.dim == 1:
            n = len(x)
            self.n_nodes = n
            self.n_cells = self.nx
            self.node_xy = np.column_stack([x])
            h = np.diff(x)
            self.cell_centres = (0.5 * (x[:-1] + x[1:]))[:, None]
            self.cell_size = h  # "area" of 1D cell = length (unit cross-section)
            i = np.arange(self.nx)
            self.cell_nodes = np.column_stack([i, i + 1])
            self.cell_corner_share = 0.5 * h
            self.ea, self.eb = i, i + 1
            self.h = h
            self.cells = np.column_stack([i, i])
            self.wc = np.column_stack([np.ones(self.nx), np.zeros(self.nx)])
            self.edge_dir = np.zeros(self.nx, dtype=int)
        else:
            nx, ny = self.nx, self.ny
            hx, hy = np.diff(x), np.diff(y)
            self.n_nodes = (nx + 1) * (ny + 1)
            self.n_cells = nx * ny
            X, Y = np.meshgrid(x, y, indexing="ij")
            self.node_xy = np.column_stack([X.ravel(), Y.ravel()])
            xc, yc = 0.5 * (x[:-1] + x[1:]), 0.5 * (y[:-1] + y[1:])
            XC, YC = np.meshgrid(xc, yc, indexing="ij")
            self.cell_centres = np.column_stack([XC.ravel(), YC.ravel()])
            area = np.outer(hx, hy)
            self.cell_size = area.ravel()
            ci, cj = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij")
            ci, cj = ci.ravel(), cj.ravel()
            nn = lambda i, j: i * (ny + 1) + j  # noqa: E731
            self.cell_nodes = np.column_stack(
                [nn(ci, cj), nn(ci + 1, cj), nn(ci, cj + 1), nn(ci + 1, cj + 1)]
            )
            self.cell_corner_share = 0.25 * self.cell_size
            cid = lambda i, j: i * ny + j  # noqa: E731
            # x-edges (i,j)-(i+1,j): adjacent cells (i,j-1) [below], (i,j) [above]
            ii, jj = np.meshgrid(np.arange(nx), np.arange(ny + 1), indexing="ij")
            ii, jj = ii.ravel(), jj.ravel()
            below, above = jj > 0, jj < ny
            c0 = np.where(below, cid(ii, np.maximum(jj - 1, 0)), 0)
            c1 = np.where(above, cid(ii, np.minimum(jj, ny - 1)), 0)
            w0 = np.where(below, 0.5 * hy[np.maximum(jj - 1, 0)], 0.0)
            w1 = np.where(above, 0.5 * hy[np.minimum(jj, ny - 1)], 0.0)
            ex = (
                nn(ii, jj),
                nn(ii + 1, jj),
                hx[ii],
                np.column_stack([c0, c1]),
                np.column_stack([w0, w1]),
            )
            # y-edges (i,j)-(i,j+1): adjacent cells (i-1,j) [left], (i,j) [right]
            ii, jj = np.meshgrid(np.arange(nx + 1), np.arange(ny), indexing="ij")
            ii, jj = ii.ravel(), jj.ravel()
            left, right = ii > 0, ii < nx
            c0 = np.where(left, cid(np.maximum(ii - 1, 0), jj), 0)
            c1 = np.where(right, cid(np.minimum(ii, nx - 1), jj), 0)
            w0 = np.where(left, 0.5 * hx[np.maximum(ii - 1, 0)], 0.0)
            w1 = np.where(right, 0.5 * hx[np.minimum(ii, nx - 1)], 0.0)
            ey = (
                nn(ii, jj),
                nn(ii, jj + 1),
                hy[jj],
                np.column_stack([c0, c1]),
                np.column_stack([w0, w1]),
            )
            self.ea = np.r_[ex[0], ey[0]]
            self.eb = np.r_[ex[1], ey[1]]
            self.h = np.r_[ex[2], ey[2]]
            self.cells = np.vstack([ex[3], ey[3]])
            self.wc = np.vstack([ex[4], ey[4]])
            self.edge_dir = np.r_[np.zeros(len(ex[0]), int), np.ones(len(ey[0]), int)]
        self.n_edges = len(self.ea)
        self.w = self.wc.sum(axis=1)
        # dual-cell (control) volumes: each cell gives an equal share to each corner
        self.vol = np.zeros(self.n_nodes)
        for m in range(self.cell_nodes.shape[1]):
            np.add.at(self.vol, self.cell_nodes[:, m], self.cell_corner_share)

    # ------------------------------------------------------------------
    def cell_coeff(self, f):
        """Cell-wise array from ``f(*cell_centre_coords)`` (scalars pass through)."""
        if np.isscalar(f):
            return np.full(self.n_cells, float(f))
        return np.asarray(f(*self.cell_centres.T), dtype=float)

    def node_values(self, f):
        return np.asarray(f(*self.node_xy.T), dtype=float)

    def cv_average(self, q_cell):
        """Exact area-weighted control-volume average of cell-wise data."""
        q = np.asarray(q_cell, dtype=float)
        acc = np.zeros(self.n_nodes)
        for m in range(self.cell_nodes.shape[1]):
            np.add.at(acc, self.cell_nodes[:, m], q * self.cell_corner_share)
        return acc / self.vol

    # ----- boundaries ---------------------------------------------------
    def _side(self, axis, lo):
        c = self.node_xy[:, axis]
        v = (self.x if axis == 0 else self.y)[0 if lo else -1]
        return np.flatnonzero(c == v)

    @property
    def left(self):
        return self._side(0, True)

    @property
    def right(self):
        return self._side(0, False)

    @property
    def bottom(self):
        if self.dim == 1:
            raise ValueError("1D mesh has no bottom/top")
        return self._side(1, True)

    @property
    def top(self):
        if self.dim == 1:
            raise ValueError("1D mesh has no bottom/top")
        return self._side(1, False)

    def boundary(self):
        nodes = [self.left, self.right]
        if self.dim == 2:
            nodes += [self.bottom, self.top]
        return np.unique(np.concatenate(nodes))
