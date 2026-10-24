"""Vertical DMOSFET half-cell for the stage-10 doping optimization.

Geometry (scaled units; x in Debye lengths, y downward negative, the semiconductor
occupies y < 0 and the gate oxide 0 < y <= t_ox).  Symmetry planes x = 0 (JFET centre)
and x = x_max (body-contact centre) carry zero flux.

* JFET column     : x in [0, x_jfet],  y in [-d_body, 0], donor level ``n_jfet * exp(theta_J)``;
* p-body          : x in [x_jfet, x_max], y in [-d_body, 0], acceptor ``a_body``;
* n+ source       : x in [x_src0, x_src1], y in [-d_src, 0], donor ``b_src``;
* drift layer     : y in [-(d_body + w_drift), -d_body], all x, donor
                    ``n0_drift * exp(sum_k theta_k P_k(eta))``, eta = 2 s / w_drift - 1 with s the
                    depth below the body junction;
* n+ drain layer  : below the drift, thickness ``t_drain``, donor ``b_drain``.

Contacts: source (top nodes strictly inside the n+ source) and body (top nodes with
x >= x_bc) are shorted (V = 0); the drain is the bottom row; the gate is the top of the
oxide on x <= gate_x1.  Design variables: theta = (theta_0..theta_3, theta_J).

The net doping on a node is the semiconductor-only control-volume average of the
cell-centre doping (``Regions.semi_average``), which is linear in the cell values, so
dC/dtheta is exact and analytic (``M @ dq/dtheta``).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
from numpy.polynomial import legendre

from widegap import mos_reference as ref
from widegap.fv.drift_diffusion import DDParams
from widegap.fv2d.mesh import TensorMesh
from widegap.fv2d.mos import Dirichlet, MOSDriftDiffusion
from widegap.fv2d.regions import Regions

N_THETA = 5  # theta_0..theta_3 (drift Legendre log-profile) and theta_J (JFET level)


def refine(x, level):
    x = np.asarray(x, dtype=float)
    for _ in range(level):
        x = np.sort(np.r_[x, 0.5 * (x[:-1] + x[1:])])
    return x


def seg(lo, hi, n, q):
    """n cells on [lo, hi], cell sizes proportional to q**i (ascending coordinate)."""
    h = q ** np.arange(n)
    h = h * (hi - lo) / h.sum()
    out = lo + np.r_[0.0, np.cumsum(h)]
    out[-1] = hi
    return out


@dataclass(frozen=True)
class DMOSParams:
    # lateral layout
    x_jfet: float = 4.0
    x_src0: float = 10.0
    x_src1: float = 13.0
    x_bc: float = 17.0
    x_max: float = 18.5
    gate_x1: float = 10.0
    src_contact_gap: float = 1.0  # source contact starts this far from the channel edge
    # vertical layout
    d_src: float = 1.5
    d_body: float = 10.0
    w_drift: float = 12.0
    t_drain: float = 3.0
    t_ox: float = 0.3
    # stage-26 extensions (defaults: absent -> every earlier device is bit-identical)
    t_csl: float = (
        0.0  # current-spreading layer between the body bottom and the drift (y > -(d_body + t_csl))
    )
    n_csl: float = 0.0  # its donor level (units of N0)
    a_body_b: float = 0.0  # > 0: N_A(y) log-linear from a_body (y = 0) to a_body_b (y = -d_body)
    x_split: tuple = ()  # extra breakpoints inside the n+ source segment [x_src0, x_src1]
    nx_sub: int = 4  # cells per sub-segment of the split source segment
    ny_csl: int = 4
    q_csl: float = 1.0
    # zero-bias source-junction exclusion box for the breakdown field index (edge midpoints)
    excl_dx0: float = 3.5  # to the left of the source edge
    excl_dx1: float = 2.0  # to the right of the source edge
    excl_dy: float = 3.5  # below the source bottom
    # doping (units of N0) and physics
    a_body: float = 20.0
    b_src: float = 200.0
    b_drain: float = 200.0
    n0_drift: float = 8.0  # drift level at theta = 0
    n_jfet: float = 16.0  # JFET level at theta_J = 0
    eps_ratio: float = 0.38
    delta: float = 1.0e-24
    mu_n: float = 1.0
    mu_p: float = 1.0
    # mesh (cells per segment; level-0 mesh, refined by bisection)
    x_ch_mid: float = 5.5  # channel split point (uniform part | graded toward the source junction)
    # cells / growth ratios on [0,x_jfet] [x_jfet,x_ch_mid] [x_ch_mid,x_src0] [x_src0,x_src1]
    # [x_src1,x_bc] [x_bc,x_max]
    nx: tuple = (6, 3, 6, 5, 4, 2)
    qx: tuple = (1.0, 1.0, 0.78, 1.4, 1.0, 1.0)
    ny_drain: int = 3
    ny_drift: int = 14
    q_drift: float = 0.92  # ascending: finer toward the body junction
    ny_body: int = 11
    q_body: float = 1.0
    ny_surf: int = 8  # [-d_src, 0]
    q_surf: float = 0.6
    ny_ox: int = 2

    @property
    def y_bottom(self):
        return -(self.d_body + self.t_csl + self.w_drift + self.t_drain)


class DMOSFET:
    def __init__(self, prm: DMOSParams | None = None, level: int = 0):
        self.prm = p = prm or DMOSParams()
        self.level = level
        xb = [0.0, p.x_jfet, p.x_ch_mid, p.x_src0, p.x_src1, p.x_bc, p.x_max]
        xs = [0.0]
        for k, (lo, hi, n, q) in enumerate(zip(xb[:-1], xb[1:], p.nx, p.qx, strict=True)):
            if k == 3 and p.x_split:  # stage 26: source segment split at gate edge / contact start
                pts = [lo] + [v for v in sorted(set(p.x_split)) if lo < v < hi] + [hi]
                for a, b in zip(pts[:-1], pts[1:], strict=True):
                    xs += list(seg(a, b, p.nx_sub, 1.0)[1:])
                continue
            xs += list(seg(lo, hi, n, q)[1:])
        self.y_dtop = -(p.d_body + p.t_csl)  # top of the drift (= body bottom without a CSL)
        y_drift0 = -(p.d_body + p.t_csl + p.w_drift)
        csl = [seg(self.y_dtop, -p.d_body, p.ny_csl, p.q_csl)[1:]] if p.t_csl > 0.0 else []
        ys = np.concatenate(
            [
                seg(p.y_bottom, y_drift0, p.ny_drain, 1.0),
                seg(y_drift0, self.y_dtop, p.ny_drift, p.q_drift)[1:],
                *csl,
                seg(-p.d_body, -p.d_src, p.ny_body, p.q_body)[1:],
                seg(-p.d_src, 0.0, p.ny_surf, p.q_surf)[1:],
                np.linspace(0.0, p.t_ox, p.ny_ox + 1)[1:],
            ]
        )
        self.mesh = m = TensorMesh(refine(xs, level), refine(ys, level))
        self.reg = Regions.from_function(m, lambda xc, yc: yc < 0)
        self.dp = DDParams(delta=p.delta, mu_n=p.mu_n, mu_p=p.mu_p, recombination=False)
        self.log_delta = float(np.log(p.delta))
        self.vfb = -ref.bulk(p.a_body, p.delta)[0]
        cc = m.cell_centres
        xc, yc = cc[:, 0], cc[:, 1]
        semi = self.reg.semi_cell
        self.m_drain = semi & (yc < y_drift0)
        self.m_drift = semi & (yc >= y_drift0) & (yc < self.y_dtop)
        self.m_csl = semi & (yc >= self.y_dtop) & (yc < -p.d_body)
        upper = semi & (yc >= -p.d_body)
        self.m_jfet = upper & (xc < p.x_jfet)
        self.m_src = upper & (xc > p.x_src0) & (xc < p.x_src1) & (yc > -p.d_src)
        self.m_body = upper & ~self.m_jfet & ~self.m_src
        eta = np.where(self.m_drift, 2.0 * (self.y_dtop - yc) / p.w_drift - 1.0, 0.0)
        self.eta = eta
        self.leg = np.stack(
            [legendre.legval(eta, np.eye(4)[k]) for k in range(4)]
        )  # (4, n_cells) P_k(eta_c)
        # fixed part of the cell doping and the semiconductor control-volume averaging operator
        self.q_fixed = np.zeros(m.n_cells)
        self.q_fixed[self.m_drain] = p.b_drain
        self.q_fixed[self.m_src] = p.b_src
        self.q_fixed[self.m_body] = -self.body_acceptor(yc[self.m_body])
        self.q_fixed[self.m_csl] = p.n_csl
        ok = self.reg.vol_s > 0
        rows = np.concatenate([m.cell_nodes[:, k] for k in range(m.cell_nodes.shape[1])])
        cols = np.tile(np.arange(m.n_cells), m.cell_nodes.shape[1])
        vals = np.concatenate([np.where(semi, m.cell_corner_share, 0.0)] * m.cell_nodes.shape[1])
        acc = sp.csr_matrix((vals, (rows, cols)), shape=(m.n_nodes, m.n_cells))
        inv = np.zeros(m.n_nodes)
        inv[ok] = 1.0 / self.reg.vol_s[ok]
        self.M = sp.diags(inv) @ acc  # c_nodes = M @ q_cell
        xy = m.node_xy
        tol = 1e-9
        top = xy[:, 1] == 0.0
        self.nodes = {
            "source": np.flatnonzero(
                top & (xy[:, 0] >= p.x_src0 + p.src_contact_gap - tol) & (xy[:, 0] < p.x_src1 - tol)
            ),
            "body": np.flatnonzero(top & (xy[:, 0] >= p.x_bc - tol)),
            "drain": np.flatnonzero(xy[:, 1] == p.y_bottom),
            "gate": np.flatnonzero((xy[:, 1] == p.t_ox) & (xy[:, 0] <= p.gate_x1 + tol)),
        }
        self.nodes["source_body"] = np.union1d(self.nodes["source"], self.nodes["body"])

    # ---- doping ------------------------------------------------------------------
    def body_acceptor(self, y):
        """Body acceptor level N_A(y)/N0: ``a_body`` (uniform) or, with ``a_body_b > 0``, log-linear
        from ``a_body`` at y = 0 to ``a_body_b`` at y = -d_body (clamped to that range)."""
        p = self.prm
        if p.a_body_b <= 0.0:
            return p.a_body
        t = np.clip(-np.asarray(y, dtype=float) / p.d_body, 0.0, 1.0)
        return p.a_body * (p.a_body_b / p.a_body) ** t

    def drift_level(self, theta):
        """Drift donor level N_D(eta)/N0 at the drift cells and the P_k factors."""
        th = np.asarray(theta, dtype=float)
        return self.prm.n0_drift * np.exp(th[:4] @ self.leg)

    def cell_doping(self, theta):
        th = np.asarray(theta, dtype=float)
        q = self.q_fixed.copy()
        q[self.m_drift] = self.drift_level(th)[self.m_drift]
        q[self.m_jfet] = self.prm.n_jfet * np.exp(th[4])
        return q

    def cell_doping_jac(self, theta):
        """dq/dtheta, shape (n_cells, 5)."""
        th = np.asarray(theta, dtype=float)
        nd = self.drift_level(th)
        d = np.zeros((self.mesh.n_cells, N_THETA))
        for k in range(4):
            d[:, k] = np.where(self.m_drift, nd * self.leg[k], 0.0)
        d[self.m_jfet, 4] = self.prm.n_jfet * np.exp(th[4])
        return d

    def doping(self, theta):
        """(C at nodes, dC/dtheta (n_nodes, 5)); semiconductor-only CV average of cell values."""
        return self.M @ self.cell_doping(theta), self.M @ self.cell_doping_jac(theta)

    def drift_profile(self, theta, y):
        """N_D(y)/N0 along the drift depth (y in [-(d_body+w_drift), -d_body])."""
        p = self.prm
        eta = 2.0 * (self.y_dtop - np.asarray(y, float)) / p.w_drift - 1.0
        P = np.stack([legendre.legval(eta, np.eye(4)[k]) for k in range(4)])
        return p.n0_drift * np.exp(np.asarray(theta, float)[:4] @ P)

    def drift_profile_jac(self, theta, y):
        p = self.prm
        eta = 2.0 * (self.y_dtop - np.asarray(y, float)) / p.w_drift - 1.0
        P = np.stack([legendre.legval(eta, np.eye(4)[k]) for k in range(4)])
        nd = p.n0_drift * np.exp(np.asarray(theta, float)[:4] @ P)
        return nd[:, None] * P.T  # (len(y), 4)

    def jfet_level(self, theta):
        return self.prm.n_jfet * float(np.exp(theta[4]))

    def blocking_edge_mask(self) -> np.ndarray:
        """Semiconductor edges that enter the breakdown field index: all edges with
        semiconductor dual-face weight except those whose midpoint lies in the source-junction
        box (the n+/p source junction is at zero bias; its built-in field, ~Vbi = 115 V_T over
        a few Debye lengths, is independent of V_D and is not a blocking field)."""
        p, m = self.prm, self.mesh
        xy = m.node_xy
        mx = 0.5 * (xy[m.ea, 0] + xy[m.eb, 0])
        my = 0.5 * (xy[m.ea, 1] + xy[m.eb, 1])
        box = (
            (mx >= p.x_src0 - p.excl_dx0)
            & (mx <= p.x_src1 + p.excl_dx1)
            & (my >= -(p.d_src + p.excl_dy))
        )
        return (self.reg.w_s > 0) & ~box

    # ---- model and boundary conditions -------------------------------------------
    def model(self, theta) -> tuple[MOSDriftDiffusion, np.ndarray]:
        c, _ = self.doping(theta)
        c = np.where(self.reg.semi_node, c, 0.0)
        return MOSDriftDiffusion(
            self.reg, c, self.dp, self.log_delta, eps_ratio=self.prm.eps_ratio
        ), c

    def bc(self, c_nodes, vg, vs, vd) -> Dirichlet:
        n = self.nodes
        return (
            Dirichlet()
            .ohmic(c_nodes, self.prm.delta, n["source_body"], vs)
            .ohmic(c_nodes, self.prm.delta, n["drain"], vd)
            .gate(n["gate"], vg, self.vfb)
        )
