r"""Stage 12c: corrected zero-current off-state model (quasi-Fermi potential field).

Spec: docs/benchmarks/stage12c_offstate_fix.md (pre-registered; AUTHORITATIVE).

The stage-12 off-state model prescribes ``phi_n = V_D`` on n-type nodes and ``phi_n = 0`` on
p-type nodes.  A p-type node inside the drain-side depletion region has ``psi >> 0``, so
``n = delta exp(psi - phi_n)`` becomes large and spurious (stage 13b diagnosis).  Corrected model:

    phi_p = 0 (everywhere),            phi_n = V_D w(r),

with ``w = 0`` on the n-type nodes of the (enlarged) source box (the n+ source, tied to the
source), ``w = 1`` on all other n-type nodes (drift, JFET, drain), and ``w`` on the remaining
nodes the **discrete harmonic extension**: Laplace's equation on the semiconductor mesh graph
with the Poisson-type edge conductances ``gw`` (dual-face length / edge length), Dirichlet data
as above, zero flux elsewhere (``mode="harmonic"``).  The alternative interpolation used by
criterion C2 is the normalized graph distance (``mode="distance"``)

    w = d0 / (d0 + d1),

``d0``/``d1`` the shortest-path distances (physical edge lengths, Dijkstra) to the w = 0 / w = 1
sets (w = 1 if there is no w = 0 set).

Implementation choices where the spec is silent (recorded in the stage-12c result):

* Only the n-type source-box nodes and the other n-type nodes are pinned.  The p-type contact
  nodes (body contact, ...) take part in the Laplace problem as ordinary nodes; their own
  ``phi_n`` is the ohmic contact value, which overrides the field (contact rows are Dirichlet
  in the solver).  In particular a uniform p+n diode has no w = 0 set, so ``w = 1`` on the p+
  side and the model *is* the stage-12 12a one-dimensional model (``phi_n = V_D`` everywhere,
  ``phi_p = 0``).
* The classification (which node is n-type / in the source box) is that of the theta = 0 device
  and is fixed, exactly as in stages 13/13b, so ``w`` does not depend on the design variables
  and the off-state adjoint keeps its form: the Dirichlet values ``V_D w`` are affine in V_D,
  ``d(vals)/dV_D = w`` on the free semiconductor rows.

Integration without editing existing code: :class:`OffstateMixin` overrides ``offstate_bc`` of
:class:`widegap.sic.SiCDMOS` / :class:`widegap.adjoint_sic_corner.RoundedSiCDMOS` (the adjoint,
BV march, dense cuts etc. all go through ``base.offstate_bc``), and
:class:`CorrectedCornerAdjoint` builds :class:`widegap.adjoint_sic_corner.CornerAdjoint` on it.
"""

from __future__ import annotations

from unittest import mock

import numpy as np
import scipy.sparse as sp
from scipy.optimize import brentq
from scipy.sparse.csgraph import connected_components, dijkstra

import widegap.adjoint_sic_corner as _corner
from widegap.adjoint_sic_corner import SRC_MARGIN, CornerAdjoint, RoundedSiCDMOS
from widegap.fv.avalanche import ionization_integrals
from widegap.fv.drift_diffusion import DDParams
from widegap.fv2d.mesh import TensorMesh
from widegap.fv2d.mos import Dirichlet
from widegap.fv2d.physics import Ionization, Mobility
from widegap.fv2d.qf_dd import EPS_MACH
from widegap.fv2d.regions import Regions
from widegap.linalg import safe_spsolve
from widegap.sic import (
    AnisoMobilityDD,
    SiC4H,
    build_pn_diode,
    tangent_predictor,
    triangular_breakdown,
)

MODES = ("harmonic", "distance")


# ------------------------------------------------------------------------- the w field
def _edges(model, semi):
    """Semiconductor graph: (a, b, conductance gw, length h) of edges with carrier weight."""
    a, b, g = model.a, model.b, model.gw
    keep = (g > 0) & semi[a] & semi[b]
    return a[keep], b[keep], g[keep], model.mesh.h[keep]


def harmonic_w(model, semi, zero, one):
    """Discrete harmonic extension: Laplace (conductances ``gw``) on the semiconductor graph with
    w = 0 on ``zero``, w = 1 on ``one`` (index arrays); returns w (nan outside the semiconductor)."""
    n = model.mesh.n_nodes
    zero, one = np.asarray(zero, dtype=int), np.asarray(one, dtype=int)
    if len(zero) + len(one) == 0:
        raise ValueError("no Dirichlet nodes for the quasi-Fermi field")
    a, b, g, _ = _edges(model, semi)
    adj = sp.coo_matrix((g, (a, b)), shape=(n, n)).tocsr()
    adj = adj + adj.T
    lap = sp.diags(np.asarray(adj.sum(axis=1)).ravel()) - adj
    w = np.full(n, np.nan)
    w[zero], w[one] = 0.0, 1.0
    dir_mask = np.zeros(n, dtype=bool)
    dir_mask[np.r_[zero, one].astype(int)] = True
    free = np.flatnonzero(semi & ~dir_mask)
    if len(free) == 0:
        return w
    lff = lap[free][:, free].tocsc()
    ldd = lap[free][:, np.flatnonzero(dir_mask)]
    touch = np.asarray(abs(ldd).sum(axis=1)).ravel() > 0
    ncomp, lab = connected_components(adj[free][:, free], directed=False)
    for c in range(ncomp):
        if not touch[lab == c].any():
            raise ValueError("a connected component of free nodes has no Dirichlet neighbour")
    rhs = -np.asarray(lap[free][:, one].sum(axis=1)).ravel() if len(one) else np.zeros(len(free))
    w[free] = np.clip(safe_spsolve(lff, rhs), 0.0, 1.0)
    return w


def distance_w(model, semi, zero, one):
    """Normalized graph distance w = d0/(d0 + d1) (Dijkstra, physical edge lengths); w = 1 when
    there is no w = 0 set, w = 0 when there is no w = 1 set."""
    n = model.mesh.n_nodes
    zero, one = np.asarray(zero, dtype=int), np.asarray(one, dtype=int)
    if len(zero) + len(one) == 0:
        raise ValueError("no Dirichlet nodes for the quasi-Fermi field")
    a, b, _, h = _edges(model, semi)
    gr = sp.coo_matrix((h, (a, b)), shape=(n, n)).tocsr()
    inf = np.full(n, np.inf)
    d0 = dijkstra(gr, directed=False, indices=zero, min_only=True) if len(zero) else inf
    d1 = dijkstra(gr, directed=False, indices=one, min_only=True) if len(one) else inf
    w = np.full(n, np.nan)
    sel = semi
    with np.errstate(invalid="ignore", divide="ignore"):
        ww = d0 / (d0 + d1)
    ww = np.where(np.isinf(d0) & np.isfinite(d1), 1.0, ww)
    ww = np.where(np.isinf(d1) & np.isfinite(d0), 0.0, ww)
    if np.any(sel & ~np.isfinite(ww)):
        raise ValueError("a semiconductor node is unreachable from the Dirichlet sets")
    w[sel] = ww[sel]
    w[zero], w[one] = 0.0, 1.0
    return w


def w_field(model, semi, zero, one, mode="harmonic"):
    if mode == "harmonic":
        return harmonic_w(model, semi, zero, one)
    if mode == "distance":
        return distance_w(model, semi, zero, one)
    raise ValueError(f"unknown interpolation mode {mode!r}")


# ------------------------------------------------------------------------- DMOSFET mixin
class OffstateMixin:
    """Replaces ``offstate_bc`` of a :class:`SiCDMOS` (sharp) or ``RoundedSiCDMOS`` by the
    corrected one.  ``w_mode`` is a class attribute (see :func:`corrected_class`)."""

    w_mode = "harmonic"

    def _offstate_sets(self):
        dev, p = self.dev, self.prm
        xy = dev.mesh.node_xy
        semi = dev.reg.semi_node
        m = SRC_MARGIN * getattr(self, "r_s", 0.0)  # rounded device: enlarged source box
        in_src = (
            (xy[:, 0] >= p.x_src0 - m - 1e-9)
            & (xy[:, 0] <= p.x_src1 + m + 1e-9)
            & (xy[:, 1] >= -p.d_src - m - 1e-9)
        )
        n_type = (self.nd > self.na) & semi
        contact = np.zeros(len(xy), dtype=bool)
        for k in ("source_body", "drain"):
            contact[dev.nodes[k]] = True
        return semi, np.flatnonzero(n_type & in_src), np.flatnonzero(n_type & ~in_src), contact

    @property
    def offstate_w(self):
        """Fixed (theta = 0 classification) quasi-Fermi weight field, nan outside the semiconductor."""
        if "_w_off" not in self.__dict__:
            semi, zero, one, contact = self._offstate_sets()
            self._w_off = w_field(self.model, semi, zero, one, self.w_mode)
            self._w_free = np.flatnonzero(semi & ~contact)
            self._w_sets = (zero, one)
        return self._w_off

    def offstate_bc(self, vd, vg=0.0) -> Dirichlet:
        w = self.offstate_w
        bc = self.bc(vg, 0.0, vd)
        free = self._w_free
        bc.phi(free, float(vd) * w[free], 0.0)
        return bc


def corrected_class(base_cls, mode="harmonic"):
    """``base_cls`` (SiCDMOS or RoundedSiCDMOS) with the corrected off-state model."""
    if mode not in MODES:
        raise ValueError(f"unknown interpolation mode {mode!r}")
    return type(f"Corrected{base_cls.__name__}_{mode}", (OffstateMixin, base_cls), {"w_mode": mode})


class CorrectedCornerAdjoint(CornerAdjoint):
    """:class:`CornerAdjoint` (rounded device, soft-min multi-cut BV) with the stage-12c off state."""

    def __init__(self, sc, design, level=0, r_um=0.1, cuts_um=(0.0,), tau_v=5.0, w_mode="harmonic",
                 **kw):  # fmt: skip
        self.w_mode = w_mode
        with mock.patch.object(_corner, "RoundedSiCDMOS", corrected_class(RoundedSiCDMOS, w_mode)):
            super().__init__(sc, design, level, r_um=r_um, cuts_um=cuts_um, tau_v=tau_v, **kw)


def n_consistency(base, X, lo=0.01, hi=0.99, tol=1e-6):
    """C3: on the nodes with lo < w < hi, n / max(N_A, N_D) (scaled densities of the theta = 0
    device).  Gated on all such nodes (contact nodes included, where phi_n is the ohmic value);
    the worst non-contact node is reported separately.  Returns worst ratios and node data."""
    w = base.offstate_w
    sel = base.dev.reg.semi_node & (w > lo) & (w < hi)
    psi, pn = np.asarray(X)[0::3], np.asarray(X)[1::3]
    ld = float(base.model.ld)
    with np.errstate(over="ignore"):
        n = np.exp(ld + np.where(sel, psi - pn, -1e3))
    big = np.maximum(base.nd, base.na)
    ratio = np.where(sel, n / np.where(big > 0, big, 1.0), 0.0)
    free = np.zeros(len(w), dtype=bool)
    free[base._w_free] = True
    ratio_free = np.where(free, ratio, 0.0)
    sc = base.sc
    um = lambda v: float(sc.length_cm(v)) * 1e4  # noqa: E731
    xy = base.dev.mesh.node_xy

    def rec(r):
        k = int(np.argmax(r))
        return {"worst_ratio": float(r[k]), "node": k, "x_um": um(xy[k, 0]), "y_um": um(xy[k, 1]),
                "w": float(w[k]), "psi": float(psi[k]), "phi_n": float(pn[k]),
                "psi_minus_phin": float(psi[k] - pn[k]), "N": float(big[k])}  # fmt: skip

    out = rec(ratio)
    out.update({"n_nodes_checked": int(sel.sum()), "tol": tol, "ok": bool(out["worst_ratio"] <= tol),
                "free_nodes": rec(ratio_free)})  # fmt: skip
    return out


def jacobian_fd_rows(model, X, rows, vals, step=1e-2):
    """C6: analytic ``jacobian_bc`` vs central differences of ``residual_bc`` (Dirichlet rows
    ``rows``/``vals``), column by column, Richardson-combined over ``step`` and ``step/2``
    ((4 D(h/2) - D(h))/3).  Per row i: err_i = max_j |J_ij - D_ij| / max_j |J_ij|, and the
    round-off floor of the difference quotient, floor_i = 3 eps S_i / (step max_j |J_ij|) with
    S_i the residual-term magnitude of row i (``row_scale``).  Returns dict of arrays."""
    X = np.asarray(X, dtype=float)
    n = len(X)
    jac = model.jacobian_bc(X, rows)
    ja = jac.toarray()

    def fd(h):
        d = np.empty((n, n))
        for k in range(n):
            xp, xm = X.copy(), X.copy()
            xp[k] += h
            xm[k] -= h
            d[:, k] = (model.residual_bc(xp, rows, vals) - model.residual_bc(xm, rows, vals)) / (
                2 * h
            )
        return d

    d1, d2 = fd(step), fd(0.5 * step)
    dd = (4.0 * d2 - d1) / 3.0
    scale = np.max(np.abs(ja), axis=1)
    err = np.max(np.abs(ja - dd), axis=1) / scale
    err_coarse = np.max(np.abs(ja - d1), axis=1) / scale
    s_i = model.row_scale(X, jac)
    floor = 3.0 * EPS_MACH * s_i / (step * scale)
    return {"err": err, "err_single_step": err_coarse, "floor": floor, "scale": scale,
            "worst_row": int(np.argmax(err)), "n_rows": n}  # fmt: skip


# ------------------------------------------------------------- laterally uniform 2D p+n
class UniformPN2D:
    """Laterally uniform 2D p+/n diode on the stage-12 12a node set (extruded: the 1D node arrays
    ``nd``/``na`` are copied onto every column, no cell averaging).  Anode (p+) at y = 0, cathode
    (n) at the bottom; ohmic contacts on the full top/bottom rows; lateral sides are symmetry
    planes.  With ``rule="new"`` the off state is the corrected model (w from the n-type set
    only: no source, hence w = 1 everywhere), with ``rule="old"`` the stage-12 DMOSFET rule
    (phi_n = V_D on n-type nodes, 0 on p-type)."""

    def __init__(self, sc: SiC4H, nd_cm, level=0, na_cm=1e19, n_lat=2, lat_width=1.0,
                 width_factor=1.25):  # fmt: skip
        self.sc = sc
        sets = {"niwa": sc.niwa(), "hatakeyama": sc.hatakeyama()}
        nd_s = float(sc.dens_s(nd_cm))
        est = {k: triangular_breakdown(nd_s, ch) for k, ch in sets.items()}
        w_dep = max(v[1] for v in est.values())
        self.v_est = max(v[2] for v in est.values())
        diode, info = build_pn_diode(sc, nd_cm, width_factor * w_dep, na_cm, level=level,
                                     rule={"w_dep": w_dep})  # fmt: skip
        self.diode, self.info = diode, info
        x = diode.x
        y = -(x - x[0])[::-1]  # ascending, anode at y = 0
        ny = len(y) - 1
        self.mesh = mesh = TensorMesh(np.linspace(0.0, lat_width * n_lat, n_lat + 1), y)
        self.reg = Regions.from_function(mesh, lambda xc, yc: yc < 0)
        idx1d = ny - (np.arange(mesh.n_nodes) % (ny + 1))  # y ascending <-> x descending
        self.nd, self.na = diode.nd[idx1d], diode.na[idx1d]
        dp = DDParams(delta=sc.delta, mu_n=1.0, mu_p=1.0, recombination=False)
        vs = float(sc.vsat_s(sc.reg.value("vsat_n_perp_c")))
        beta = sc.reg.value("ct_beta_n_perp_c")
        self.model = AnisoMobilityDD(
            self.reg, self.nd - self.na, dp, float(np.log(sc.delta)), eps_ratio=1.0,
            ionization=Ionization(self.nd, self.na, sc.donor_levels(), sc.acceptor_levels(na_cm)),
            mobility=Mobility(vsat_n=vs, vsat_p=vs, beta_n=beta, beta_p=beta),
            interface=None, mu_edge=np.ones(mesh.n_edges),
        )  # fmt: skip
        xy = mesh.node_xy
        self.anode = np.flatnonzero(xy[:, 1] == 0.0)
        self.cathode = np.flatnonzero(xy[:, 1] == y[0])
        contact = np.zeros(mesh.n_nodes, dtype=bool)
        contact[self.anode] = contact[self.cathode] = True
        self.free = np.flatnonzero(self.reg.semi_node & ~contact)
        n_type = self.nd > self.na
        self.n_set = np.flatnonzero(n_type)
        self.col = np.flatnonzero(xy[:, 0] == 0.0)  # column x = 0, y ascending
        self.col = self.col[np.argsort(-xy[self.col, 1])]  # anode first (depth increasing)
        self.depth = -xy[self.col, 1]
        self._w = {}

    def w(self, mode="harmonic"):
        if mode not in self._w:
            self._w[mode] = w_field(self.model, self.reg.semi_node, np.array([], dtype=int),
                                    self.n_set, mode)  # fmt: skip
        return self._w[mode]

    def bc(self, vd) -> Dirichlet:
        bc = Dirichlet()
        self.model.ohmic(bc, self.anode, 0.0)
        self.model.ohmic(bc, self.cathode, float(vd))
        return bc

    def offstate_bc(self, vd, rule="new", mode="harmonic") -> Dirichlet:
        bc = self.bc(vd)
        if rule == "new":
            phi = float(vd) * self.w(mode)[self.free]
        elif rule == "old":
            phi = np.where((self.nd > self.na)[self.free], float(vd), 0.0)
        else:
            raise ValueError(rule)
        bc.phi(self.free, phi, 0.0)
        return bc

    def integral(self, X, ch):
        psi = np.asarray(X)[0::3][self.col]
        e = -np.diff(psi) / np.diff(self.depth)
        return ionization_integrals(self.depth, e, ch)[0]

    def bv(self, ch, rule="new", mode="harmonic", xtol=1e-12, v_max_factor=6.0, max_step=20.0):
        """Scaled BV: root of I_n(V_D) = 1 along the column (as stage-12 12a ``avalanche_bv``).
        March in steps of at most ``v_est / 3`` (halved on a failed Newton solve), then Brent on
        the bracket with the tangent predictor from the lower bracket end."""
        mod = self.model
        bc_of = lambda v: self.offstate_bc(v, rule, mode)  # noqa: E731
        s = mod.solve_bc(mod.initial_guess(bc_of(0.0)), bc_of(0.0), max_step=max_step)
        if not s.report.converged:
            raise RuntimeError("zero-bias solve did not converge")
        lo_v, lo_x = 0.0, s.X.copy()
        step = self.v_est / 3.0
        while True:
            v = lo_v + step
            if v > v_max_factor * self.v_est:
                raise ValueError("I = 1 not bracketed")
            pred = tangent_predictor(mod, lo_x, bc_of, lo_v, step)
            s = mod.solve_bc(pred, bc_of(v), max_step=max_step)
            if not s.report.converged:
                step *= 0.5
                if step < 1e-6 * self.v_est:
                    raise RuntimeError(f"march failed at V = {lo_v}")
                continue
            if self.integral(s.X, ch) >= 1.0:
                hi_v = v
                break
            lo_v, lo_x = v, s.X.copy()

        def g(v):
            pred = tangent_predictor(mod, lo_x, bc_of, lo_v, v - lo_v)
            sol = mod.solve_bc(pred, bc_of(v), max_step=max_step)
            if not sol.report.converged:
                raise RuntimeError(f"Brent solve failed at V = {v}")
            return self.integral(sol.X, ch) - 1.0

        return float(brentq(g, lo_v, hi_v, xtol=xtol, rtol=4 * EPS_MACH))


__all__ = [
    "MODES",
    "CorrectedCornerAdjoint",
    "OffstateMixin",
    "UniformPN2D",
    "corrected_class",
    "distance_w",
    "harmonic_w",
    "jacobian_fd_rows",
    "n_consistency",
    "w_field",
]
