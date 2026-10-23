r"""Stage 16: sourced inversion-layer (channel) mobility of the 4H-SiC DMOSFET (forward model).

Spec: docs/benchmarks/stage16_channel_mobility.md (pre-registered; AUTHORITATIVE).

Noguchi et al. (IEDM 2017 = N17, JJAP 2020 = N20), Si-face 4H-SiC, NO-nitrided oxide, 300 K:

    1/mu_ch = 1/mu_ph + 1/mu_SR + 1/mu_C                                   (N20 Eq. 4)
    mu_ph   = 66.5 E_eff^-0.39,   mu_SR = 76.7 E_eff^-2                    (N17 Eqs. 1, 2)
    mu_C    = mu_C0 (N_S/N_S0)^gamma,  mu_C0 = A (N_dep/N_dep0)^B,
              gamma = C log10(N_dep/N_dep0) + D                            (N20 Eqs. 6-8)
    E_eff   = (q/eps_SiC)(N_dep + eta N_S)   [MV/cm],   eta = 1/3           (N20 Eq. 2)
    N_dep   = sqrt((2 eps_SiC N_A/q)[(2kT/q) ln(N_A/n_i) - V_B]),  V_B = 0  (N20 Eq. 3)

mu in cm^2/(V s), N_S and N_dep in cm^-2.  Every Noguchi number is read from the registry
(``noguchi_*``); the lattice quantities (eps_SiC, n_i, kT/q) are the registry-derived ones of
:class:`widegap.sic.SiC4H`.

Discretization choices recorded in section 2 of the spec (and the ones the spec leaves open,
recorded here BEFORE the first dev run and not revisited after seeing H5/H6):

* ``N_S(x_i)``: trapezoid on the mesh NODES of the column with depth <= d_ch (no interpolation
  to exactly d_ch; this equals the box-method control-volume electron charge), mobile electrons
  only.  A column is a *channel column* when its surface node is p-type (N_A > N_D).
* The mobility of an edge between two channel columns is the HARMONIC mean of the two column
  values (series conduction); lateral (x) electron edges with both nodes within d_ch of the
  interface, both columns p-type.  Everything else keeps the bulk mobility.
* Clamps: ``N_S <- max(N_S, 1e11)``; E_eff is computed from the clamped N_S; then
  ``E_eff <- max(E_eff, 1e-2)`` MV/cm.  Each clamp is counted separately, for all channel columns
  and for the columns under the gate.
* Picard: plain lagged fixed point, no under-relaxation.  Each iteration is one warm-started
  Newton solve of the drift-diffusion system with FIXED per-edge mobility (the analytic Jacobian
  of :class:`widegap.sic.AnisoMobilityDD` is exact for fixed mu), then the mobility update;
  converged when ``max |dmu_ch/mu_ch| < 1e-10`` over the channel columns.
* "Channel share" of R_on: the quasi-Fermi-potential drop of the surface electrons across the
  gated p-type surface span divided by V_DS (a property of one solution, defined identically for
  both mobility models), and the model delta ``1 - R_bulk/R_channel``.

The hook into the solver is the existing ``AnisoMobilityDD.mu_edge`` attribute (a per-edge
low-field mobility factor): this module overwrites it between Newton solves, so ``sic.py`` is not
modified and the default path is bit-identical.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from widegap.fv.drift_diffusion import DDParams
from widegap.fv2d.mesh import TensorMesh
from widegap.fv2d.mos import Dirichlet
from widegap.fv2d.physics import InterfaceCharge, Ionization, Mobility
from widegap.fv2d.regions import Regions
from widegap.materials.ionization import Q
from widegap.sic import EPS0, AnisoMobilityDD, SiC4H

NOGUCHI_NAMES = (
    "noguchi_mu_ph_coeff",
    "noguchi_mu_ph_exp",
    "noguchi_mu_sr_coeff",
    "noguchi_mu_sr_exp",
    "noguchi_coul_A",
    "noguchi_coul_B",
    "noguchi_coul_C",
    "noguchi_coul_D",
    "noguchi_coul_NS0",
    "noguchi_coul_Ndep0",
    "noguchi_eta",
)
# validity clamps (model choices, spec section 2; lower ends of the plotted fit ranges)
ACC_NAME = "mu_acc_n_zhao2007"  # stage 18: accumulation-layer mobility, pessimistic bracket end
NS_MIN_CM2 = 1.0e11
EEFF_MIN_MV_CM = 1.0e-2
NS_EXTRAP_CM2 = 1.0e13  # columns above this are an extrapolation (counted, not clamped)
NA_FIT_RANGE_CM3 = (1.0e16, 4.0e17)  # measured N_A range of N20 (extrapolation flag)
PICARD_TOL = 1.0e-10
PICARD_MAX = 60


@dataclass(frozen=True)
class NoguchiParams:
    """The sourced Noguchi constants (read through a registry, so provenance is enforced)."""

    ph_coeff: float
    ph_exp: float
    sr_coeff: float
    sr_exp: float
    A: float
    B: float
    C: float
    D: float
    ns0: float
    ndep0: float
    eta: float

    @classmethod
    def from_registry(cls, reg) -> NoguchiParams:
        v = reg.value
        return cls(
            v("noguchi_mu_ph_coeff"), v("noguchi_mu_ph_exp"),
            v("noguchi_mu_sr_coeff"), v("noguchi_mu_sr_exp"),
            v("noguchi_coul_A"), v("noguchi_coul_B"), v("noguchi_coul_C"), v("noguchi_coul_D"),
            v("noguchi_coul_NS0"), v("noguchi_coul_Ndep0"), v("noguchi_eta"),
        )  # fmt: skip


def acc_mobility_from_registry(reg) -> float:
    """Zhao (2007) accumulation-layer mobility [cm^2/Vs] through a registry (provenance enforced)."""
    return float(reg.value(ACC_NAME))


def eps_sic_f_per_cm(sc: SiC4H) -> float:
    return EPS0 * sc.eps_r * 1e-2


def depletion_sheet(sc: SiC4H, na_cm3, v_b=0.0):
    """N_dep [cm^-2] (N20 Eq. 3)."""
    na = np.asarray(na_cm3, dtype=float)
    phi = 2.0 * sc.vt * np.log(na / sc.n_i) - v_b
    return np.sqrt(2.0 * eps_sic_f_per_cm(sc) * na / Q * phi)


def noguchi_mobility(sc: SiC4H, p: NoguchiParams, na_cm3, ns_cm2, clamp=True):
    """Total channel mobility and its parts for arrays of (N_A [cm^-3], N_S [cm^-2]).

    Returns a dict of arrays: mu, mu_ph, mu_sr, mu_c, e_eff [MV/cm], n_dep, gamma, ns_used and the
    boolean flags ``clamp_ns``, ``clamp_eeff``, ``ns_extrap``, ``gamma_neg``, ``na_out_of_range``.
    """
    na = np.asarray(na_cm3, dtype=float)
    ns_in = np.asarray(ns_cm2, dtype=float)
    ns = np.maximum(ns_in, NS_MIN_CM2) if clamp else ns_in
    ndep = depletion_sheet(sc, na)
    e_raw = (Q / eps_sic_f_per_cm(sc)) * (ndep + p.eta * ns) * 1e-6  # MV/cm
    e = np.maximum(e_raw, EEFF_MIN_MV_CM) if clamp else e_raw
    mu_ph = p.ph_coeff * e**p.ph_exp
    mu_sr = p.sr_coeff * e**p.sr_exp
    r = ndep / p.ndep0
    gamma = p.C * np.log10(r) + p.D
    mu_c = p.A * r**p.B * (ns / p.ns0) ** gamma
    mu = 1.0 / (1.0 / mu_ph + 1.0 / mu_sr + 1.0 / mu_c)
    return {
        "mu": mu, "mu_ph": mu_ph, "mu_sr": mu_sr, "mu_c": mu_c, "e_eff": e, "n_dep": ndep,
        "gamma": gamma, "ns_used": ns,
        "clamp_ns": ns_in < NS_MIN_CM2, "clamp_eeff": e_raw < EEFF_MIN_MV_CM,
        "ns_extrap": ns_in > NS_EXTRAP_CM2, "gamma_neg": gamma < 0.0,
        "na_out_of_range": (na < NA_FIT_RANGE_CM3[0]) | (na > NA_FIT_RANGE_CM3[1]),
    }  # fmt: skip


# ------------------------------------------------------------------------- geometry / selection
class ChannelMap:
    """Channel columns, their node sets and the selected lateral electron edges of one device.

    ``model`` is an :class:`AnisoMobilityDD`; ``nd``/``na`` are the node doping arrays (scaled, the
    ones the model was built from).  ``gate_x_max`` (scaled) limits the *gated* subset used for
    reporting; the mobility replacement itself applies to all p-type surface columns (spec).

    Stage 18 (docs/benchmarks/stage18_accumulation_bracket.md): ``mu_acc_cm2`` (default ``None`` =
    stage-16 behaviour, bit-identical) sets a CONSTANT mobility [cm^2/Vs] on the lateral electron
    edges (both nodes within d_ch) of the JFET accumulation surface: the contiguous run of n-type
    surface columns starting at the first column (x = 0) up to the first p-type surface column,
    restricted to x <= ``gate_x_max`` (under the gate; required with ``mu_acc_cm2``).  This
    excludes the n+ source side of the channel.

    Stage 18b (docs/benchmarks/stage18b_accumulation_excess.md): ``acc_mode="all"`` (default) is the
    stage-18 behaviour (constant ``mu_acc`` on every accumulation edge, bit-identical).
    ``acc_mode="excess"`` applies ``mu_acc`` to the gate-induced charge only:

        mu_edge = mu_bulk (1 - f) + mu_acc f,   f = max(nbar - n_fb, 0) / nbar,

    ``nbar`` = mean of the two node electron densities of the edge, ``n_fb`` = mean of the two node
    flat-band densities (``f = 0`` and mu_edge = mu_bulk in depletion).  The node flat-band density
    is ``N_D^+ - N_A^-`` of the model's own incomplete-ionization levels (``model.ion``) at the
    node's neutral-bulk potential psi_0 (root of rho(psi_0) = 0 for that node's scaled (N_D, N_A),
    ``widegap.materials.ionization.neutral_potential``), computed once here from the doping (floor
    1e-30).  ``mu_acc_cm2="bulk"`` (test hook for criterion B1) uses each edge's own bulk mobility.
    The edge mobility then depends on n: :func:`picard_solve` lags it and includes its relative
    change in the convergence measure.
    """

    def __init__(self, sc: SiC4H, model, nd, na, d_ch_um=20.0, gate_x_max=None, x_gate_min=None,
                 mu_acc_cm2=None, acc_mode="all"):  # fmt: skip
        self.sc, self.model = sc, model
        mesh = model.mesh
        self.mesh = mesh
        self.d_ch_um = float(d_ch_um)
        self.d_ch_s = float(sc.length_s(d_ch_um * 1e-4))
        xy = mesh.node_xy
        ny1 = mesh.ny + 1
        semi = model.reg.semi_node
        ix = np.arange(mesh.n_nodes) // ny1
        self.node_col = ix
        tol = 1e-9 * self.d_ch_s
        self.in_ch = semi & (xy[:, 1] <= 0.0) & (xy[:, 1] >= -self.d_ch_s - tol)
        surf = semi & (xy[:, 1] == 0.0)
        self.surf_node = np.full(mesh.nx + 1, -1, dtype=int)
        self.surf_node[ix[surf]] = np.flatnonzero(surf)
        has = self.surf_node >= 0
        s = np.where(has, self.surf_node, 0)
        self.col_p = has & (np.asarray(na)[s] > np.asarray(nd)[s])
        self.cols = np.flatnonzero(self.col_p)
        self.col_x = mesh.x
        self.col_na_cm3 = np.asarray(na)[s] * sc.n0
        # gated subset
        gm = np.ones(mesh.nx + 1, dtype=bool)
        if gate_x_max is not None:
            gm &= mesh.x <= gate_x_max + 1e-9
        if x_gate_min is not None:
            gm &= mesh.x >= x_gate_min - 1e-9
        self.col_gated = self.col_p & gm
        # per-column node lists (ordered from the surface down)
        self.col_nodes = {}
        for i in self.cols:
            nodes = np.flatnonzero(self.in_ch & (ix == i))
            self.col_nodes[int(i)] = nodes[np.argsort(-xy[nodes, 1])]
        # lateral electron edges: x-direction, both nodes in the layer, both columns p-type
        a, b = model.a, model.b
        sel = (
            (mesh.edge_dir == 0)
            & self.in_ch[a]
            & self.in_ch[b]
            & self.col_p[ix[a]]
            & self.col_p[ix[b]]
            & (model.gw > 0.0)
        )
        self.edges = np.flatnonzero(sel)
        self.edge_ca, self.edge_cb = ix[a[self.edges]], ix[b[self.edges]]
        self.mu_bulk = np.asarray(model.mu_edge, dtype=float).copy()  # snapshot (low-field, scaled)
        # ---- stage 18: JFET accumulation surface (n-type columns under the gate)
        if acc_mode not in ("all", "excess"):
            raise ValueError("acc_mode must be 'all' or 'excess'")
        self.acc_mode = acc_mode
        self.mu_acc_bulk = isinstance(mu_acc_cm2, str)
        if self.mu_acc_bulk and (mu_acc_cm2 != "bulk" or acc_mode != "excess"):
            raise ValueError("mu_acc_cm2='bulk' is only defined for acc_mode='excess'")
        self.mu_acc_cm2 = (
            None if mu_acc_cm2 is None else (mu_acc_cm2 if self.mu_acc_bulk else float(mu_acc_cm2))
        )
        self.nfb_node = None
        self.col_acc = np.zeros(mesh.nx + 1, dtype=bool)
        self.acc_edges = np.empty(0, dtype=int)
        if self.mu_acc_cm2 is not None:
            if gate_x_max is None:
                raise ValueError("mu_acc_cm2 requires gate_x_max (columns under the gate)")
            n_type = has & (np.asarray(nd)[s] > np.asarray(na)[s])
            first_not_n = np.flatnonzero(~n_type)
            run = np.arange(mesh.nx + 1) < (first_not_n[0] if len(first_not_n) else mesh.nx + 1)
            self.col_acc = n_type & run & (mesh.x <= gate_x_max + 1e-9)
            acc = (
                (mesh.edge_dir == 0)
                & self.in_ch[a]
                & self.in_ch[b]
                & self.col_acc[ix[a]]
                & self.col_acc[ix[b]]
                & (model.gw > 0.0)
            )
            self.acc_edges = np.flatnonzero(acc)
            self._ix = ix
            if acc_mode == "excess":
                self._init_excess(nd, na, ix)

    def _init_excess(self, nd, na, ix):
        """Node flat-band densities (scaled) and psi_0 on the in-layer nodes of accumulation
        columns; per-edge n_fb (mean of its two nodes) and the scaled mu_acc."""
        from widegap.materials.ionization import (
            ionized_acceptors,
            ionized_donors,
            neutral_potential,
        )

        model = self.model
        if model.ion is None:
            raise ValueError("acc_mode='excess' needs the model's incomplete-ionization levels")
        delta = float(np.exp(model.ld))
        nodes = np.flatnonzero(self.in_ch & self.col_acc[ix])
        nd, na = np.asarray(nd, float), np.asarray(na, float)
        pairs = np.unique(np.c_[nd[nodes], na[nodes]], axis=0)
        psi0 = {
            (a_, b_): neutral_potential(a_, b_, delta, model.ion.donors, model.ion.acceptors)
            for a_, b_ in pairs
        }
        self.psi0_node = np.full(self.mesh.n_nodes, np.nan)
        self.nfb_node = np.full(self.mesh.n_nodes, np.nan)
        for i in nodes:
            ps = psi0[(nd[i], na[i])]
            nfb = ionized_donors(ps, nd[i], model.ion.donors) - ionized_acceptors(
                ps, na[i], model.ion.acceptors
            )
            self.psi0_node[i], self.nfb_node[i] = ps, max(float(nfb), 1e-30)
        e = self.acc_edges
        self.nfb_edge = 0.5 * (self.nfb_node[model.a[e]] + self.nfb_node[model.b[e]])
        self.mu_acc_edge_s = (
            self.mu_bulk[e] if self.mu_acc_bulk else np.full(len(e), self.sc.mu_s(self.mu_acc_cm2))
        )

    def refresh_doping(self, nd, na):
        """Recompute the excess-mode flat-band densities after the doping changed (theta moves);
        the column / edge sets stay frozen.  Call after ``mu_bulk`` was updated."""
        if self.acc_mode == "excess" and self.mu_acc_cm2 is not None:
            self._init_excess(nd, na, self._ix)

    def excess_fraction(self, n_scaled):
        """f = max(nbar - n_fb, 0)/nbar on the accumulation edges (excess mode)."""
        e, m = self.acc_edges, self.model
        nbar = 0.5 * (np.asarray(n_scaled)[m.a[e]] + np.asarray(n_scaled)[m.b[e]])
        return np.maximum(nbar - self.nfb_edge, 0.0) / nbar

    def sheet_density(self, n_scaled):
        """N_S [cm^-2] for every p-type surface column (nan elsewhere): trapezoid on nodes."""
        sc, xy = self.sc, self.mesh.node_xy
        out = np.full(self.mesh.nx + 1, np.nan)
        for i, nodes in self.col_nodes.items():
            depth = -xy[nodes, 1] * sc.l_d  # cm
            dens = np.asarray(n_scaled)[nodes] * sc.n0  # cm^-3
            out[i] = float(np.sum(0.5 * (dens[1:] + dens[:-1]) * np.diff(depth)))
        return out

    def column_mobility(self, p: NoguchiParams, n_scaled, clamp=True):
        """(per-column result dict restricted to channel columns, ns array over all columns)."""
        ns_all = self.sheet_density(n_scaled)
        c = self.cols
        res = noguchi_mobility(self.sc, p, self.col_na_cm3[c], ns_all[c], clamp=clamp)
        return res, ns_all

    def edge_mobility(self, mu_col_cm2, n_scaled=None):
        """New scaled per-edge mu array: bulk everywhere except the selected channel edges.
        ``n_scaled`` (node electron densities) is needed only for ``acc_mode='excess'``."""
        mu = self.mu_bulk.copy()
        full = np.full(self.mesh.nx + 1, np.nan)
        full[self.cols] = mu_col_cm2
        ma, mb = full[self.edge_ca], full[self.edge_cb]
        harm = 2.0 * ma * mb / (ma + mb)
        mu[self.edges] = self.sc.mu_s(harm)
        if self.mu_acc_cm2 is not None and self.acc_mode == "all":  # stage 18
            mu[self.acc_edges] = self.sc.mu_s(self.mu_acc_cm2)
        elif self.mu_acc_cm2 is not None:  # stage 18b: accumulation mobility on the excess charge
            if n_scaled is None:
                raise ValueError("acc_mode='excess' needs n_scaled")
            f = self.excess_fraction(n_scaled)
            mu[self.acc_edges] = self.mu_bulk[self.acc_edges] * (1.0 - f) + self.mu_acc_edge_s * f
        return mu

    def counts(self, res):
        gated = self.col_gated[self.cols]
        out = {"n_channel_columns": int(len(self.cols)), "n_gated_columns": int(gated.sum()),
               "n_channel_edges": int(len(self.edges))}  # fmt: skip
        for k in ("clamp_ns", "clamp_eeff", "ns_extrap", "gamma_neg", "na_out_of_range"):
            out[k] = {"all": int(res[k].sum()), "gated": int((res[k] & gated).sum())}
        return out

    def drop_share(self, pn_scaled, vds_scaled):
        """Quasi-Fermi drop of the surface electrons across the gated p-type span / V_DS."""
        g = np.flatnonzero(self.col_gated)
        if len(g) < 2 or vds_scaled == 0.0:
            return float("nan")
        n0, n1 = self.surf_node[g[0]], self.surf_node[g[-1]]
        return float(abs(pn_scaled[n0] - pn_scaled[n1]) / abs(vds_scaled))


# ------------------------------------------------------------------------------- Picard driver
def picard_solve(cmap: ChannelMap, p: NoguchiParams, X0, bc, tol=PICARD_TOL, max_iter=PICARD_MAX,
                 mode="channel", max_step=20.0, mu_init=None):  # fmt: skip
    """Lagged-mobility solve: DD solve with fixed mu_edge -> update mu_ch -> repeat.

    ``mode='bulk'``: mu_ch is replaced by the bulk mobility on the channel edges (code-path limit,
    criterion H2); ``mode='channel'``: Noguchi.  The starting mobility is evaluated from the
    densities of ``X0`` (or ``mu_init`` if given).  Leaves ``model.mu_edge`` at the last used
    mobility.  Returns a dict (never raises on non-convergence; every flag is recorded).
    """
    model = cmap.model
    n_start = model.densities(np.asarray(X0, float))[0]

    def update(n_scaled):
        res, ns_all = cmap.column_mobility(p, n_scaled)
        if mode == "bulk":
            return cmap.mu_bulk.copy(), res, ns_all, None
        return cmap.edge_mobility(res["mu"], n_scaled), res, ns_all, res["mu"]

    if mu_init is not None:
        mu_edge, mu_cols = mu_init
    else:
        mu_edge, _, _, mu_cols = update(n_start)
    model.mu_edge = np.asarray(mu_edge, dtype=float)
    X = np.asarray(X0, dtype=float).copy()
    hist, newton_its, newton_ok = [], [], []
    sol, converged, rel = None, False, np.inf
    for _ in range(max_iter):
        sol = model.solve_bc(X, bc, max_step=max_step)
        newton_its.append(int(sol.report.iterations))
        newton_ok.append(bool(sol.report.converged))
        if not sol.report.converged:
            break
        X = sol.X.copy()
        mu_edge_new, res, ns_all, mu_cols_new = update(sol.n)
        if mode == "bulk":
            rel = 0.0
        else:
            rel = float(np.max(np.abs(mu_cols_new / mu_cols - 1.0)))
            if cmap.acc_mode == "excess" and len(cmap.acc_edges):  # n-dependent accumulation mu
                e = cmap.acc_edges
                rel = max(rel, float(np.max(np.abs(mu_edge_new[e] / model.mu_edge[e] - 1.0))))
        hist.append(rel)
        model.mu_edge = np.asarray(mu_edge_new, dtype=float)
        mu_cols = mu_cols_new
        if rel < tol:
            converged = True
            break
    out = {
        "converged": bool(converged and all(newton_ok)), "picard_converged": bool(converged),
        "newton_all_converged": bool(all(newton_ok)), "iterations": len(newton_its),
        "final_rel": rel, "rel_history": hist, "newton_its": newton_its, "sol": sol,
        "mu_cols_cm2": None if mu_cols is None else np.asarray(mu_cols),
        "mode": mode,
    }  # fmt: skip
    if sol is not None and sol.report.converged:
        res, ns_all = cmap.column_mobility(p, sol.n)
        out["columns"] = res
        out["ns_all"] = ns_all
        out["counts"] = cmap.counts(res)
        out["mu_edge"] = np.asarray(model.mu_edge).copy()
    return out


# ------------------------------------------------------------ long-channel lateral n-MOSFET (H3)
def _graded(a, b, h0, hmax, g=1.25):
    """Points from a to b (a < b) with spacing h0 at a growing by g to at most hmax."""
    pts, h = [a], h0
    while pts[-1] + h < b - 0.5 * h:
        pts.append(pts[-1] + h)
        h = min(h * g, hmax)
    pts.append(b)
    return np.array(pts)


class LateralNMOS:
    """Uniform-body long-channel n-MOSFET (spec H3), built on the stage-9a physics machinery.

    Layout (um): n+ source [0, ext], channel [ext, ext + L] (metallurgical junctions at the n+
    edges, so L is the distance between them), n+ drain [ext + L, 2 ext + L]; gate from
    ``ext - overlap`` to ``ext + L + overlap`` (overlaps both n+ regions: no ungated gap).
    The p-body has uniform ``N_A`` and depth ``depth_um``; the body contact is the bottom row.
    Same physics as the DMOSFET: registry D_it, ``sc.vfb``, Caughey-Thomas saturation.
    """

    def __init__(self, sc: SiC4H, na_cm3, L_um=10.0, ext_um=1.0, overlap_um=0.5, d_nplus_um=0.3,
                 depth_um=1.0, t_ox_um=0.05, n_plus_cm3=1e19, h_junction_um=0.02,
                 h_surface_nm=0.1, qf_cm2=0.0, dit_scale=1.0, dit_width_ev=0.2):  # fmt: skip
        self.sc, self.na_cm3, self.L_um = sc, float(na_cm3), float(L_um)
        um = 1e-4
        ls = lambda v: float(sc.length_s(v * um))  # noqa: E731
        x_s0, x_d0 = ext_um, ext_um + L_um
        x_end = 2 * ext_um + L_um
        mid = 0.5 * (x_s0 + x_d0)
        gate0, gate1 = x_s0 - overlap_um, x_d0 + overlap_um
        hmax = 0.4
        left = _graded(x_s0, mid, h_junction_um, hmax)
        right = (x_s0 + x_d0) - left  # mirror about the channel midpoint
        src = x_s0 - _graded(0.0, x_s0, h_junction_um, 0.2)  # fine at the junction, coarse outward
        drn = x_d0 + _graded(0.0, x_end - x_d0, h_junction_um, 0.2)
        xs = np.unique(np.round(np.r_[src, left, right, drn, gate0, gate1], 12))
        z = [0.0]
        h = h_surface_nm * 1e-3
        z_end = depth_um
        while z[-1] + h < z_end - 0.5 * h:
            z.append(z[-1] + h)
            h *= 1.3
        z.append(z_end)
        z = np.unique(np.round(np.r_[np.array(z), d_nplus_um], 12))
        ys = np.r_[-z[::-1], np.linspace(0.0, t_ox_um, 3)[1:]]
        self.mesh = mesh = TensorMesh(ls(1.0) * xs, ls(1.0) * ys)
        self.reg = reg = Regions.from_function(mesh, lambda xc, yc: yc < 0)
        cc = mesh.cell_centres
        xc, yc = cc[:, 0] / ls(1.0), cc[:, 1] / ls(1.0)
        semi = reg.semi_cell
        nplus = semi & (yc > -d_nplus_um) & ((xc < x_s0) | (xc > x_d0))
        nd_cell = np.where(nplus, sc.dens_s(n_plus_cm3), 0.0)
        na_cell = np.where(semi & ~nplus, float(sc.dens_s(na_cm3)), 0.0)
        self.nd, self.na = reg.semi_average(nd_cell), reg.semi_average(na_cell)
        self.donors, self.acceptors = sc.donor_levels(), sc.acceptor_levels(na_cm3)
        n_tot = 0.5 * (self.nd[mesh.ea] + self.nd[mesh.eb] + self.na[mesh.ea] + self.na[mesh.eb])
        mu_perp = sc.mobility_cm2(n_tot * sc.n0, "n", "perp")
        mu_par = sc.mobility_cm2(n_tot * sc.n0, "n", "par")
        mu_edge = sc.mu_s(np.where(mesh.edge_dir == 0, mu_perp, mu_par))
        mu_p = float(sc.mu_s(sc.mobility_cm2(na_cm3, "p", "perp")))
        vsat_cm, beta = sc.saturation()
        vs = float(sc.vsat_s(vsat_cm))
        self.mob = Mobility(vsat_n=vs, vsat_p=vs, beta_n=beta, beta_p=beta)
        self.dp = DDParams(delta=sc.delta, mu_n=1.0, mu_p=mu_p, recombination=False)
        ec_ei = sc.vt * np.log(sc.n_c / sc.n_i)
        dit_cm = sc.reg.value("Dit_Ec_minus_0p2eV_NO_0001") * dit_scale
        self.itf = InterfaceCharge(
            q_f=float(sc.sheet_s(qf_cm2)), dit=float(sc.sheet_s(dit_cm * sc.vt)),
            e1=float(sc.energy_s(ec_ei - dit_width_ev)), e2=float(sc.energy_s(ec_ei)),
        )  # fmt: skip
        self.eps_ratio = sc.reg.value("eps_r_SiO2_thermal") / sc.eps_r
        self.v_fb_phys, self.vfb_code, self.psi_body = sc.vfb(na_cm3)
        self.model = AnisoMobilityDD(
            reg, self.nd - self.na, self.dp, float(np.log(sc.delta)), eps_ratio=self.eps_ratio,
            ionization=Ionization(self.nd, self.na, self.donors, self.acceptors),
            mobility=self.mob, interface=self.itf, mu_edge=mu_edge,
        )  # fmt: skip
        xy = mesh.node_xy
        sx = lambda v: ls(v)  # noqa: E731
        tol = 1e-9
        self.nodes = {
            "source": np.flatnonzero(
                (np.abs(xy[:, 0]) < tol) & (xy[:, 1] <= 0) & (xy[:, 1] >= -sx(d_nplus_um) - tol)
            ),
            "drain": np.flatnonzero(
                (np.abs(xy[:, 0] - sx(x_end)) < tol)
                & (xy[:, 1] <= 0)
                & (xy[:, 1] >= -sx(d_nplus_um) - tol)
            ),
            "body": np.flatnonzero(np.abs(xy[:, 1] + sx(depth_um)) < tol),
            "gate": np.flatnonzero(
                (np.abs(xy[:, 1] - sx(t_ox_um)) < tol)
                & (xy[:, 0] >= sx(gate0) - tol) & (xy[:, 0] <= sx(gate1) + tol)
            ),
        }  # fmt: skip
        self.x_mid_s = sx(mid)
        self.x_s0_s, self.x_d0_s = sx(x_s0), sx(x_d0)
        self.gate_x = (sx(gate0), sx(gate1))
        self.n_cols = mesh.nx + 1

    def bc(self, vg, vs, vd) -> Dirichlet:
        mod, n = self.model, self.nodes
        bc = Dirichlet()
        mod.ohmic(bc, n["source"], vs)
        mod.ohmic(bc, n["drain"], vd)
        mod.ohmic(bc, n["body"], 0.0)
        bc.gate(n["gate"], vg, self.vfb_code)
        return bc

    def equilibrium(self):
        bc0 = self.bc(0.0, 0.0, 0.0)
        mod = self.model
        return mod.solve_bc(mod.initial_guess(bc0), mod.equilibrium_bc(bc0), max_step=20.0)

    def drain_current_a_per_cm(self, sol):
        i = self.model.contact_current(sol, self.nodes["drain"])
        return float(abs(i) * Q * self.sc.n0 * self.sc.mu0 * self.sc.vt)
