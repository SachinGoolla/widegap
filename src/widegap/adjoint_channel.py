r"""Stage 17: exact adjoint of the channel-model R_on,sp at the Picard fixed point.

Spec: docs/benchmarks/stage17_channel_aware_opt.md (pre-registered; AUTHORITATIVE).

The stage-16 forward model (:func:`widegap.channel_mobility.picard_solve`) converges to a fixed
point ``F(X, mu(X), theta) = 0`` where ``mu`` is the per-edge low-field electron mobility factor
of :class:`widegap.sic.AnisoMobilityDD` (``mu_edge``).  On the selected channel edges
(:class:`~widegap.channel_mobility.ChannelMap`) ``mu`` is a function of the electron densities of
the nodes within ``d_ch`` of the two columns through the Noguchi mobility::

    N_S,i   = trapezoid_k  w_k n0 n_k                    (cm^-2, mesh nodes of column i)
    mu_i    = mu_Noguchi(N_A,i, N_S,i)                   (clamps: N_S >= 1e11, E_eff >= 1e-2)
    mu_edge = mu_s( 2 mu_a mu_b / (mu_a + mu_b) )        (harmonic mean of the two columns)

so ``dmu/dX`` (``mu_X``) is sparse and analytic: ``dn_k/dpsi_k = n_k``, ``dn_k/dphi_n,k = -n_k``,
the Noguchi chain rule below, and the harmonic-mean partials.  With the on-state current
``I = sum of the drain-node rows of F`` and ``R = V_DS x / I`` the exact gradient is

    (F_X + F_mu mu_X)^T lam = dI/dX|_total,       dI/dtheta = (e_drain - lam)^T F_theta,
    dR/dtheta = -(R / I) dI/dtheta,

``F_mu`` is the derivative of the electron-continuity rows with respect to ``mu_edge`` at fixed
potentials (the velocity-saturation factor included: ``d/dmu [mu (1 + (mu r f)^beta)^(-1/beta)] =
(1 + x)^(-1/beta - 1)``).  ``F_theta`` is the existing bulk-mobility ``dF/dnd`` with
``d mu/d nd`` set to ZERO on the channel edges: their mobility depends on theta only through X
(``N_dep`` uses the surface ``N_A``, which the five design variables do not change; checked and
reported).  The *lagged* gradient solves the same system without ``F_mu mu_X`` (R2).

Frozen classification (J4): the channel columns (p-type surface at theta = 0), their node sets and
the selected edges are those of the theta = 0 device; the per-theta guards (surface sign of
``N_D - N_A``, ``N_S`` / ``E_eff`` distance from the clamps in the gated columns) are recorded at
every evaluation.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from widegap.adjoint2d import solve_transposed
from widegap.adjoint_sic import EvalError
from widegap.channel_mobility import (
    EEFF_MIN_MV_CM,
    NS_MIN_CM2,
    ChannelMap,
    NoguchiParams,
    depletion_sheet,
    eps_sic_f_per_cm,
    noguchi_mobility,
    picard_solve,
)
from widegap.materials.ionization import Q

GUARD_FACTOR = 2.0
D_CH_UM = 0.020  # 20 nm (stage 16 / spec)


def noguchi_dmu_dns(sc, p: NoguchiParams, na_cm3, ns_cm2):
    """(mu, d mu / d N_S) [cm^2/Vs, cm^2/Vs per cm^-2] of :func:`noguchi_mobility` (clamped).

    The derivative is zero where the N_S clamp is active (``ns_in < 1e11``) and the E_eff clamp
    contributes nothing where it is active."""
    na = np.asarray(na_cm3, dtype=float)
    ns_in = np.asarray(ns_cm2, dtype=float)
    r = noguchi_mobility(sc, p, na, ns_in)
    ns, e = r["ns_used"], r["e_eff"]
    act_ns = (ns_in >= NS_MIN_CM2).astype(float)
    de_dns = (Q / eps_sic_f_per_cm(sc)) * p.eta * 1e-6 * act_ns * (~r["clamp_eeff"])
    dmph = p.ph_exp * r["mu_ph"] / e * de_dns
    dmsr = p.sr_exp * r["mu_sr"] / e * de_dns
    dmc = r["gamma"] * r["mu_c"] / ns * act_ns
    mu = r["mu"]
    dmu = mu**2 * (dmph / r["mu_ph"] ** 2 + dmsr / r["mu_sr"] ** 2 + dmc / r["mu_c"] ** 2)
    return mu, dmu, r


class MuJacobian:
    """Sparse ``d mu_edge / dX`` (selected channel edges x 3N) at a converged state."""

    def __init__(self, cmap: ChannelMap, p: NoguchiParams):
        self.cmap, self.p = cmap, p
        sc, xy = cmap.sc, cmap.mesh.node_xy
        pos, node, wt = [], [], []
        for k, i in enumerate(cmap.cols):
            nodes = cmap.col_nodes[int(i)]
            depth = -xy[nodes, 1] * sc.l_d  # cm, increasing
            w = np.zeros(len(nodes))
            if len(nodes) > 1:
                dd = np.diff(depth)
                w[:-1] += 0.5 * dd
                w[1:] += 0.5 * dd
            pos.append(np.full(len(nodes), k))
            node.append(nodes)
            wt.append(w)
        self.pos, self.node = np.concatenate(pos), np.concatenate(node)
        self.wt = np.concatenate(wt) * sc.n0  # cm^-2 per scaled n
        colpos = np.full(cmap.mesh.nx + 1, -1, dtype=int)
        colpos[cmap.cols] = np.arange(len(cmap.cols))
        self.ka, self.kb = colpos[cmap.edge_ca], colpos[cmap.edge_cb]

    def matrix(self, n_scaled):
        """(mu_X, dict) with mu_X of shape (n_channel_edges, 3 N)."""
        cmap, sc, p = self.cmap, self.cmap.sc, self.p
        nn = 3 * cmap.mesh.n_nodes
        ns_all = cmap.sheet_density(n_scaled)
        c = cmap.cols
        mu, dmu, res = noguchi_dmu_dns(sc, p, cmap.col_na_cm3[c], ns_all[c])
        wn = self.wt * np.asarray(n_scaled)[self.node]  # dN_S/dpsi (= -dN_S/dphi_n)
        rows = np.r_[self.pos, self.pos]
        cols = np.r_[3 * self.node, 3 * self.node + 1]
        d_ns = sp.csr_matrix((np.r_[wn, -wn], (rows, cols)), shape=(len(c), nn))
        ma, mb = mu[self.ka], mu[self.kb]
        s2 = (ma + mb) ** 2
        ca, cb = 2.0 * mb**2 / s2 / sc.mu0, 2.0 * ma**2 / s2 / sc.mu0
        ne = len(self.ka)
        idx = np.arange(ne)
        e_mat = sp.csr_matrix(
            (np.r_[ca * dmu[self.ka], cb * dmu[self.kb]],
             (np.r_[idx, idx], np.r_[self.ka, self.kb])),
            shape=(ne, len(c)),
        )  # fmt: skip
        return sp.csr_matrix(e_mat @ d_ns), {"mu": mu, "ns": ns_all[c], "res": res}


def fmu_matrix(mod, X, edges):
    """``F_mu``: (3N x n_edges) derivative of the electron rows with respect to ``mu_edge[edges]``."""
    mesh, mob = mod.mesh, mod.mob
    a, b = mod.a[edges], mod.b[edges]
    _, pn, _ = mod.split(X)
    _, B, _, A, E, _, _ = mod._flux_parts(X)
    f = np.abs(pn[mod.b] - pn[mod.a]) / mesh.h
    x = (mod.prm.mu_n * mod.mu_edge * f / mob.vsat_n) ** mob.beta_n
    dm = (1.0 + x) ** (-1.0 / mob.beta_n - 1.0)
    dj = (dm * mod.prm.mu_n * mod.gw * B * A * E)[edges]
    k = np.arange(len(edges))
    return sp.csr_matrix(
        (np.r_[dj, -dj], (np.r_[3 * a + 1, 3 * b + 1], np.r_[k, k])),
        shape=(3 * mesh.n_nodes, len(edges)),
    )


def _bc_matrix(J, rows):
    """Dirichlet-row modification used by ``jacobian_bc`` (rows -> identity)."""
    J = sp.csr_matrix(J)
    mask = np.ones(J.shape[0])
    mask[rows] = 0.0
    ident = sp.csr_matrix((np.ones(len(rows)), (rows, rows)), shape=J.shape)
    return sp.csc_matrix(sp.diags(mask) @ J + ident)


class ChannelRon:
    """Channel-model R_on,sp of a theta-parametrized device with exact and lagged adjoint gradients.

    ``ev`` is a :class:`widegap.adjoint_sic.SiCAdjoint` (the 13d ``OxideAdjoint`` at the on-state
    mesh level).  The channel map is frozen at theta = 0.  Every evaluation follows the stage-16
    protocol: bulk on-state at theta, then the Picard loop from that state (no mobility carried
    across thetas).  ``guards`` collects the J4 record of every evaluation.
    """

    def __init__(self, ev, params: NoguchiParams, d_ch_um=D_CH_UM):
        self.ev, self.P = ev, params
        ev.set_theta(np.zeros(5))
        self.cmap = ChannelMap(ev.sc, ev.mod, ev.nd, ev.na, d_ch_um, gate_x_max=ev.prm.gate_x1)
        self.mu_jac = MuJacobian(self.cmap, params)
        self.col_na0 = self.cmap.col_na_cm3.copy()
        self.sign0 = self._surface_sign(ev.nd, ev.na)
        self.guards: list = []
        self.stats = {"evals": 0, "picard_iterations": 0, "failed_picard": 0, "failed_evals": 0}
        self.picard_log: list = []

    # -------------------------------------------------------------------------- helpers
    def _surface_sign(self, nd, na):
        cm = self.cmap
        s = np.where(cm.surf_node >= 0, cm.surf_node, 0)
        return np.where(cm.surf_node >= 0, np.sign(np.asarray(nd)[s] - np.asarray(na)[s]), 0.0)

    def guard_record(self, theta, r):
        cm, ev = self.cmap, self.ev
        sign = self._surface_sign(ev.nd, ev.na)
        s = np.where(cm.surf_node >= 0, cm.surf_node, 0)
        na_now = np.asarray(ev.na)[s] * ev.sc.n0
        res = r["columns"]
        gated = cm.col_gated[cm.cols]
        ns_g = r["ns_all"][cm.cols][gated]
        e_raw = (
            (Q / eps_sic_f_per_cm(ev.sc))
            * (res["n_dep"] + self.P.eta * np.maximum(r["ns_all"][cm.cols], NS_MIN_CM2))
            * 1e-6
        )
        e_g = e_raw[gated]
        ns_all_c = r["ns_all"][cm.cols]
        near_ns = (ns_g > NS_MIN_CM2 / GUARD_FACTOR) & (ns_g < NS_MIN_CM2 * GUARD_FACTOR)
        near_e = (e_g > EEFF_MIN_MV_CM / GUARD_FACTOR) & (e_g < EEFF_MIN_MV_CM * GUARD_FACTOR)
        # ungated p columns (body contact): kink of the clamp inside the active edge set
        ung = ~gated
        kink_ung = (
            (ns_all_c[ung] > NS_MIN_CM2 / GUARD_FACTOR)
            & (ns_all_c[ung] < NS_MIN_CM2 * GUARD_FACTOR)
        ).sum()
        rec = {
            "theta": np.asarray(theta, float).tolist(),
            "sign_changes": int(np.sum((sign != self.sign0) & (self.sign0 != 0))),
            "gated_near_ns_clamp": int(near_ns.sum()),
            "gated_near_eeff_clamp": int(near_e.sum()),
            "gated_ns_min": float(ns_g.min()), "gated_eeff_raw_min": float(e_g.min()),
            "n_gated": int(gated.sum()),
            "ungated_near_ns_clamp": int(kink_ung),
            "ungated_ns_clamped": int((ns_all_c[ung] < NS_MIN_CM2).sum()),
            "na_surface_max_rel_change": float(
                np.max(np.abs(na_now[cm.cols] / self.col_na0[cm.cols] - 1.0))
            ),
        }  # fmt: skip
        rec["ok"] = bool(
            rec["sign_changes"] == 0
            and rec["gated_near_ns_clamp"] == 0
            and rec["gated_near_eeff_clamp"] == 0
        )
        return rec

    # ----------------------------------------------------------------------- forward
    def solve(self, theta):
        """Forward channel on-state at theta.  Returns (picard result, bc, i_scaled)."""
        ev, cmap = self.ev, self.cmap
        ev.set_theta(theta)
        cmap.mu_bulk = np.asarray(ev.mod.mu_edge, dtype=float).copy()
        sol_b, bc = ev.on_state()
        ev.mod.mu_edge = cmap.mu_bulk.copy()
        r = picard_solve(cmap, self.P, sol_b.X, bc)
        self.stats["evals"] += 1
        self.stats["picard_iterations"] += r["iterations"]
        self.picard_log.append({"iterations": r["iterations"], "converged": r["converged"],
                                "final_rel": r["final_rel"]})  # fmt: skip
        if not r["converged"]:
            self.stats["failed_picard"] += 1
            self.stats["failed_evals"] += 1
            ev.stats["failed_evals"] += 1
            raise EvalError(
                f"channel Picard/Newton did not converge (picard={r['picard_converged']}, "
                f"newton={r['newton_all_converged']}, its={r['iterations']})"
            )
        i_sc = ev.mod.contact_current(r["sol"], ev.base.dev.nodes["drain"])
        return r, bc, i_sc

    def ron(self, theta, grad=True, lagged=False, detail=False):
        """R_on,ch [ohm cm^2] (and exact ``grad`` [dR/dtheta], ``grad_lagged``) at theta."""
        ev = self.ev
        r, bc, i_sc = self.solve(theta)
        i_a = ev.base.current_a_per_cm(abs(i_sc))
        R = ev.V_DS * ev.x_cm / i_a
        guard = self.guard_record(theta, r)
        self.guards.append(guard)
        out = {"R": float(R), "I_scaled": float(i_sc), "picard": r, "guard": guard,
               "iterations": r["iterations"]}  # fmt: skip
        if not grad:
            ev.mod.mu_edge = self.cmap.mu_bulk.copy()
            return out
        sol = r["sol"]
        X = sol.X
        mod, cmap = ev.mod, self.cmap
        rows, _ = mod._all_bc(bc)
        dn = ev.dev.nodes["drain"]
        sel = np.r_[3 * dn + 1, 3 * dn + 2]
        jcore = sp.csr_matrix(mod._jacobian_core(X))
        mu_x, mu_info = self.mu_jac.matrix(sol.n)
        f_mu = fmu_matrix(mod, X, cmap.edges)
        jtot = sp.csr_matrix(jcore + f_mu @ mu_x)
        fac = -R / i_sc

        def grad_of(jmat):
            g_i = np.asarray(jmat[sel].sum(axis=0)).ravel()
            lam = solve_transposed(_bc_matrix(jmat, rows), g_i)
            lam[rows] = 0.0
            w = -lam
            w[sel] += 1.0
            return lam, w

        saved = ev.dmu_dnd
        masked = saved.copy()
        masked[cmap.edges] = 0.0  # channel edges: mu depends on theta only through X
        ev.dmu_dnd = masked
        try:
            lam, w = grad_of(jtot)
            g_nd = ev.dF_dnd(w, X)
            out["lam"], out["w"] = lam, w
            out["g_nd"] = fac * g_nd
            out["grad"] = ev.dnd_dtheta.T @ out["g_nd"]
            if lagged or detail:
                lam_l, w_l = grad_of(jcore)
                g_nd_l = ev.dF_dnd(w_l, X)
                out["g_nd_lagged"] = fac * g_nd_l
                out["grad_lagged"] = ev.dnd_dtheta.T @ out["g_nd_lagged"]
                out["w_lagged"] = w_l
        finally:
            ev.dmu_dnd = saved
        out["mu_info"] = mu_info
        out["state"] = {"X": X, "bc": bc, "rows": rows, "sel": sel, "jcore": jcore,
                        "f_mu": f_mu, "mu_x": mu_x}  # fmt: skip
        ev.mod.mu_edge = self.cmap.mu_bulk.copy()
        return out


def mobility_part_theta(ev, zero_edges, w, X, k):
    """Per-edge ``w``-weighted ``dF/dtheta_k`` of the electron rows through the bulk mobility
    (edges in ``zero_edges`` excluded, e.g. the channel edges; None = all edges), and the per-node
    ionization part, for the R2 JFET-sensitivity map.

    Returns (edge_contrib [n_edges], node_contrib [n_nodes]) with
    ``sum(edge) + sum(node) = w^T dF/dtheta_k`` (explicit part, core rows)."""
    mod, mesh = ev.mod, ev.mesh
    a, b = mod.a, mod.b
    d = ev.dnd_dtheta[:, k]
    _, pn, _ = mod.split(X)
    _, B, _, A, E, _, _ = mod._flux_parts(X)
    f = np.abs(pn[b] - pn[a]) / mesh.h
    mob = mod.mob
    x = (mod.prm.mu_n * mod.mu_edge * f / mob.vsat_n) ** mob.beta_n
    dm = (1.0 + x) ** (-1.0 / mob.beta_n - 1.0)
    dmu = ev.dmu_dnd.copy()
    if zero_edges is not None:
        dmu[zero_edges] = 0.0
    dj = dm * mod.prm.mu_n * mod.gw * B * A * E * dmu
    s = w[1::3][a] - w[1::3][b]
    with np.errstate(all="ignore"):
        edge = np.where(s != 0.0, s * dj, 0.0) * (d[a] + d[b])
    # ionization part: dF_psi,i/dnd_i = -vol_i * (ionized fraction), reproduced via dF_dnd on
    # a zero-mobility-derivative copy
    saved = ev.dmu_dnd
    ev.dmu_dnd = np.zeros_like(saved)
    try:
        node = ev.dF_dnd(w, X) * d
    finally:
        ev.dmu_dnd = saved
    return edge, node


__all__ = [
    "ChannelRon",
    "MuJacobian",
    "fmu_matrix",
    "mobility_part_theta",
    "noguchi_dmu_dns",
    "depletion_sheet",
]
