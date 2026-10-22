"""Stage 9a: additional drift-diffusion physics on the stage-8.1 MOS core.

Three models, added to :class:`widegap.fv2d.mos.MOSDriftDiffusion` by the subclass
:class:`PhysicsMOSDD` (scaled units of docs/equations.md; X = (psi, phi_n, phi_p)
interleaved per node, n = delta exp(psi - phi_n), p = delta exp(phi_p - psi)).

1. Incomplete ionization out of equilibrium (log-space, stage-7 ``Level``s):

       N_D^+ = sum_l f_l N_D sigma(-(psi - phi_n + l_D)),
       N_A^- = sum_l f_l N_A sigma(-(phi_p - psi + l_A)),     sigma(z) = 1/(1+e^-z).

   (``expit`` is overflow-free; at phi = 0 this is the stage-7 ``IonizedJunction``.)
   An empty level tuple means "completely ionized".  With s = sigma(psi-phi_n+l),
   t = sigma(phi_p-psi+l_A) and rho_i = N_D^+ - N_A^-:
       d rho_i/d psi = -N_D s(1-s) - N_A t(1-t),
       d rho_i/d phi_n = +N_D s(1-s),   d rho_i/d phi_p = +N_A t(1-t).

2. Field-dependent mobility (Caughey-Thomas), per edge and carrier,

       m(F) = mu/mu_0 = [1 + (mu_0 F / v_sat)^beta]^(-1/beta),   F = |phi_b - phi_a| / h_e,

   so the edge flux is  j = m(F) j_0  with j_0 the constant-mobility SG flux.  The exact
   derivative with respect to the two endpoint quasi-Fermi potentials is
   dj/dphi_b = m dj_0/dphi_b + j_0 dm/dphi_b, dm/dphi_a = -dm/dphi_b, with

       dm/dphi_b = -(mu_0/v_sat)^beta F^(beta-1) [1 + (mu_0 F/v_sat)^beta]^(-1/beta-1) sign(Dphi)/h_e.

   F = 0: for beta > 1 the expression is 0; for beta = 1 the factor F^0 = 1 times
   sign(0) = 0 is taken (kink of |Dphi|), which is harmless because j_0 vanishes linearly
   with Dphi (expm1), so j_0 dm/dphi is continuous and tends to 0.  beta < 1 is rejected.
   This is the analytic limit and the Jacobian is continuous for all beta >= 1.

3. Interface charge: a sheet density (charge per unit interface length, same sign as
   rho) at semiconductor/oxide interface nodes,  Q = Q_f + Q_it(s), s = psi - phi_n,
   enters the Poisson row i as  F_psi,i -= L_i Q(s_i),  with L_i = half the sum of the
   lengths of the interface edges at node i (2D) or 1 (1D, unit cross-section).

   Trap charge, uniform D_it on [E1, E2] (energies in V_T relative to midgap, midgap =
   intrinsic level so that n = n_i e^{s}).  A trap at energy E above midgap has
   occupancy f(E) = 1/(1+exp(E - s)) (Fermi-Dirac w.r.t. the electron quasi-Fermi
   level, electron energy -psi).  States with E > 0 are acceptor-like (charge -1 when
   occupied), E < 0 donor-like (charge +1 when empty):

       Q_it = D_it [ -int_a^b f dE + int_{E1}^{c} (1 - f) dE ],  a = max(E1,0), b = E2,
                                                                  c = min(E2,0).
   With sp(x) = ln(1+e^x): int f dE = s - sp(s-E) (since d/dE[s - sp(s-E)] = 1/(1+e^{E-s}))
   and int (1-f) dE = sp(E-s), hence the closed form

       Q_it(s) = D_it [ sp(c-s) - sp(E1-s) - sp(s-a) + sp(s-b) ],
       dQ_it/ds = D_it [ -sig(c-s) + sig(E1-s) - sig(s-a) + sig(s-b) ]  (<= 0).
   (For a window straddling midgap, a=c=0; a window on one side of midgap is clipped:
   a = max(E1,0), b = max(E2,0), lo = min(E1,0), c = min(E2,0), with E1 -> lo above.)  sp is evaluated with ``logaddexp``.
   In the Poisson row  d F/d psi = -L Q',  d F/d phi_n = +L Q'.

Also here: exact MOS reference with Q_f + Q_it(psi_s) (``gate_voltage_it`` and
``psi_s_of_vg_it``), extending ``widegap.mos_reference`` without editing it.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
from scipy.optimize import brentq
from scipy.special import expit

from widegap import mos_reference as ref
from widegap.fv.drift_diffusion import DDParams, bernoulli_prime
from widegap.fv2d.mos import Dirichlet, MOSDriftDiffusion
from widegap.fv2d.regions import Regions
from widegap.materials.ionization import Level, neutral_potential

# ---------------------------------------------------------------- model functions


def softplus(x):
    return np.logaddexp(0.0, x)


def trap_charge(s, dit, e1, e2):
    """(Q_it, dQ_it/ds) for uniform D_it on [e1, e2] (closed form, see module doc)."""
    if e2 < e1:
        raise ValueError("require e1 <= e2")
    s = np.asarray(s, dtype=float)
    a, b = max(e1, 0.0), max(e2, 0.0)  # acceptor-like range [a, b] (empty if e2 <= 0)
    lo, c = min(e1, 0.0), min(e2, 0.0)  # donor-like range [lo, c] (empty if e1 >= 0)
    q = dit * (softplus(c - s) - softplus(lo - s) - softplus(s - a) + softplus(s - b))
    dq = dit * (-expit(c - s) + expit(lo - s) - expit(s - a) + expit(s - b))
    return q, dq


def ct_factor(dphi, h, mu0, vsat, beta):
    """Caughey-Thomas m = mu/mu0 and dm/d(dphi), dphi = phi_b - phi_a (module doc)."""
    if beta < 1.0:
        raise ValueError("beta must be >= 1")
    dphi = np.asarray(dphi, dtype=float)
    f = np.abs(dphi) / h
    r = mu0 / vsat
    x = (r * f) ** beta
    m = (1.0 + x) ** (-1.0 / beta)
    dmdf = -(r**beta) * f ** (beta - 1.0) * (1.0 + x) ** (-1.0 / beta - 1.0)
    return m, dmdf * np.sign(dphi) / h


@dataclass(frozen=True)
class Mobility:
    vsat_n: float
    vsat_p: float | None = None
    beta_n: float = 2.0
    beta_p: float = 2.0

    def __post_init__(self):
        if self.vsat_p is None:
            object.__setattr__(self, "vsat_p", self.vsat_n)


@dataclass(frozen=True)
class InterfaceCharge:
    q_f: float = 0.0
    dit: float = 0.0
    e1: float = -20.0
    e2: float = 20.0

    def charge(self, s):
        """(Q_f + Q_it(s), dQ/ds)."""
        q, dq = trap_charge(s, self.dit, self.e1, self.e2)
        return self.q_f + q, dq

    @property
    def magnitude(self):
        return abs(self.q_f) + abs(self.dit) * (self.e2 - self.e1)


@dataclass(eq=False)
class Ionization:
    nd: np.ndarray
    na: np.ndarray
    donors: tuple[Level, ...] = ()
    acceptors: tuple[Level, ...] = ()


def interface_length(reg: Regions) -> np.ndarray:
    """Interface face length per node (semiconductor/oxide boundary)."""
    m = reg.mesh
    out = np.zeros(m.n_nodes)
    if m.dim == 1:
        for k in range(2):
            ox = ~reg.semi_cell
            nodes = m.cell_nodes[ox, k]
            out[nodes[reg.semi_node[nodes]]] = 1.0
        return out
    sc = reg.semi_cell[m.cells]
    iface = (m.wc[:, 0] > 0) & (m.wc[:, 1] > 0) & (sc[:, 0] != sc[:, 1])
    np.add.at(out, m.ea[iface], 0.5 * m.h[iface])
    np.add.at(out, m.eb[iface], 0.5 * m.h[iface])
    return out


def neutral_psi(nd, na, delta, donors, acceptors):
    """Bulk charge-neutral psi with incomplete ionization (phi = 0)."""
    return neutral_potential(nd, na, delta, donors, acceptors)


# ---------------------------------------------------------------- the solver class


class PhysicsMOSDD(MOSDriftDiffusion):
    """Stage-8.1 MOS DD plus the three stage-9a models (each optional)."""

    def __init__(
        self,
        regions: Regions,
        c_nodes,
        params: DDParams,
        log_delta,
        eps_ratio=0.38,
        sources=None,
        ionization: Ionization | None = None,
        mobility: Mobility | None = None,
        interface: InterfaceCharge | None = None,
    ):
        super().__init__(regions, c_nodes, params, log_delta, eps_ratio, sources)
        self.ion, self.mob, self.itf = ionization, mobility, interface
        semi = regions.semi_node
        if ionization is not None:
            self.nd = np.where(semi, ionization.nd, 0.0)
            self.na = np.where(semi, ionization.na, 0.0)
            self.c = np.zeros(self.mesh.n_nodes)
        self.Lif = interface_length(regions) if interface is not None else None
        self.iface_nodes = np.flatnonzero(self.Lif > 0) if interface is not None else None

    # ----- model pieces ----------------------------------------------------------
    def _ion_terms(self, X):
        """(rho_i, d/dpsi, d/dphi_n, d/dphi_p) of N_D^+ - N_A^- (zeros if complete)."""
        psi, pn, pp = self.split(X)
        psi = self._psi_c(psi)
        ion = self.ion
        rho = np.zeros_like(psi)
        d0, d1, d2 = np.zeros_like(psi), np.zeros_like(psi), np.zeros_like(psi)
        if not ion.donors:
            rho += self.nd
        for lv in ion.donors:
            z = psi - pn + lv.log_gamma
            s = expit(z)
            rho += lv.fraction * self.nd * expit(-z)
            sd = lv.fraction * self.nd * s * (1.0 - s)
            d0 -= sd
            d1 += sd
        if not ion.acceptors:
            rho -= self.na
        for lv in ion.acceptors:
            w = pp - psi + lv.log_gamma
            t = expit(w)
            rho -= lv.fraction * self.na * expit(-w)
            td = lv.fraction * self.na * t * (1.0 - t)
            d0 -= td
            d2 += td
        return rho, d0, d1, d2

    def _mob_factors(self, X):
        _, pn, pp = self.split(X)
        a, b = self.a, self.b
        one, zero = np.ones(len(a)), np.zeros(len(a))
        if self.mob is None:
            return one, one, zero, zero
        mob, prm, h = self.mob, self.prm, self.mesh.h
        mn, kn = ct_factor(pn[b] - pn[a], h, prm.mu_n, mob.vsat_n, mob.beta_n)
        mp, kp = ct_factor(pp[b] - pp[a], h, prm.mu_p, mob.vsat_p, mob.beta_p)
        return mn, mp, kn, kp

    def fluxes(self, X):
        _, B, Bm, A, E, P, Fp = self._flux_parts(X)
        gw, prm = self.gw, self.prm
        mn, mp, _, _ = self._mob_factors(X)
        return mn * prm.mu_n * gw * B * A * E, -mp * prm.mu_p * gw * Bm * P * Fp

    # ----- residual / Jacobian (BC-free cores) --------------------------------------
    def _residual_core(self, X):
        psi, pn, pp = self.split(X)
        n, p = self.densities(X)
        a, b = self.a, self.b
        N = len(psi)
        fpsi = -(psi[b] - psi[a]) * self.G
        jn, jp = self.fluxes(X)
        F = np.zeros_like(X)
        Fpsi, Fn, Fp = F[0::3], F[1::3], F[2::3]
        for Fb, fe in ((Fpsi, fpsi), (Fn, jn), (Fp, jp)):
            Fb += np.bincount(a, fe, minlength=N)
            Fb -= np.bincount(b, fe, minlength=N)
        rho_fixed = self.c if self.ion is None else self._ion_terms(X)[0]
        Fpsi -= self.vol * (p - n + rho_fixed + self.s_psi)
        Fn -= self.vol * self.s_n
        Fp -= self.vol * self.s_p
        if self.itf is not None:
            i = self.iface_nodes
            q, _ = self.itf.charge(psi[i] - pn[i])
            Fpsi[i] -= self.Lif[i] * q
        return F

    def _jacobian_core(self, X):
        psi, pn, pp = self.split(X)
        N = len(psi)
        prm, g, G = self.prm, self.gw, self.G
        mu_n, mu_p = prm.mu_n, prm.mu_p
        a, b = self.a, self.b
        D, B, Bm, A, E, P, Fp = self._flux_parts(X)
        Bp_, Bmp = bernoulli_prime(D), bernoulli_prime(-D)
        mn, mp, kn, kp = self._mob_factors(X)
        jn0 = mu_n * g * B * A * E
        jp0 = -mu_p * g * Bm * P * Fp
        rows, cols, vals = [], [], []

        def add(r_node, r_var, c_node, c_var, v):
            rows.append(3 * r_node + r_var)
            cols.append(3 * c_node + c_var)
            vals.append(np.broadcast_to(v, r_node.shape))

        jn_psib = mn * mu_n * g * (Bp_ + B) * A * E
        jn_psia = -mn * mu_n * g * Bp_ * A * E
        jn_pnb = -mn * mu_n * g * B * A * (E + 1.0) + jn0 * kn
        jn_pna = mn * mu_n * g * B * A - jn0 * kn
        jp_psib = mp * mu_p * g * (Bmp + Bm) * P * Fp
        jp_psia = -mp * mu_p * g * Bmp * P * Fp
        jp_ppb = -mp * mu_p * g * Bm * P * (Fp + 1.0) + jp0 * kp
        jp_ppa = mp * mu_p * g * Bm * P - jp0 * kp
        for sign, row in ((1.0, a), (-1.0, b)):
            add(row, 0, a, 0, sign * G)
            add(row, 0, b, 0, -sign * G)
            add(row, 1, b, 0, sign * jn_psib)
            add(row, 1, a, 0, sign * jn_psia)
            add(row, 1, b, 1, sign * jn_pnb)
            add(row, 1, a, 1, sign * jn_pna)
            add(row, 2, b, 0, sign * jp_psib)
            add(row, 2, a, 0, sign * jp_psia)
            add(row, 2, b, 2, sign * jp_ppb)
            add(row, 2, a, 2, sign * jp_ppa)

        i = np.arange(N)
        n, p = self.densities(X)
        vol = self.vol
        d0, d1, d2 = vol * (p + n), -vol * n, -vol * p
        if self.ion is not None:
            _, r0, r1, r2 = self._ion_terms(X)
            d0, d1, d2 = d0 - vol * r0, d1 - vol * r1, d2 - vol * r2
        add(i, 0, i, 0, d0)
        add(i, 0, i, 1, d1)
        add(i, 0, i, 2, d2)
        if self.itf is not None:
            k = self.iface_nodes
            _, dq = self.itf.charge(psi[k] - pn[k])
            add(k, 0, k, 0, -self.Lif[k] * dq)
            add(k, 0, k, 1, self.Lif[k] * dq)
        return sp.csc_matrix(
            (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
            shape=(3 * N, 3 * N),
        )

    def residual_bc(self, X, rows, vals):
        F = self._residual_core(X)
        F[rows] = X[rows] - vals
        return F

    def jacobian_bc(self, X, rows):
        J = sp.csr_matrix(self._jacobian_core(X))
        mask = np.ones(J.shape[0])
        mask[rows] = 0.0
        ident = sp.csr_matrix((np.ones(len(rows)), (rows, rows)), shape=J.shape)
        return sp.csc_matrix(sp.diags(mask) @ J + ident)

    def residual(self, X, dir_nodes=None, dir_vals=None):  # pragma: no cover
        raise NotImplementedError("use residual_bc / solve_bc")

    def jacobian(self, X, dir_nodes=None):  # pragma: no cover
        raise NotImplementedError("use jacobian_bc / solve_bc")

    def solve(self, *a, **k):  # pragma: no cover
        raise NotImplementedError("use solve_bc")

    # ----- scaling, guesses, helpers ---------------------------------------------------
    def _volume_terms(self, X):
        out = super()._volume_terms(X)
        if self.ion is not None:
            out[0::3] += self.vol * (self.nd + self.na)
        if self.itf is not None:
            out[0::3] += self.Lif * self.itf.magnitude
        return out

    def initial_guess(self, bc: Dirichlet | None = None):
        if self.ion is None:
            return super().initial_guess(bc)
        n = self.mesh.n_nodes
        psi = np.zeros(n)
        pairs = np.unique(np.c_[self.nd, self.na][self.reg.semi_node], axis=0)
        for nd, na in pairs:
            sel = self.reg.semi_node & (self.nd == nd) & (self.na == na)
            psi[sel] = neutral_psi(nd, na, np.exp(self.ld), self.ion.donors, self.ion.acceptors)
        X = np.zeros(3 * n)
        X[0::3] = psi
        if bc is not None:
            rows, vals = self._all_bc(bc)
            X[rows] = vals
        return X

    def ohmic(self, bc: Dirichlet, nodes, v=0.0):
        """Ohmic contact with the (incomplete-ionization) neutral psi: psi_neutral + V."""
        nodes = np.asarray(nodes, dtype=int)
        if self.ion is None:
            return bc.ohmic(self.c, np.exp(self.ld), nodes, v)
        X = self.initial_guess()
        return bc.psi(nodes, X[0::3][nodes] + v).phi(nodes, v)

    def space_charge(self, sol):
        rho = self.c if self.ion is None else self._ion_terms(sol.X)[0]
        return float(np.sum(self.vol * (sol.p - sol.n + rho)))


# ---------------------------------------------------------------- exact MOS reference


def gate_voltage_it(psi_s, a, delta, t_ox, eps_r, vfb, itf: InterfaceCharge):
    """V_G for a surface potential: Gauss across the oxide with sheet charge
    Q_i = Q_f + Q_it(psi_surf), psi_surf = psi_b + psi_s (phi_n = 0):

        eps_s psi'(0-) - eps_ox psi'(0+) = Q_i  ->
        V_G - V_FB = psi_b + psi_s - (t_ox/eps_r) [Q_s(psi_s) + Q_i]."""
    psi_b = ref.bulk(a, delta)[0]
    qi = float(itf.charge(psi_b + psi_s)[0])
    return vfb + psi_b + psi_s - (t_ox / eps_r) * (ref.Qs(psi_s, a, delta) + qi)


def psi_s_of_vg_it(vg, a, delta, t_ox, eps_r, vfb, itf: InterfaceCharge):
    """psi_s(V_G) by bracketing (RHS strictly increasing: Q_s' <= 0, Q_it' <= 0)."""

    def f(u):
        return gate_voltage_it(u, a, delta, t_ox, eps_r, vfb, itf) - vg

    lo, hi = -1.0, 1.0
    while f(lo) > 0:
        lo *= 2.0
    while f(hi) < 0:
        hi *= 2.0
    return brentq(f, lo, hi, xtol=1e-14, rtol=4 * np.finfo(float).eps, maxiter=500)


__all__ = [
    "InterfaceCharge",
    "Ionization",
    "Level",
    "Mobility",
    "PhysicsMOSDD",
    "ct_factor",
    "gate_voltage_it",
    "interface_length",
    "neutral_psi",
    "psi_s_of_vg_it",
    "softplus",
    "trap_charge",
]
