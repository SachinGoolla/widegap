r"""Stage 13b: corner-aware avalanche-BV objective with a rounded body junction.

Spec: docs/benchmarks/stage13b_corner_aware.md (pre-registered; AUTHORITATIVE).  It fixes the
stage-13 specification flaw (a single planar cut was gamed while the JFET/body corner BV fell).
Everything of stage 13 (:mod:`widegap.adjoint_sic`, not modified here) is reused:

* **Rounded junctions.**  The p-body acceptor profile and the n+ source donor profile are
  products of Gaussian-straggle erfc steps with rounding length ``r`` (the standard deviation of
  the smoothed edge, ``F(s) = erfc(-s / (sqrt(2) r)) / 2``):

      N_A(x, y) = a_body  F(x - x_jfet) F(y + d_body)                       (lateral + vertical edge)
      N_D,src   = b_src   F(x - x_src0) F(x_src1 - x) F(y + d_src)

  (the top ``y = 0`` and the symmetry plane ``x = x_max`` are not edges).  The donors are NOT
  netted with the acceptors: nodes carry ``nd`` and ``na`` separately (compensated material with
  incomplete ionization), so the theta dependence stays ``nd(theta) = M q_D(theta)`` with the
  stage-13 analytic ``dnd/dtheta``; ``na`` does not depend on theta.  For consistency the JFET
  donor layer is the complement of the body lateral edge, ``n_jfet e^{theta_J} F(x_jfet - x)``
  (otherwise a sliver of abrupt lateral n/p wall would survive at the corner).  The node
  doping is the semiconductor control-volume average of the cell-centre values, as before.
* **Objective.**  ``BV_obj = -tau log sum_k exp(-BV_k / tau)`` over ``K`` fixed vertical cuts, each
  ``BV_k`` the root of the electron-initiated ionization integral ``I_n,k(V_D, theta) = 1`` and
  ``dBV_obj/dtheta = sum_k w_k dBV_k/dtheta`` with the soft-min weights
  ``w_k = exp(-BV_k/tau) / sum_j exp(-BV_j/tau)`` and the stage-13 implicit-function adjoint
  gradient per cut.
* **Dense independent check** (:func:`dense_bv`): the avalanche BV along many vertical lines
  (bilinear interpolation of the mesh potential, so the same physical positions on every mesh).

The zero-current node classification (which nodes carry ``phi_n = V_D``) is that of the theta = 0
device; the n-type nodes in the (rounded) source box, enlarged by ``2.5 r``, are tied to the
source (``phi_n = 0``).
"""

from __future__ import annotations

from unittest import mock

import numpy as np
from scipy.special import erfc, logsumexp

import widegap.adjoint_sic as _asic
from widegap.adjoint_sic import EvalError, SiCAdjoint, _AffineBC
from widegap.devices.dmosfet import DMOSFET
from widegap.fv.avalanche import ionization_integrals
from widegap.fv2d.mos import Dirichlet
from widegap.fv2d.physics import Ionization, neutral_psi
from widegap.sic import AnisoMobilityDD, DMOSDesign, SiC4H, SiCDMOS

SRC_MARGIN = 2.5  # source box enlargement (in units of r) for the phi_n = 0 classification


def smooth_step(s, r):
    """F(s) = erfc(-s / (sqrt(2) r)) / 2: 0 -> 1 across s = 0 with standard deviation r."""
    return 0.5 * erfc(-np.asarray(s, dtype=float) / (np.sqrt(2.0) * r))


# ------------------------------------------------------------------------------ device
class RoundedDMOSFET(DMOSFET):
    """:class:`DMOSFET` with erfc-product p-body and n+ source (rounding length ``r``, scaled)."""

    def __init__(self, prm=None, level: int = 0, r: float = 1.0):
        super().__init__(prm, level)
        p = self.prm
        self.r = float(r)
        cc = self.mesh.cell_centres
        xc, yc = cc[:, 0], cc[:, 1]
        semi = self.reg.semi_cell
        st = lambda s: smooth_step(s, self.r)  # noqa: E731
        self.m_upper = semi & (yc >= -p.d_body)
        a_c = p.a_body if p.a_body_b <= 0.0 else self.body_acceptor(yc)
        self.acc_cell = np.where(semi, a_c * st(xc - p.x_jfet) * st(yc + p.d_body), 0.0)
        self.src_cell = np.where(
            semi, p.b_src * st(xc - p.x_src0) * st(p.x_src1 - xc) * st(yc + p.d_src), 0.0
        )
        self.g_jfet = np.where(self.m_upper, st(p.x_jfet - xc), 0.0)

    def cell_donors(self, theta):
        th = np.asarray(theta, dtype=float)
        d = np.zeros(self.mesh.n_cells)
        d[self.m_drain] = self.prm.b_drain
        d[self.m_drift] = self.drift_level(th)[self.m_drift]
        d[self.m_csl] = self.prm.n_csl
        d += self.prm.n_jfet * np.exp(th[4]) * self.g_jfet
        d += self.src_cell
        return d

    def cell_acceptors(self):
        return self.acc_cell.copy()

    def cell_doping(self, theta):  # net (report/compat); the solver uses donors/acceptors apart
        return self.cell_donors(theta) - self.cell_acceptors()

    def cell_doping_jac(self, theta):
        th = np.asarray(theta, dtype=float)
        nd = self.drift_level(th)
        d = np.zeros((self.mesh.n_cells, 5))
        for k in range(4):
            d[:, k] = np.where(self.m_drift, nd * self.leg[k], 0.0)
        d[:, 4] = self.prm.n_jfet * np.exp(th[4]) * self.g_jfet
        return d


class RoundedSiCDMOS(SiCDMOS):
    """Stage-12 device on the rounded-junction doping (compensated nd / na node fields)."""

    def __init__(self, sc: SiC4H, design: DMOSDesign, level=0, r_um=0.1, dit_scale=1.0, qf_cm2=0.0):
        self.r_um = float(r_um)
        self.r_s = float(sc.length_s(r_um * 1e-4))
        factory = lambda prm, lev: RoundedDMOSFET(prm, lev, self.r_s)  # noqa: E731
        with mock.patch("widegap.sic.DMOSFET", factory):
            super().__init__(sc, design, level, dit_scale=dit_scale, qf_cm2=qf_cm2)
        # rebuild everything that was derived from the netted doping in the parent constructor
        dev = self.dev
        z = np.zeros(5)
        self.nd = dev.M @ dev.cell_donors(z)
        self.na = dev.M @ dev.cell_acceptors()
        mesh = dev.mesh
        nd_cm = (
            sc.n0
            * 0.5
            * (self.nd[mesh.ea] + self.nd[mesh.eb] + self.na[mesh.ea] + self.na[mesh.eb])
        )
        mu_perp = sc.mobility_cm2(nd_cm, "n", "perp")
        mu_par = sc.mobility_cm2(nd_cm, "n", "par")
        self.mu_edge = sc.mu_s(np.where(mesh.edge_dir == 0, mu_perp, mu_par))
        self.model = AnisoMobilityDD(
            dev.reg, self.nd - self.na, self.dp, float(np.log(sc.delta)),
            eps_ratio=self.prm.eps_ratio,
            ionization=Ionization(self.nd, self.na, self.donors, self.acceptors),
            mobility=self.mob, interface=self.itf, mu_edge=self.mu_edge,
        )  # fmt: skip
        self._count_solves()

    def _contact_psi(self, nodes):
        """Neutral (incomplete-ionization) potential at contact nodes, cached: the rounded doping
        makes every node's (nd, na) pair unique, so the parent's per-unique-pair root find in
        ``initial_guess`` (one brentq per node, ~0.2 s) would dominate every boundary-condition
        build.  Same values as ``model.ohmic``."""
        cache = self.__dict__.setdefault("_psi_cache", {})
        out = np.empty(len(nodes))
        for j, i in enumerate(np.asarray(nodes, dtype=int)):
            if i not in cache:
                cache[i] = neutral_psi(self.nd[i], self.na[i], np.exp(self.model.ld),
                                       self.donors, self.acceptors)  # fmt: skip
            out[j] = cache[i]
        return out

    def bc(self, vg, vs, vd) -> Dirichlet:
        n = self.dev.nodes
        bc = Dirichlet()
        for key, v in (("source_body", vs), ("drain", vd)):
            nodes = n[key]
            bc.psi(nodes, self._contact_psi(nodes) + v).phi(nodes, v)
        bc.gate(n["gate"], vg, self.vfb_code)
        return bc

    def offstate_bc(self, vd, vg=0.0):
        """As :meth:`SiCDMOS.offstate_bc`, but the source box (phi_n = 0 for n-type nodes) is
        enlarged by ``SRC_MARGIN r`` to hold the rounded n+ source."""
        dev, p = self.dev, self.prm
        bc = self.bc(vg, 0.0, vd)
        xy = dev.mesh.node_xy
        semi = dev.reg.semi_node
        m = SRC_MARGIN * self.r_s
        in_src = (
            (xy[:, 0] >= p.x_src0 - m - 1e-9)
            & (xy[:, 0] <= p.x_src1 + m + 1e-9)
            & (xy[:, 1] >= -p.d_src - m - 1e-9)
        )
        n_type = self.nd > self.na
        contact = np.zeros(len(xy), dtype=bool)
        for k in ("source_body", "drain"):
            contact[dev.nodes[k]] = True
        free = np.flatnonzero(semi & ~contact)
        bc.phi(free, np.where(n_type[free] & ~in_src[free], vd, 0.0), 0.0)
        return bc


# ------------------------------------------------------------------------- soft minimum
def softmin(bv, tau):
    """(-tau log sum exp(-bv/tau), weights w_k = softmax(-bv/tau))."""
    z = -np.asarray(bv, dtype=float) / tau
    lse = float(logsumexp(z))
    return -tau * lse, np.exp(z - lse)


# --------------------------------------------------------------------------- adjoint model
class CornerAdjoint(SiCAdjoint):
    """:class:`SiCAdjoint` on the rounded device with a soft-min multi-cut BV objective."""

    def __init__(self, sc, design, level=0, r_um=0.1, cuts_um=(0.0,), tau_v=5.0, **kw):
        def factory(sc_, design_, level_, dit_scale=1.0, qf_cm2=0.0):
            return RoundedSiCDMOS(sc_, design_, level_, r_um, dit_scale=dit_scale, qf_cm2=qf_cm2)

        with mock.patch.object(_asic, "SiCDMOS", factory):
            super().__init__(sc, design, level, **kw)
        self.r_um, self.tau_v = float(r_um), float(tau_v)
        self.cuts_um = [float(c) for c in cuts_um]
        self.cuts, self.cut_x = [], []
        for c in self.cuts_um:
            x, nodes = self.base.column_nodes(float(sc.length_s(c * 1e-4)))
            self.cut_x.append(x)
            self.cuts.append(nodes)
        self.cut_x_um = [float(sc.length_cm(x)) * 1e4 for x in self.cut_x]
        self.caches: list = [None] * len(self.cuts)

    def set_theta(self, theta):
        th = np.asarray(theta, dtype=float).copy()
        dev, sc, mesh, mod = self.dev, self.sc, self.mesh, self.mod
        nd = dev.M @ dev.cell_donors(th)
        na = dev.M @ dev.cell_acceptors()
        self.nd, self.na = nd, na
        mod.nd, mod.na = nd.copy(), na.copy()
        mod.ion.nd, mod.ion.na = mod.nd, mod.na
        self.dnd_dtheta = dev.M @ dev.cell_doping_jac(th)
        n_cm = sc.n0 * 0.5 * (nd[mesh.ea] + nd[mesh.eb] + na[mesh.ea] + na[mesh.eb])
        mp_, dmp = self._mobility_cm2_and_deriv(n_cm, "perp")
        mq_, dmq = self._mobility_cm2_and_deriv(n_cm, "par")
        xdir = mesh.edge_dir == 0
        mod.mu_edge = np.asarray(sc.mu_s(np.where(xdir, mp_, mq_)), dtype=float)
        self.dmu_dnd = 0.5 * sc.n0 / sc.mu0 * np.where(xdir, dmp, dmq)
        self.flips = int(np.sum((nd > na) != self.n_type0))
        self.stats["flips_max"] = max(self.stats["flips_max"], self.flips)
        self.theta = th

    def bv_cuts(self, theta, grad=True, warm=True):
        """BV_k (and dBV_k/dtheta, dBV_k/dnd) for every objective cut."""
        outs = []
        for k, nodes in enumerate(self.cuts):
            cache = self.caches[k]
            if (
                cache is None and warm
            ):  # start from the nearest already-solved cut (saves a cold march)
                done = [j for j in range(len(self.cuts)) if self.caches[j] is not None]
                if done:
                    cache = self.caches[min(done, key=lambda j: abs(j - k))]
            self.cut_nodes, self.bv_cache = nodes, cache
            outs.append(self.bv(theta, grad=grad, warm=warm))
            self.caches[k] = self.bv_cache
        return outs

    def objective(self, theta, grad=True, warm=True):
        """Soft-min BV objective [V] and its gradient dBV_obj/dtheta (and dBV_obj/dnd)."""
        outs = self.bv_cuts(theta, grad=grad, warm=warm)
        bvk = np.array([o["BV"] for o in outs])
        obj, w = softmin(bvk, self.tau_v)
        res = {"BV": obj, "bv_k": bvk, "w": w, "outs": outs, "flips": self.flips,
               "limiting_cut": int(np.argmin(bvk))}  # fmt: skip
        if grad:
            res["grad"] = sum(wk * o["grad"] for wk, o in zip(w, outs, strict=True))
            res["g_nd"] = sum(wk * o["g_nd"] for wk, o in zip(w, outs, strict=True))
            res["grad_k"] = np.array([o["grad"] for o in outs])
        return res

    def reset_caches(self):
        self.caches = [None] * len(self.cuts)
        self.bv_cache = None
        self.X_on = None
        self.X_eq = None


# --------------------------------------------------------------- dense independent cut set
def line_operator(mesh, nodes_fn_x, x):
    """Rows (left_nodes, right_nodes, w) to interpolate nodal data bilinearly on the vertical line
    at coordinate ``x`` (columns of the tensor mesh): value = (1-w) f[left] + w f[right]."""
    xs = mesh.x
    i = int(np.clip(np.searchsorted(xs, x, side="right") - 1, 0, len(xs) - 2))
    w = float((x - xs[i]) / (xs[i + 1] - xs[i]))
    return xs[i], xs[i + 1], w


class DenseCuts:
    """Many vertical cut lines at fixed physical x (bilinear interpolation of psi)."""

    def __init__(self, ev: CornerAdjoint, x_um):
        sc, mesh = ev.sc, ev.mesh
        self.x_um = np.asarray(x_um, dtype=float)
        self.ev = ev
        self.lines = []
        for xu in self.x_um:
            x = float(sc.length_s(xu * 1e-4))
            xl, xr, w = line_operator(mesh, None, x)
            _, nl = ev.base.column_nodes(xl)
            _, nr = ev.base.column_nodes(xr)
            depth = -mesh.node_xy[nl, 1]
            self.lines.append((nl, nr, w, depth))

    def integrals(self, X):
        """I_n of every line from the full state vector ``X`` (psi = X[0::3])."""
        psi = np.asarray(X)[0::3]
        out = np.empty(len(self.lines))
        ch = self.ev.ch
        for j, (nl, nr, w, depth) in enumerate(self.lines):
            p = (1.0 - w) * psi[nl] + w * psi[nr]
            e = -np.diff(p) / np.diff(depth)
            out[j] = ionization_integrals(depth, e, ch)[0]
        return out


def march_states(ev: CornerAdjoint, theta, i_fun, v_step_v=20.0, v_max_v=2500.0, extend=1.0):
    """Zero-current off-state march (tangent predictor + Newton, step halving) from V_D = 0.

    Stores (v_s, X, tangent) per state; stops ``extend`` times the voltage at which
    ``i_fun(X).max()`` first reaches 1."""
    ev.set_theta(theta)
    sc = ev.sc
    eq = ev._eq()
    rows0, vals0 = ev._off_bc(0.0)
    s0 = ev.mod.solve_bc(eq.X, _AffineBC(rows0, vals0), max_step=20.0)
    if not s0.report.converged:
        raise EvalError("off-state zero-bias solve did not converge")
    v, x = 0.0, s0.X.copy()
    h = float(sc.volt_s(v_step_v))
    states = []
    stop_s = None
    while True:
        rows, dvals = ev._dvals(v)
        t = ev._tangent(x, rows, dvals)
        states.append((v, x.copy(), t))
        if stop_s is None and float(np.max(i_fun(x))) >= 1.0:
            stop_s = v * extend
        if stop_s is not None and v >= stop_s - 1e-12:
            return states
        if v > float(sc.volt_s(v_max_v)):
            raise EvalError("I = 1 not reached below the voltage ceiling")
        while True:
            v_try = v + h
            s = ev._off_solve(v_try, x + h * t)
            if s is None:
                h *= 0.5
                if h < 1e-3:
                    raise EvalError(f"off-state march failed at V_D = {v:.1f} V_T")
                continue
            break
        v, x = v_try, s.X.copy()


def state_at(ev: CornerAdjoint, states, v_s):
    """Converged state at ``v_s`` from the nearest lower stored state (tangent predictor)."""
    k = max(i for i, s in enumerate(states) if s[0] <= v_s + 1e-12)
    v0, x0, t0 = states[k]
    if abs(v_s - v0) < 1e-12:
        return x0
    s = ev._off_solve(v_s, x0 + (v_s - v0) * t0)
    if s is None:
        raise EvalError(f"refinement solve failed at V_D = {v_s:.3f} V_T")
    return s.X


def dense_bv(
    ev: CornerAdjoint,
    theta,
    dense: DenseCuts,
    v_step_v=20.0,
    tol_v=0.02,
    cap=2.0,
    loose_tol_v=1.0,
    near_frac=0.03,
):
    """BV of every dense vertical cut (root of I_n,k = 1) and their minimum.

    A march to the first crossing (and on to ``cap`` times that voltage so most cuts bracket),
    then bracket refinement by log-linear interpolation of ln I_k on additional solves; each
    solve updates the brackets of all cuts.  Cuts that never reach I = 1 below the cap are
    returned as NaN (reported).  Cuts within ``near_frac`` of the smallest bracketed BV are
    refined to ``tol_v`` [V], the others to ``loose_tol_v`` [V]."""
    sc = ev.sc
    states = march_states(ev, theta, dense.integrals, v_step_v=v_step_v, extend=cap)
    n = len(dense.lines)
    pts = [(v, np.log(np.maximum(dense.integrals(x), 1e-300))) for v, x, _ in states]
    blo = np.zeros((n, 2))  # bracket ends (v_s, ln I): ln I < 0 below, >= 0 above
    bhi = np.zeros((n, 2))
    blo[:] = np.nan
    bhi[:] = 0.0
    for v, lni in pts:
        for j in range(n):
            if lni[j] < 0.0:
                if np.isnan(blo[j, 0]) or v > blo[j, 0]:
                    blo[j] = (v, lni[j])
            else:
                if bhi[j, 0] == 0.0 or v < bhi[j, 0]:
                    bhi[j] = (v, lni[j])
    tol_s, loose_s = float(sc.volt_s(tol_v)), float(sc.volt_s(loose_tol_v))
    resolved_all = False
    n_solves = 0
    for _ in range(400):
        ok = (bhi[:, 0] > 0.0) & ~np.isnan(blo[:, 0])
        width = np.where(ok, bhi[:, 0] - blo[:, 0], 0.0)
        vmin = bhi[ok, 0].min()
        tol_j = np.where(bhi[:, 0] <= vmin * (1.0 + near_frac), tol_s, loose_s)
        todo = np.flatnonzero(ok & (width > tol_j))
        if len(todo) == 0:
            resolved_all = True
            break
        j = todo[np.argmin(bhi[todo, 0])]  # smallest BV first
        (v0, l0), (v1, l1) = blo[j], bhi[j]
        vn = v0 + (0.0 - l0) * (v1 - v0) / (l1 - l0)
        vn = float(np.clip(vn, v0 + 0.02 * (v1 - v0), v1 - 0.02 * (v1 - v0)))
        x = state_at(ev, states, vn)
        n_solves += 1
        lni = np.log(np.maximum(dense.integrals(x), 1e-300))
        for q in range(n):
            if lni[q] < 0.0:
                if np.isnan(blo[q, 0]) or vn > blo[q, 0]:
                    blo[q] = (vn, lni[q])
            else:
                if bhi[q, 0] == 0.0 or vn < bhi[q, 0]:
                    bhi[q] = (vn, lni[q])
    bv = np.full(n, np.nan)
    ok = (bhi[:, 0] > 0.0) & ~np.isnan(blo[:, 0])
    vv = blo[ok, 0] + (0.0 - blo[ok, 1]) * (bhi[ok, 0] - blo[ok, 0]) / (bhi[ok, 1] - blo[ok, 1])
    bv[ok] = sc.volt_phys(vv)
    jm = int(np.nanargmin(bv))
    return {"bv": bv, "bv_min": float(bv[jm]), "argmin": jm, "x_argmin_um": float(dense.x_um[jm]),
            "resolved_all": bool(resolved_all), "n_refine_solves": n_solves,
            "n_unbracketed": int((~ok).sum()), "states": states}  # fmt: skip


def solve_at_voltage(ev: CornerAdjoint, theta, v_volt, v_step_v=20.0):
    """Electrostatic state (X) at drain voltage ``v_volt`` [V]: march, then a final solve."""
    ev.set_theta(theta)
    sc = ev.sc
    target = float(sc.volt_s(v_volt))
    eq = ev._eq()
    rows0, vals0 = ev._off_bc(0.0)
    s0 = ev.mod.solve_bc(eq.X, _AffineBC(rows0, vals0), max_step=20.0)
    if not s0.report.converged:
        raise EvalError("off-state zero-bias solve did not converge")
    v, x = 0.0, s0.X.copy()
    h = float(sc.volt_s(v_step_v))
    while v < target - 1e-12:
        h_try = min(h, target - v)
        rows, dvals = ev._dvals(v)
        t = ev._tangent(x, rows, dvals)
        s = ev._off_solve(v + h_try, x + h_try * t)
        if s is None:
            h *= 0.5
            if h < 1e-3:
                raise EvalError(f"off-state march failed at V_D = {v:.1f} V_T")
            continue
        v, x = v + h_try, s.X.copy()
    return x


__all__ = [
    "CornerAdjoint",
    "DenseCuts",
    "RoundedDMOSFET",
    "RoundedSiCDMOS",
    "dense_bv",
    "march_states",
    "smooth_step",
    "softmin",
    "solve_at_voltage",
    "state_at",
]
