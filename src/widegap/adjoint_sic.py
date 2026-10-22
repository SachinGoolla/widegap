r"""Stage 13: full-physics adjoint of the stage-12 4H-SiC DMOSFET (R_on,sp and avalanche BV).

Spec: docs/benchmarks/stage13_sic_dmosfet_optimization.md (pre-registered; AUTHORITATIVE).
Reuses the stage-12 device (:class:`widegap.sic.SiCDMOS`: geometry, sourced parameters,
bias conditions, planar cut) and the stage-9a analytic Jacobian.  Design variables
``theta = (theta_0..theta_3, theta_J)`` as in :class:`widegap.devices.dmosfet.DMOSFET`:
``N_D(y) = N_D0 exp(sum_k theta_k P_k(eta))`` in the drift and ``N_JFET = N_J0 exp(theta_J)``.

Dependence of the discrete system ``F(X; nd) = 0`` on the donor density ``nd`` (node values,
units of N_0; ``nd = M q_+(theta)``, so ``dnd/dtheta = M dq/dtheta`` is exact):

* Poisson rows through the *ionized* fraction (incomplete ionization, stage 9a),
  ``F_psi,i = ... - vol_i (p - n + N_D^+ - N_A^-)``,
  ``dF_psi,i/dnd_i = -vol_i sum_l f_l sigma(-(psi_i - phi_n,i + l_l))``;
* electron continuity rows through the per-edge low-field mobility factor ``u_e``
  (Caughey-Thomas/Ishikawa fit evaluated at ``N_e = n_0 (nd_a + nd_b + na_a + na_b)/2``,
  perp-c on x-edges and par-c on y-edges).  The edge flux is ``j_n = M_e j_0`` with
  ``M_e = u [1 + (u F / v_sat)^beta]^(-1/beta)`` (``F = |phi_b - phi_a|/h``), whose exact
  derivative is ``dM_e/du = (1 + x)^(-1/beta - 1)``, ``x = (u F/v_sat)^beta``; and
  ``du/dnd_a = du/dnd_b = (n_0/(2 mu_0)) dmu/dN`` with the analytic derivative of the fit.
  Hole fluxes do not depend on ``theta``.

Everything above is *analytic* (no finite difference in theta anywhere in the production path).
The unit tests compare the analytic ``dF/dtheta`` with a residual-only central difference
(step 1e-6) to ~1e-8 relative.

On-state (R_on,sp at V_GS = 15 V, V_DS = 0.1 V): ``R = V_DS x_max / I_D`` with ``I_D`` the drain
current per cm of depth.  The drain current is the sum of the drain-node rows of the core
electron + hole continuity residuals (all other edges cancel, no sources), hence
``dI/dX`` is the same row sum of the analytic Jacobian and
``dI/dtheta = dI/dtheta|_X - lam^T dF/dtheta`` with ``F_X^T lam = dI/dX`` (BC-modified Jacobian,
``dF/dtheta`` = 0 on Dirichlet rows).  The explicit term is nonzero because the flux contains
``u_e(theta)``.

Off-state (avalanche BV): zero-current 2D electrostatics with incomplete ionization (stage 12);
the BV is the root of ``I_n(V_D, theta) = 1`` of the electron-initiated ionization integral along
a FIXED vertical grid column (the stage-12 planar column ``x = x_bc``).  By the implicit
function theorem ``dBV/dtheta = -(dI_n/dtheta)/(dI_n/dV_D)``; both total derivatives come from
the adjoint of the electrostatic system (BC-modified Jacobian, ``F_X^T lam = dI_n/dpsi``,
``dI_n/dV_D = lam^T d(vals)/dV_D`` over the Dirichlet rows, because ``V_D`` enters through the
drain, ``phi_n`` of the n-type nodes, ...), and ``dI_n/dpsi`` along the cut is analytic
(:func:`ionization_n_grad`).  The root is found by a safeguarded Newton iteration in ``V_D`` on
``ln I_n`` using the tangent ``dX/dV_D`` (warm-started from the previous design, cold march
from the equilibrium otherwise).

Modelling choice: the node classification of ``offstate_bc`` (which nodes carry ``phi_n = V_D``)
is that of the theta = 0 device and is kept fixed, so BV is a smooth function of theta; the
number of nodes whose ``nd > na`` status differs from theta = 0 is reported at every evaluation.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
from scipy.special import expit as _expit

from widegap.adjoint2d import solve_transposed
from widegap.fv2d.mos import continue_bias
from widegap.fv2d.physics import Ionization
from widegap.sic import AnisoMobilityDD, DMOSDesign, SiC4H, SiCDMOS


# ----------------------------------------------------------------------- ionization integral
def _phi(w):
    """(e^w - 1)/w, series near 0 (same branch rule as ``avalanche._expm1_over``)."""
    w = np.asarray(w, dtype=float)
    out = np.empty_like(w)
    s = np.abs(w) < 1e-5
    ws = w[s]
    out[s] = 1.0 + ws / 2.0 + ws**2 / 6.0 + ws**3 / 24.0
    wl = w[~s]
    out[~s] = np.expm1(wl) / wl
    return out


def _dphi(w):
    """phi'(w) = (e^w (w - 1) + 1)/w^2 ; Taylor sum_{k>=1} k w^(k-1)/(k+1)! for |w| < 0.1."""
    w = np.asarray(w, dtype=float)
    out = np.empty_like(w)
    s = np.abs(w) < 0.1
    ws = w[s]
    acc = np.zeros_like(ws)
    term = np.ones_like(ws)  # w^(k-1)
    fact = 1.0  # (k+1)!
    for k in range(1, 16):
        fact *= k + 1
        acc += k * term / fact
        term = term * ws
    out[s] = acc
    wl = w[~s]
    with np.errstate(over="ignore"):
        out[~s] = (np.exp(wl) * (wl - 1.0) + 1.0) / wl**2
    return out


def ionization_n_grad(x, e_el, ch):
    """(I_n, dI_n/de) of ``widegap.fv.avalanche.ionization_integrals`` (electron-initiated
    integral with exact element-wise exponential quadrature) with respect to the element fields.

    With z_j = (a_n,j - a_p,j) h_j, A_j = sum_{i<=j} z_i, T_j = a_n,j h_j phi(-z_j) exp(-A_{j-1}),
    I = sum_j T_j:
        dI/dz_i = -a_n,i h_i phi'(-z_i) exp(-A_{i-1}) - sum_{j>i} T_j,
        dI/da_n,i = h_i phi(-z_i) exp(-A_{i-1}) + h_i dI/dz_i,   dI/da_p,i = -h_i dI/dz_i,
    and da/de = sign(e) da/d|e| (``ch.dalpha_*``)."""
    x = np.asarray(x, dtype=float)
    e = np.asarray(e_el, dtype=float)
    h = np.diff(x)
    an, ap = ch.alpha_n(e), ch.alpha_p(e)
    z = (an - ap) * h
    big_a = np.concatenate(([0.0], np.cumsum(z)))
    ea = np.exp(-big_a[:-1])
    ph = _phi(-z)
    t = an * h * ph * ea
    i_n = float(t.sum())
    tail = np.cumsum(t[::-1])[::-1] - t  # sum_{j>i} T_j
    dz = -an * h * _dphi(-z) * ea - tail
    dan = h * ph * ea + h * dz
    dap = -h * dz
    sg = np.sign(e)
    return i_n, dan * ch.dalpha_n(e) * sg + dap * ch.dalpha_p(e) * sg


# ------------------------------------------------------------------------------ helpers
class _AffineBC:
    """Minimal Dirichlet-like object (``arrays()``) carrying fixed rows/values."""

    def __init__(self, rows, vals):
        self.rows, self.vals = rows, vals

    def arrays(self):
        return self.rows, self.vals


def richardson_fd(fun, theta, u, eps=1e-3):
    """Central FD of ``fun`` along ``u`` at ``eps`` and ``eps/2`` and the Richardson value
    (4 FD(eps/2) - FD(eps))/3 (error O(eps^4)).  Returns (richardson, fd_eps, fd_half)."""
    theta = np.asarray(theta, dtype=float)

    def fd(e):
        return (fun(theta + e * u) - fun(theta - e * u)) / (2.0 * e)

    f1, f2 = fd(eps), fd(eps / 2.0)
    return (4.0 * f2 - f1) / 3.0, f1, f2


class EvalError(RuntimeError):
    """A forward solve did not converge (counted as a failed evaluation)."""


@dataclass
class OffState:
    theta: np.ndarray
    v_s: float  # BV in V_T
    bv_V: float
    X: np.ndarray  # electrostatic state at the last Newton iterate
    v_iter: float
    i_n: float
    di_dv: float
    rows: np.ndarray
    dvals: np.ndarray
    g_psi: np.ndarray
    n_iter: int
    cold: bool


def _counted(fn):
    """Count evaluations that ultimately fail (an EvalError escaping to the caller)."""

    def wrapper(self, *a, **k):
        try:
            return fn(self, *a, **k)
        except EvalError:
            self.stats["failed_evals"] += 1
            raise

    wrapper.__name__, wrapper.__doc__ = fn.__name__, fn.__doc__
    return wrapper


# --------------------------------------------------------------------------------- device
class SiCAdjoint:
    """theta-parametrized stage-12 DMOSFET with adjoint gradients of R_on,sp and the cut BV."""

    V_GS, V_DS = 15.0, 0.1  # V (spec section 1)

    def __init__(
        self,
        sc: SiC4H,
        design: DMOSDesign,
        level: int = 0,
        cut_x: float | None = None,
        ch=None,
        dit_scale=1.0,
        qf_cm2=0.0,
    ):
        self.sc, self.design, self.level = sc, design, level
        self.base = b = SiCDMOS(sc, design, level, dit_scale=dit_scale, qf_cm2=qf_cm2)
        self.dev, self.prm = b.dev, b.prm
        self.mesh = b.dev.mesh
        self.ch = ch if ch is not None else sc.niwa()
        self.col_x, self.cut_nodes = b.column_nodes(b.prm.x_bc if cut_x is None else cut_x)
        self.n_type0 = b.nd > b.na
        # theta-dependent model (separate from ``base.model`` which keeps the theta = 0 doping
        # for the contact boundary values: contacts only touch fixed-doping cells)
        m0 = b.model
        self.mod = AnisoMobilityDD(
            b.dev.reg, b.nd - b.na, b.dp, float(np.log(sc.delta)), eps_ratio=b.prm.eps_ratio,
            ionization=Ionization(b.nd.copy(), b.na.copy(), b.donors, b.acceptors),
            mobility=b.mob, interface=b.itf, mu_edge=m0.mu_edge.copy(),
        )  # fmt: skip
        self.stats = {
            "attempts": 0,
            "failed_attempts": 0,
            "failed_evals": 0,
            "warm_fallbacks": 0,
            "flips_max": 0,
        }
        orig = self.mod.solve_bc

        def counted(x0, bc, **kw):
            sol = orig(x0, bc, **kw)
            self.stats["attempts"] += 1
            self.stats["failed_attempts"] += int(not sol.report.converged)
            return sol

        self.mod.solve_bc = counted
        self.live_class = False
        self.theta = None
        self.X_eq = None
        self.X_on = None
        self.bv_cache: tuple[float, np.ndarray] | None = None
        self.x_cm = float(sc.length_cm(b.prm.x_max))
        self.set_theta(np.zeros(5))

    # ------------------------------------------------------------------ theta -> model
    def _mobility_cm2_and_deriv(self, n_cm3, orient):
        """Registry Ishikawa fit (electrons) and d mu/dN, identical to ``SiC4H.mobility_cm2``."""
        sc = self.sc
        r, sfx, t = sc.reg, f"{orient}_c", sc.t_k / 300.0
        mx, mn = r.value(f"mu_n_max_{sfx}"), r.value(f"mu_n_min_{sfx}")
        nref, al = r.value(f"N_n_ref_{sfx}"), r.value(f"alpha_n_{sfx}")
        d1, d2, d3 = (r.value(f"delta{k}_n_{sfx}") for k in (1, 2, 3))
        n = np.asarray(n_cm3, dtype=float)
        x = (n / nref) ** al * t**d2
        mu = (mx - mn) * t**d1 / (1.0 + x) + mn * t**d3
        with np.errstate(divide="ignore", invalid="ignore"):
            dmu = np.where(n > 0, -(mx - mn) * t**d1 * al * x / (n * (1.0 + x) ** 2), 0.0)
        return mu, dmu

    def set_theta(self, theta):
        th = np.asarray(theta, dtype=float).copy()
        dev, sc, mesh, mod = self.dev, self.sc, self.mesh, self.mod
        q = dev.cell_doping(th)
        nd = dev.M @ np.maximum(q, 0.0)
        na = dev.M @ np.maximum(-q, 0.0)
        self.nd, self.na = nd, na
        mod.nd = nd.copy()
        mod.na = na.copy()
        mod.ion.nd, mod.ion.na = mod.nd, mod.na
        self.dnd_dtheta = dev.M @ dev.cell_doping_jac(th)  # (n_nodes, 5)
        n_cm = sc.n0 * 0.5 * (nd[mesh.ea] + nd[mesh.eb] + na[mesh.ea] + na[mesh.eb])
        mp_, dmp = self._mobility_cm2_and_deriv(n_cm, "perp")
        mq_, dmq = self._mobility_cm2_and_deriv(n_cm, "par")
        xdir = mesh.edge_dir == 0
        mu_e = sc.mu_s(np.where(xdir, mp_, mq_))
        mod.mu_edge = np.asarray(mu_e, dtype=float)
        # d u_e / d nd_a = d u_e / d nd_b
        self.dmu_dnd = 0.5 * sc.n0 / sc.mu0 * np.where(xdir, dmp, dmq)
        self.flips = int(np.sum((nd > na) != self.n_type0))
        self.stats["flips_max"] = max(self.stats["flips_max"], self.flips)
        self.theta = th

    # -------------------------------------------------------------------- dF/dnd (analytic)
    def dF_dnd(self, w, X):
        """g_i = sum_r w_r dF_r/dnd_i of the CORE residual (no Dirichlet masking)."""
        mod, mesh = self.mod, self.mesh
        n = mesh.n_nodes
        psi, pn, _ = mod.split(X)
        psi_c = mod._psi_c(psi)
        out = np.zeros(n)
        frac = np.ones(n)
        if mod.ion.donors:
            frac = np.zeros(n)
            for lv in mod.ion.donors:
                frac += lv.fraction * _expit(-(psi_c - pn + lv.log_gamma))
        out += w[0::3] * (-mod.vol * frac)
        a, b = mod.a, mod.b
        _, B, _, A, E, _, _ = mod._flux_parts(X)
        f = np.abs(pn[b] - pn[a]) / mesh.h
        mu_n, mob = mod.prm.mu_n, mod.mob
        x = (mu_n * mod.mu_edge * f / mob.vsat_n) ** mob.beta_n
        dm = (1.0 + x) ** (-1.0 / mob.beta_n - 1.0)
        dj = dm * mu_n * mod.gw * B * A * E * self.dmu_dnd
        s = w[1::3][a] - w[1::3][b]
        with np.errstate(all="ignore"):
            c = np.where(s != 0.0, s * dj, 0.0)  # off-state: lam = 0 on these (Dirichlet) rows
        out += np.bincount(a, c, minlength=n) + np.bincount(b, c, minlength=n)
        return out

    # ------------------------------------------------------------------------ equilibrium
    def _eq(self):
        mod, base = self.mod, self.base
        bc = base.bc(0.0, 0.0, 0.0)
        eqbc = mod.equilibrium_bc(bc)
        for guess in (self.X_eq, None):
            x0 = mod.initial_guess(bc) if guess is None else guess
            s = mod.solve_bc(x0, eqbc, max_step=20.0)
            if s.report.converged:
                self.X_eq = s.X.copy()
                return s
        raise EvalError("equilibrium solve did not converge")

    # --------------------------------------------------------------------------- on-state
    def on_state(self):
        sc, base, mod = self.sc, self.base, self.mod
        vg, vd = float(sc.volt_s(self.V_GS)), float(sc.volt_s(self.V_DS))
        bc = base.bc(vg, 0.0, vd)
        s = None
        if self.X_on is not None:
            s = mod.solve_bc(self.X_on, bc, max_step=20.0)
        if s is None or not s.report.converged:
            eq = self._eq()
            out = continue_bias(
                mod, eq.X, lambda p: base.bc(p[0], 0.0, p[1]), [0.0, 0.0],
                [[0.0, vd], [vg, vd]], max_step=20.0,
            )  # fmt: skip
            s = out[-1][1]
            if not s.report.converged:
                raise EvalError("on-state solve did not converge")
        self.X_on = s.X.copy()
        return s, bc

    @_counted
    def ron(self, theta, grad=True, nodes=False):
        """R_on,sp [ohm cm^2] (and dR/dtheta, optionally dR/dnd at nodes) at ``theta``."""
        self.set_theta(theta)
        s, bc = self.on_state()
        mod, base = self.mod, self.base
        i_sc = mod.contact_current(s, base.dev.nodes["drain"])
        i_a = base.current_a_per_cm(abs(i_sc))
        r = self.V_DS * self.x_cm / i_a
        out = {"R": float(r), "I_scaled": float(i_sc), "I_A_per_cm": i_a, "X": s.X,
               "newton_its": s.report.iterations}  # fmt: skip
        if not grad:
            return out
        X = s.X
        rows, _ = mod._all_bc(bc)
        dn = self.dev.nodes["drain"]
        sel = np.r_[3 * dn + 1, 3 * dn + 2]
        jcore = sp.csr_matrix(mod._jacobian_core(X))
        g_i = np.asarray(jcore[sel].sum(axis=0)).ravel()
        lam = solve_transposed(mod.jacobian_bc(X, rows), g_i)
        lam[rows] = 0.0
        w = -lam
        w[sel] += 1.0
        g_nd = self.dF_dnd(w, X)  # dI_scaled/dnd
        fac = -r / i_sc  # dR = -(R/I) dI
        out["g_nd"] = fac * g_nd
        out["grad"] = self.dnd_dtheta.T @ out["g_nd"]
        return out

    # ------------------------------------------------------------------------- off-state
    def _cut_i(self, X):
        """(I_n along the fixed column, dI_n/dX)."""
        nodes = self.cut_nodes
        psi = X[0::3]
        depth = -self.mesh.node_xy[nodes, 1]
        dd = np.diff(depth)
        e = -np.diff(psi[nodes]) / dd
        i_n, d_e = ionization_n_grad(depth, e, self.ch)
        g = np.zeros(len(X))
        g[3 * nodes[1:]] += -d_e / dd
        g[3 * nodes[:-1]] += d_e / dd
        return i_n, g

    def _offstate_bc_live(self, v):
        """``SiCDMOS.offstate_bc`` with the n-type node set taken from the CURRENT doping
        (report-only diagnostic; the production path keeps the theta = 0 classification)."""
        base, p = self.base, self.prm
        bc = base.bc(0.0, 0.0, float(v))
        xy = self.mesh.node_xy
        semi = self.dev.reg.semi_node
        in_src = (
            (xy[:, 0] >= p.x_src0 - 1e-9)
            & (xy[:, 0] <= p.x_src1 + 1e-9)
            & (xy[:, 1] >= -p.d_src - 1e-9)
        )
        contact = np.zeros(len(xy), dtype=bool)
        for k in ("source_body", "drain"):
            contact[self.dev.nodes[k]] = True
        free = np.flatnonzero(semi & ~contact)
        n_type = self.nd > self.na
        bc.phi(free, np.where(n_type[free] & ~in_src[free], float(v), 0.0), 0.0)
        return bc

    def _off_bc(self, v):
        bc = self._offstate_bc_live(v) if self.live_class else self.base.offstate_bc(float(v))
        return self.mod._all_bc(bc)

    def _dvals(self, v):
        r0, v0 = self._off_bc(v)
        r1, v1 = self._off_bc(v + 1.0)
        if not np.array_equal(r0, r1):
            raise RuntimeError("Dirichlet row sets differ")
        return r0, v1 - v0

    def _tangent(self, X, rows, dvals):
        jac = self.mod.jacobian_bc(X, rows)
        rs = np.asarray(abs(jac).max(axis=1).todense()).ravel()
        rhs = np.zeros(len(X))
        rhs[rows] = dvals
        return spsolve(sp.csc_matrix(sp.diags(1.0 / rs) @ jac), rhs / rs)

    def _off_solve(self, v, xg):
        rows, vals = self._off_bc(v)
        s = self.mod.solve_bc(xg, _AffineBC(rows, vals), max_step=20.0)
        return s if s.report.converged else None

    def bv_scaled(self, theta, tol_v=1e-8, max_iter=40, warm=True, v_step_v=20.0):
        """BV in V_T: root of I_n(V_D) = 1 on the fixed column.  Returns an :class:`OffState`."""
        self.set_theta(theta)
        th = self.theta
        sc = self.sc
        if warm and self.bv_cache is not None:
            try:
                return self._newton_bv(
                    th, self.bv_cache[0], self.bv_cache[1], None, None, tol_v, max_iter, False
                )
            except EvalError:
                self.stats["warm_fallbacks"] += 1
        # cold: march from equilibrium with the tangent predictor
        eq = self._eq()
        rows0, vals0 = self._off_bc(0.0)
        s0 = self.mod.solve_bc(eq.X, _AffineBC(rows0, vals0), max_step=20.0)
        if not s0.report.converged:
            raise EvalError("off-state zero-bias solve did not converge")
        v, x = 0.0, s0.X.copy()
        lo = (0.0, x)
        h = float(sc.volt_s(v_step_v))
        hi = None
        while True:
            v_try = v + h
            rows, dvals = self._dvals(v)
            pred = x + h * self._tangent(x, rows, dvals)
            s = self._off_solve(v_try, pred)
            if s is None:
                h *= 0.5
                if h < 1e-3:
                    raise EvalError(f"off-state march failed at V_D = {v:.1f} V_T")
                continue
            v, x = v_try, s.X.copy()
            i_n, _ = self._cut_i(x)
            if i_n >= 1.0:
                hi = (v, x)
                break
            lo = (v, x)
            if v > float(sc.volt_s(5000.0)):
                raise EvalError("I_n = 1 not reached below 5 kV")
        # polish from the upper (I >= 1) state with the bracket
        return self._newton_bv(th, hi[0], hi[1], lo[0], hi[0], tol_v, max_iter, True)

    def _newton_bv(self, th, v, x, v_lo, v_hi, tol_v, max_iter, cold):
        sc = self.sc
        xg = x
        cap = 0.15
        for it in range(max_iter):
            s = self._off_solve(v, xg)
            if s is None:
                raise EvalError(f"off-state solve failed at V_D = {v:.3f} V_T")
            X = s.X
            i_n, g = self._cut_i(X)
            rows, dvals = self._dvals(v)
            t = self._tangent(X, rows, dvals)
            di_dv = float(g @ t)
            if i_n <= 0.0 or di_dv <= 0.0:
                raise EvalError("ionization integral vanished / non-monotone")
            if i_n < 1.0:
                v_lo = v if v_lo is None else max(v_lo, v)
            else:
                v_hi = v if v_hi is None else min(v_hi, v)
            dv = -np.log(i_n) / (di_dv / i_n)
            dv = float(np.clip(dv, -cap * v, cap * v))
            v_new = v + dv
            if v_lo is not None and v_hi is not None and not (v_lo < v_new < v_hi):
                v_new = 0.5 * (v_lo + v_hi)
                dv = v_new - v
            if abs(dv) < tol_v and abs(i_n - 1.0) < 1e-9:
                bv = v + dv
                self.bv_cache = (bv, X.copy())
                return OffState(th, bv, float(sc.volt_phys(bv)), X.copy(), v, i_n, di_dv,
                                rows, dvals, g, it + 1, cold)  # fmt: skip
            xg = X + dv * t
            v = v_new
        raise EvalError("BV Newton iteration did not converge")

    @_counted
    def bv(self, theta, grad=True, nodes=False, warm=True):
        """Avalanche BV [V] on the fixed column (and dBV/dtheta [V], optionally dBV/dnd)."""
        st = self.bv_scaled(theta, warm=warm)
        out = {"BV": st.bv_V, "v_s": st.v_s, "state": st, "n_iter": st.n_iter, "cold": st.cold,
               "flips": self.flips}  # fmt: skip
        if not grad:
            return out
        mod = self.mod
        # adjoint at the last Newton iterate (distance to the root < tol_v)
        lam = solve_transposed(mod.jacobian_bc(st.X, st.rows), st.g_psi)
        di_dv = float(lam[st.rows] @ st.dvals)  # dI_n/dV_D = lam^T d(vals)/dV_D
        lam[st.rows] = 0.0
        g_nd = self.dF_dnd(-lam, st.X)  # dI_n/dnd = -lam^T dF/dnd
        sc_v = float(self.sc.vt)
        out["g_nd"] = -g_nd / di_dv * sc_v
        out["grad"] = self.dnd_dtheta.T @ out["g_nd"]
        out["di_dv_adjoint"], out["di_dv_tangent"] = di_dv, st.di_dv
        return out


__all__ = [
    "EvalError",
    "OffState",
    "SiCAdjoint",
    "ionization_n_grad",
    "richardson_fd",
]
