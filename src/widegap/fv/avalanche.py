"""Impact ionization and avalanche breakdown, 1D (stage 9b).

Spec: docs/benchmarks/stage9b_avalanche.md. Chynoweth coefficients
alpha_{n,p}(E) = a_{n,p} exp(-b_{n,p}/|E|) are *study parameters* here (scaled units:
alpha in 1/L_D, E in V_T/L_D); material values enter later via the provenance registry.

Orientation (as in stage 4): anode (p) at x[0], cathode (n) at x[-1]. In reverse bias
electrons drift towards +x and holes towards -x; with the code's current convention
(jn = mu_n(n' - n psi'), jp = -mu_p(p' + p psi')) both jn and jp are <= 0 on that branch.

Ionization integral on element-wise constant alpha
--------------------------------------------------
With d = alpha_n - alpha_p, A(x) = int_0^x d, and element k of width h_k carrying
constant (alpha_n, alpha_p) (the element field |E_k|):

    I_n = sum_k alpha_n,k e^{-A_k}        h_k (1 - e^{-z_k})/z_k ,
    I_p = sum_k alpha_p,k e^{A_W - A_{k+1}} h_k (e^{ z_k}-1)/z_k ,   z_k = d_k h_k ,

where A_k = sum_{j<k} z_j (A_W = A_N). Each factor is the *exact* element integral
(expm1-based, no cancellation); I_n is the electron-initiated form
int alpha_n exp(-int_0^x d) dx and I_p the hole-initiated form
int alpha_p exp(-int_x^W (alpha_p - alpha_n)) dx.
Because e^{-A_{k+1}} (e^{z_k}-1) = e^{-A_k}(1 - e^{-z_k}), every element term of I_p equals
e^{A_W} alpha_p,k/alpha_n,k times the corresponding term of I_n, and the telescoping sum
sum_k d_k h_k phi_k e^{-A_k} = 1 - e^{-A_W} holds *exactly in the discrete scheme*. Hence

    1 - I_n = e^{-A_W} (1 - I_p)      (discrete identity),

so I_n = 1 <=> I_p = 1 for any element-wise field: both forms have the same root
(B2). The identity would fail for a quadrature that treats alpha by node values or by a
midpoint rule on the exponent, which is why the exact exponential element factor is used.

DD generation (QF solver, stage 7.1)
-------------------------------------
G = alpha_n|jn| + alpha_p|jp| lives on *elements*: alpha_e = alpha(|E_e|) with the
element field E_e = -(psi_b - psi_a)/h_e, and jn_e, jp_e are the element (edge)
Scharfetter-Gummel fluxes, so g_e = alpha_n,e |jn_e| + alpha_p,e |jp_e| is element
constant. The box method integrates the continuity equation over the control volume
[x_{i-1/2}, x_{i+1/2}], which contains the half of each adjacent element, hence

    vol_i G_i = sum_{e adjacent to i} (h_e/2) g_e

(vol_i = sum of the same half-widths). Continuity: jn' = R - G, jp' = -R + G, i.e.
sources s_n = -G, s_p = +G in the stage-7.1 convention. The total flux jn + jp is
conserved exactly element to element. |j| guard: |j| = abs(j), the derivative uses
sign(j) with sign(0) := -1 (the reverse-bias orientation); ``sign_violation`` reports
the largest wrong-signed (positive) flux relative to the largest flux.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
from scipy.optimize import brentq
from scipy.sparse.linalg import spsolve

from widegap.fv.drift_diffusion import DDParams, bernoulli_prime
from widegap.fv.qf_drift_diffusion import QFDriftDiffusion1D, QFReport, QFSolution, pack
from widegap.fv.reverse_bias import Diode1D, RBReport

EPS_MACH = np.finfo(float).eps


@dataclass(frozen=True)
class Chynoweth:
    a_n: float
    b_n: float
    a_p: float
    b_p: float

    @staticmethod
    def _alpha(e, a, b):
        ae = np.abs(np.asarray(e, dtype=float))
        out = np.zeros_like(ae)
        m = ae > 0.0
        out[m] = a * np.exp(-b / ae[m])
        return out

    @staticmethod
    def _dalpha(e, a, b):
        """d alpha / d|E| = alpha b / |E|^2 (0 where alpha underflows)."""
        ae = np.abs(np.asarray(e, dtype=float))
        al = Chynoweth._alpha(ae, a, b)
        out = np.zeros_like(ae)
        m = al > 0.0
        out[m] = al[m] * b / ae[m] ** 2
        return out

    def alpha_n(self, e):
        return self._alpha(e, self.a_n, self.b_n)

    def alpha_p(self, e):
        return self._alpha(e, self.a_p, self.b_p)

    def dalpha_n(self, e):
        return self._dalpha(e, self.a_n, self.b_n)

    def dalpha_p(self, e):
        return self._dalpha(e, self.a_p, self.b_p)


def _expm1_over(z):
    """(e^z - 1)/z, series near 0."""
    z = np.asarray(z, dtype=float)
    out = np.empty_like(z)
    s = np.abs(z) < 1e-5
    zs = z[s]
    out[s] = 1.0 + zs / 2.0 + zs**2 / 6.0 + zs**3 / 24.0
    zl = z[~s]
    out[~s] = np.expm1(zl) / zl
    return out


def ionization_integrals(x, e_el, ch: Chynoweth):
    """(I_n, I_p, A_W): discrete electron- and hole-initiated ionization integrals
    on element fields ``e_el`` with exact element-wise exponential quadrature."""
    x = np.asarray(x, dtype=float)
    h = np.diff(x)
    an, ap = ch.alpha_n(e_el), ch.alpha_p(e_el)
    z = (an - ap) * h
    big_a = np.concatenate(([0.0], np.cumsum(z)))
    a_w = big_a[-1]
    i_n = float(np.sum(an * h * _expm1_over(-z) * np.exp(-big_a[:-1])))
    i_p = float(np.sum(ap * h * _expm1_over(z) * np.exp(a_w - big_a[1:])))
    return i_n, i_p, float(a_w)


def uniform_field_integral(alpha_n, alpha_p, width):
    """Closed form I_n for a uniform field: alpha_n (1 - e^{-dW})/d, or alpha W if d = 0."""
    z = (alpha_n - alpha_p) * width
    return float(alpha_n * width * _expm1_over(np.array([-z]))[0])


def uniform_field_critical(ch: Chynoweth, width, e_lo=None, e_hi=1e3):
    """E* solving uniform_field_integral(alpha_n(E), alpha_p(E), W) = 1 (Brent)."""
    e_lo = 1e-3 * min(ch.b_n, ch.b_p) if e_lo is None else e_lo

    def f(e):
        return (
            uniform_field_integral(
                ch.alpha_n(np.array([e]))[0], ch.alpha_p(np.array([e]))[0], width
            )
            - 1.0
        )

    return float(brentq(f, e_lo, e_hi, xtol=1e-14, rtol=4 * EPS_MACH))


@dataclass
class AvalancheBV:
    bv: float
    psi: np.ndarray
    integral: float
    form: str
    bracket: tuple[float, float]
    reports: list[RBReport] = field(default_factory=list)


def avalanche_bv(
    diode: Diode1D, ch: Chynoweth, form="n", v_step=2.0, v_max=200.0, xtol=1e-12, level=1.0
) -> AvalancheBV:
    """BV_av: root of I_form(V_R) = ``level`` (1 = breakdown) on the stage-4 zero-current field.

    Bias continuation in steps of ``v_step`` until I >= level, then Brent on the bracket
    with warm starts from the lower bracket end. ``level = 1 - 1/M`` locates the bias
    of a target multiplication M.
    """
    if form not in ("n", "p"):
        raise ValueError("form must be 'n' or 'p'")
    k = 0 if form == "n" else 1
    reports: list[RBReport] = []

    def integral(psi):
        return ionization_integrals(diode.x, diode.element_field(psi), ch)[k]

    psi_lo, v_lo = None, 0.0
    psi_lo, rep = diode.solve(0.0)
    reports.append(rep)
    if integral(psi_lo) >= level:
        raise ValueError("ionization integral already >= level at zero bias")
    v_hi = v_lo
    while True:
        v_hi = v_lo + v_step
        if v_hi > v_max:
            raise ValueError("I = level not bracketed below v_max")
        psi, rep = diode.solve(v_hi, psi0=psi_lo)
        reports.append(rep)
        if integral(psi) >= level:
            break
        psi_lo, v_lo = psi, v_hi

    def g(v):
        p, r = diode.solve(v, psi0=psi_lo)
        reports.append(r)
        return integral(p) - level

    bv = float(brentq(g, v_lo, v_hi, xtol=xtol, rtol=4 * EPS_MACH))
    psi_bv, rep = diode.solve(bv, psi0=psi_lo)
    reports.append(rep)
    return AvalancheBV(bv, psi_bv, integral(psi_bv), form, (v_lo, v_hi), reports)


class AvalancheQFDD1D(QFDriftDiffusion1D):
    """Stage-7.1 QF drift-diffusion with impact-ionization generation (see module doc).

    ``gen_scale`` multiplies G (0 reproduces the base solver; used for continuation).
    """

    def __init__(
        self, x, c_cv, params: DDParams, log_delta, chynoweth, sources=None, gen_scale=1.0
    ):
        super().__init__(x, c_cv, params, log_delta, sources)
        self.ch = chynoweth
        self.gen_scale = float(gen_scale)

    # ----- element quantities ----------------------------------------------------
    def _jn_jp_partials(self, X):
        """Partial derivatives of the element fluxes (same formulas as the base Jacobian)."""
        mu_n, mu_p = self.prm.mu_n, self.prm.mu_p
        D, g, B, Bm, A, E, P, Fp = self._flux_parts(X)
        Bp_, Bmp = bernoulli_prime(D), bernoulli_prime(-D)
        jn = (
            mu_n * g * (Bp_ + B) * A * E,  # psi_b
            -mu_n * g * Bp_ * A * E,  # psi_a
            -mu_n * g * B * A * (E + 1.0),  # phin_b
            mu_n * g * B * A,  # phin_a
        )
        jp = (
            mu_p * g * (Bmp + Bm) * P * Fp,  # psi_b
            -mu_p * g * Bmp * P * Fp,  # psi_a
            -mu_p * g * Bm * P * (Fp + 1.0),  # phip_b
            mu_p * g * Bm * P,  # phip_a
        )
        return jn, jp

    def element_generation(self, X):
        """Element-constant generation g_e and the ingredients for its derivative."""
        psi = X[0::3]
        eel = -np.diff(psi) / self.h
        an, ap = self.ch.alpha_n(eel), self.ch.alpha_p(eel)
        jn, jp = self.fluxes(X)
        g_e = self.gen_scale * (an * np.abs(jn) + ap * np.abs(jp))
        return g_e, eel, an, ap, jn, jp

    def node_generation(self, X):
        """vol_i G_i = sum over adjacent elements of (h_e/2) g_e."""
        g_e = self.element_generation(X)[0]
        out = np.zeros(len(self.x))
        half = 0.5 * self.h * g_e
        out[:-1] += half
        out[1:] += half
        return out

    def sign_violation(self, X):
        """Largest wrong-signed (positive) flux over the largest flux magnitude."""
        jn, jp = self.fluxes(X)
        top = max(float(np.max(np.abs(jn))), float(np.max(np.abs(jp))), 1e-300)
        return float(max(np.max(jn), np.max(jp), 0.0)) / top

    # ----- residual / Jacobian -------------------------------------------------------
    def generation_residual(self, X):
        """The generation part of the residual alone (zero on the Dirichlet rows)."""
        gv = self.node_generation(X)
        gv[0] = gv[-1] = 0.0
        F = np.zeros_like(X)
        F[1::3] += gv
        F[2::3] -= gv
        return F

    def residual(self, X, bc):
        return super().residual(X, bc) + self.generation_residual(X)

    def jacobian(self, X):
        return (super().jacobian(X) + self.generation_jacobian(X)).tocsc()

    def generation_jacobian(self, X):
        """Analytic Jacobian of ``generation_residual``."""
        n_nodes = len(self.x)
        g_e, eel, an, ap, jn, jp = self.element_generation(X)
        del g_e
        s = self.gen_scale
        sn = np.where(jn > 0.0, 1.0, -1.0)
        sp_ = np.where(jp > 0.0, 1.0, -1.0)
        se = np.where(eel >= 0.0, 1.0, -1.0)
        dan, dap = self.ch.dalpha_n(eel), self.ch.dalpha_p(eel)
        (jn_pb, jn_pa, jn_nb, jn_na), (jp_pb, jp_pa, jp_pb_, jp_pa_) = self._jn_jp_partials(X)
        de = (dan * np.abs(jn) + dap * np.abs(jp)) * se / self.h  # dg/d|E| * sign(E) / h
        # dg/dX (columns): psi_b, psi_a, phin_b, phin_a, phip_b, phip_a
        dg = {
            (1, 0): s * (-de + an * sn * jn_pb + ap * sp_ * jp_pb),
            (0, 0): s * (de + an * sn * jn_pa + ap * sp_ * jp_pa),
            (1, 1): s * an * sn * jn_nb,
            (0, 1): s * an * sn * jn_na,
            (1, 2): s * ap * sp_ * jp_pb_,
            (0, 2): s * ap * sp_ * jp_pa_,
        }
        a = np.arange(n_nodes - 1)
        b = a + 1
        node = {0: a, 1: b}
        half = 0.5 * self.h
        rows, cols, vals = [], [], []
        for r_side in (0, 1):  # row node: a or b
            rn = node[r_side]
            keep = (rn != 0) & (rn != n_nodes - 1)
            for r_var, sgn in ((1, 1.0), (2, -1.0)):
                for (c_side, c_var), dv in dg.items():
                    rows.append(3 * rn[keep] + r_var)
                    cols.append(3 * node[c_side][keep] + c_var)
                    vals.append(sgn * half[keep] * dv[keep])
        return sp.csc_matrix(
            (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
            shape=(len(X), len(X)),
        )

    def _volume_terms(self, X):
        out = super()._volume_terms(X)
        gv = self.node_generation(X)
        out[1::3] += gv
        out[2::3] += gv
        return out

    def solve(self, X0, bc, max_iter=100, max_step=2.0):
        """Newton as in the base class, but the step limit acts on *all* unknowns.

        The base solver limits only the psi update; with impact ionization and
        injected minority currents the quasi-Fermi updates can overshoot into
        exp overflow, so the whole update vector is damped by max_step.
        """
        X = np.array(X0, dtype=float)
        n_nodes = len(X) // 3
        X[:3] = bc[0]
        X[3 * n_nodes - 3 :] = bc[-1]
        rep = QFReport()
        for _ in range(max_iter + 1):
            F = self.residual(X, bc)
            J = self.jacobian(X)
            S = self.row_scale(X, J)
            ratio = np.abs(F[3:-3]) / S[3:-3]
            rep.residual_history.append(float(np.max(ratio)))
            if np.all(ratio <= 1e2 * EPS_MACH) and rep.last_dX <= 1e-8:
                rep.converged = True
                break
            if rep.iterations == max_iter:
                break
            rs = np.asarray(abs(J).max(axis=1).todense()).ravel()
            Js = sp.diags(1.0 / rs) @ J
            dX = spsolve(sp.csc_matrix(Js), -F / rs)
            step = float(np.max(np.abs(dX)))
            if step > max_step:
                dX *= max_step / step
                rep.limited_steps += 1
            rep.last_dX = float(np.max(np.abs(dX)))
            X += dX
            rep.iterations += 1
        psi, pn, pp = self.split(X)
        n, p = self.densities(X)
        jn, jp = self.fluxes(X)
        return QFSolution(self.x, X.copy(), psi.copy(), pn.copy(), pp.copy(), n, p, jn, jp, rep)

    def multiplication_integral(self, sol: QFSolution):
        """(I_n, I_p, A_W) on the computed element field of a DD solution."""
        return ionization_integrals(self.x, -np.diff(sol.psi) / self.h, self.ch)


def injected_reverse_solutions(
    diode: Diode1D, ch: Chynoweth, v_target, n_inj=0.01, v0=24.0, dv=3.0
):
    """Reverse-bias DD with minority-electron injection, without and with impact ionization.

    R = 0 (no recombination), so the only currents are the electron current injected at the
    anode contact (n = n_inj there, electron diffusion through the neutral p side) and the
    O(delta^2) hole saturation current. Bias continuation from ``v0`` to ``v_target`` in
    steps ``dv`` (G off), then gen_scale 0 -> 0.5 -> 0.9 -> 1. Returns
    (model_off, sol_off, model_on, sol_on); every Newton report is checked for convergence.
    """
    x, c, delta = diode.x, diode.c_cv, diode.delta
    ld = float(np.log(delta))
    prm = DDParams(delta=delta, recombination=False)
    x_j = float(x[int(np.argmax(c > 0.0))])  # first n-type control volume
    lp = x_j - x[0]

    def bcs(v):
        pa = diode.psi_anode_eq - v
        return [(pa, pa + ld - np.log(n_inj), -v), (diode.psi_cathode_eq, 0.0, 0.0)]

    def guess(v, psi):
        e = np.abs(diode.element_field(psi))
        xm = 0.5 * (x[1:] + x[:-1])
        x_e = xm[(xm > x_j) & (e < 1e-2 * e.max())][0]  # depletion edge on the n side
        nlin = np.where(x < x_j, n_inj * np.clip((x_j - x) / lp, 0.0, 1.0), 0.0)
        phin = np.minimum(0.0, psi + ld - np.log(np.maximum(nlin, delta * np.exp(psi))))
        t = np.clip((x - x_e) / 20.0, 0.0, 1.0)
        return pack(psi, phin, -v * (1.0 - t * t * (3.0 - 2.0 * t)))

    v0 = min(float(v0), float(v_target))
    psi_old = None
    for v in np.arange(0.0, v0 + 1e-9, 2.0):
        psi_old, _ = diode.solve(v, psi0=psi_old)
    if abs(v0 - np.floor(v0 / 2.0) * 2.0) > 1e-9:
        psi_old, _ = diode.solve(v0, psi0=psi_old)
    m0 = AvalancheQFDD1D(x, c, prm, ld, ch, gen_scale=0.0)
    sol = m0.solve(guess(v0, psi_old), bcs(v0))
    if not sol.report.converged:
        raise RuntimeError(f"DD baseline did not converge at V={v0}")
    v = v0
    while v < v_target - 1e-12:
        vn = min(v + dv, float(v_target))
        pn, _ = diode.solve(vn, psi0=psi_old)
        xo = sol.X
        xn = xo.copy()
        xn[0::3] = pn
        xn[1::3] = xo[1::3] + np.where(x < x_j, pn - psi_old, 0.0)
        xn[2::3] = xo[2::3] * (vn / v)
        sol = m0.solve(xn, bcs(vn))
        if not sol.report.converged:
            raise RuntimeError(f"DD baseline did not converge at V={vn}")
        v, psi_old = vn, pn
    sol_off = sol
    X = sol.X
    m1, sol_on = m0, sol
    for gs in (0.5, 0.9, 1.0):
        m1 = AvalancheQFDD1D(x, c, prm, ld, ch, gen_scale=gs)
        sol_on = m1.solve(X, bcs(v_target))
        if not sol_on.report.converged:
            raise RuntimeError(
                f"DD with generation did not converge at V={v_target}, gen_scale={gs}"
            )
        X = sol_on.X.copy()
    return m0, sol_off, m1, sol_on


__all__ = [
    "AvalancheBV",
    "AvalancheQFDD1D",
    "Chynoweth",
    "avalanche_bv",
    "injected_reverse_solutions",
    "ionization_integrals",
    "uniform_field_critical",
    "uniform_field_integral",
]
