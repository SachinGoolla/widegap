"""Generic 2D discrete adjoint on the stage-8.1 MOS quasi-Fermi drift-diffusion model.

Independent implementation of the adjoint doping-sensitivity method of Zhu & Andrei,
AIP Adv. 9, 095301 (2019), and Brungi & Andrei (2026), on the verified solver of
``widegap.fv2d.mos`` (docs/benchmarks/stage10_mosfet_optimization.md).

For ``F(X, C; p) = 0`` (X = (psi, phi_n, phi_p) interleaved per node, Dirichlet rows
included) and an objective ``J(X)``:

* solve ``F_X^T lam = dJ/dX`` with the BC-modified Jacobian of the forward Newton solver;
* the doping enters the Poisson rows only as ``-vol_s,i C_i`` (``vol_s`` = semiconductor
  part of the control volume), so ``dJ/dC_i = lam_psi,i vol_s,i`` on free nodes;
* the drain-current objective ``I_D`` is the contact current of
  ``QFDriftDiffusion.contact_current`` (sum of ``jn + jp`` over the edges crossing the
  drain node set); its ``dI_D/dX`` follows from the analytic Scharfetter-Gummel flux
  derivatives used in the Jacobian;
* the smoothed field index ``B_s = kappa^-1 log sum_e exp(kappa |E_e| / E_crit)`` over a
  set of semiconductor edges (stage-6 smoothing) has an explicit ``dB_s/dX`` through the
  edge fields ``E_e = (psi_b - psi_a) / h_e``;
* ``BV_s`` solves ``B_s(X(V_D, C)) = 1`` (drain-bias continuation bracketing + Brent);
  by the implicit function theorem ``dBV_s/dC = -(dB_s/dC) / (dB_s/dV_D)`` with
  ``dB_s/dV_D = -lam^T dF/dV_D`` and ``dF/dV_D = -1`` on the three Dirichlet rows of every
  drain node, i.e. ``dB_s/dV_D = sum_{drain rows} lam``.

Sensitivities are with respect to ``C`` (net donor, scaled); ZA19/BA26 use
``D = N_A - N_D`` so ``gamma_D = -gamma_C``.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
from scipy.optimize import brentq
from scipy.sparse.linalg import spsolve

from widegap.fv.drift_diffusion import bernoulli_prime
from widegap.fv2d.mos import MOSDriftDiffusion, continue_bias


# ----- transposed solve --------------------------------------------------------------
def solve_transposed(J: sp.spmatrix, g: np.ndarray, refine: int = 2) -> np.ndarray:
    """Solve ``J^T lam = g`` with the row equilibration of the forward Newton solver
    (J = diag(rs) J_s  =>  J_s^T (rs * lam) = g) and iterative refinement on the raw J."""
    J = sp.csc_matrix(J)
    rs = np.asarray(abs(J).max(axis=1).todense()).ravel()
    Js_t = (sp.diags(1.0 / rs) @ J).T.tocsc()
    lam = spsolve(Js_t, g) / rs
    Jt = J.T.tocsr()
    for _ in range(refine):
        r = g - Jt @ lam
        lam = lam + spsolve(Js_t, r) / rs
    return lam


def adjoint_solve(mod: MOSDriftDiffusion, X: np.ndarray, bc, g: np.ndarray) -> np.ndarray:
    rows, _ = mod._all_bc(bc)
    return solve_transposed(mod.jacobian_bc(X, rows), g)


def dj_dc(mod: MOSDriftDiffusion, lam: np.ndarray, bc) -> np.ndarray:
    """dJ/dC_i = lam_psi,i vol_s,i on free (non-Dirichlet-psi) nodes, 0 on contact nodes."""
    rows, _ = mod._all_bc(bc)
    out = lam[0::3] * mod.vol
    psi_rows = rows[rows % 3 == 0] // 3
    out[psi_rows] = 0.0
    return out


# ----- drain-current objective ----------------------------------------------------------
def flux_partials(mod: MOSDriftDiffusion, X: np.ndarray):
    """Per-edge derivatives of (jn, jp) w.r.t. (psi_a, psi_b, phi_a, phi_b) (same formulas as
    the analytic Jacobian): returns dict with keys 'n' and 'p', each a 4-tuple."""
    D, B, Bm, A, E, P, Fp = mod._flux_parts(X)
    Bp_, Bmp = bernoulli_prime(D), bernoulli_prime(-D)
    g, mu_n, mu_p = mod.gw, mod.prm.mu_n, mod.prm.mu_p
    n = (
        -mu_n * g * Bp_ * A * E,  # psi_a
        mu_n * g * (Bp_ + B) * A * E,  # psi_b
        mu_n * g * B * A,  # phi_n,a
        -mu_n * g * B * A * (E + 1.0),  # phi_n,b
    )
    p = (
        -mu_p * g * Bmp * P * Fp,
        mu_p * g * (Bmp + Bm) * P * Fp,
        mu_p * g * Bm * P,
        -mu_p * g * Bm * P * (Fp + 1.0),
    )
    return {"n": n, "p": p}


def contact_current_grad_x(mod: MOSDriftDiffusion, X: np.ndarray, nodes):
    """(I, dI/dX) for I = ``mod.contact_current`` of the node set (drain-current objective)."""
    N = mod.mesh.n_nodes
    inset = np.zeros(N, dtype=bool)
    inset[np.asarray(nodes, dtype=int)] = True
    a, b = mod.a, mod.b
    sgn = (inset[a] & ~inset[b]).astype(float) - (~inset[a] & inset[b]).astype(float)
    jn, jp = mod.fluxes(X)
    cur = float(np.sum(sgn * (jn + jp)))
    fp = flux_partials(mod, X)
    g = np.zeros(3 * N)

    def acc(node, var, val):
        g[var::3] += np.bincount(node, sgn * val, minlength=N)

    npa, npb, nna, nnb = fp["n"]
    ppa, ppb, pna, pnb = fp["p"]
    acc(a, 0, npa + ppa)
    acc(b, 0, npb + ppb)
    acc(a, 1, nna)
    acc(b, 1, nnb)
    acc(a, 2, pna)
    acc(b, 2, pnb)
    return cur, g


def drain_current_grad_c(mod: MOSDriftDiffusion, X: np.ndarray, bc, nodes):
    """(I_D, dI_D/dC at nodes) at fixed bias by the adjoint."""
    cur, gx = contact_current_grad_x(mod, X, nodes)
    lam = adjoint_solve(mod, X, bc, gx)
    return cur, dj_dc(mod, lam, bc)


# ----- field index ---------------------------------------------------------------------
@dataclass
class FieldIndex:
    """B(X) over the edges flagged in ``mask``.  ``kappa=None`` -> exact max (no gradient)."""

    e_crit: float
    mask: np.ndarray
    kappa: float | None = 200.0

    def edge_fields(self, mod, X):
        psi = X[0::3]
        return (psi[mod.b] - psi[mod.a]) / mod.mesh.h

    def __call__(self, mod: MOSDriftDiffusion, X: np.ndarray, grad: bool = False):
        e = self.edge_fields(mod, X)[self.mask]
        be = np.abs(e) / self.e_crit
        m = float(be.max())
        if self.kappa is None:
            if grad:
                raise ValueError("exact index has no gradient")
            return m
        w = np.exp(self.kappa * (be - m))
        bs = m + np.log(w.sum()) / self.kappa
        if not grad:
            return float(bs)
        sm = w / w.sum()
        idx = np.flatnonzero(self.mask)
        coef = sm * np.sign(e) / (mod.mesh.h[idx] * self.e_crit)
        N = mod.mesh.n_nodes
        g = np.zeros(3 * N)
        g[0::3] += np.bincount(mod.b[idx], coef, minlength=N) - np.bincount(
            mod.a[idx], coef, minlength=N
        )
        return float(bs), g


def semiconductor_edges(mod: MOSDriftDiffusion) -> np.ndarray:
    """Edges carrying semiconductor dual-face weight (excludes pure-oxide edges)."""
    return mod.reg.w_s > 0


# ----- BV_s by drain-bias continuation + Brent -------------------------------------------
@dataclass
class BVResult:
    bv: float
    X: np.ndarray
    bc: object
    cache: list = field(default_factory=list)  # [(V, X)] converged states for warm starts
    solves: int = 0
    newton_its: int = 0


def _guess(cache, v, fallback):
    """Linear interpolation (in V) of cached states; nearest if outside; else fallback."""
    if not cache:
        return fallback
    vs = np.array([c[0] for c in cache])
    k = int(np.searchsorted(vs, v))
    if k == 0:
        return cache[0][1] if abs(vs[0] - v) < 1.0 else fallback
    if k == len(vs):
        return cache[-1][1] if abs(vs[-1] - v) < 1.0 else fallback
    (v0, x0), (v1, x1) = cache[k - 1], cache[k]
    t = (v - v0) / (v1 - v0)
    return x0 + t * (x1 - x0)


def breakdown_voltage(
    mod: MOSDriftDiffusion,
    bc_of,
    index: FieldIndex,
    X0: np.ndarray,
    guesses=(),
    v_step: float = 25.0,
    vmax: float = 1000.0,
    xtol: float = 1e-10,
    target: float = 1.0,
) -> BVResult:
    """Smallest drain bias V with ``index(X(V)) = target`` (bracket by marching in
    ``v_step``, then Brent).  ``bc_of(V)`` gives the Dirichlet data at drain bias V;
    ``X0`` is the converged state at V = 0; ``guesses`` is a list of (V, X) from earlier
    (nearby) designs used as Newton starting points."""
    guesses = sorted(guesses, key=lambda t: t[0])
    res = BVResult(0.0, X0, None)

    def solve(v, x_start, x_fallback):
        out = continue_bias(mod, x_start, lambda p: bc_of(p[0]), [v_from[0]], [[v]])
        s = out[-1][1]
        res.solves += out[-1][3]
        res.newton_its += out[-1][2]
        if not s.report.converged and x_fallback is not None:
            out = continue_bias(mod, x_fallback, lambda p: bc_of(p[0]), [v_from[0]], [[v]])
            s = out[-1][1]
            res.solves += out[-1][3]
            res.newton_its += out[-1][2]
        return s

    v_from = [0.0]
    b0 = index(mod, X0)
    if b0 >= target:
        raise ValueError("index >= target at zero bias")
    cache = [(0.0, np.array(X0))]
    v_lo, x_lo = 0.0, np.array(X0)
    v_hi = x_hi = None
    v = 0.0
    while v < vmax:
        v = min(v + v_step, vmax)
        v_from[0] = v_lo
        s = solve(v, _guess(guesses, v, x_lo), x_lo)
        if not s.report.converged:
            raise RuntimeError(f"forward solve did not converge at V_D={v}")
        b = index(mod, s.X)
        if b >= target:
            v_hi, x_hi = v, s.X.copy()
            break
        cache.append((v, s.X.copy()))
        v_lo, x_lo = v, s.X.copy()
    if v_hi is None:
        raise ValueError("index = target not reached by vmax")

    local = [(v_lo, x_lo), (v_hi, x_hi)]

    def f(vq):
        v_from[0] = v_lo
        xg = _guess(sorted(local, key=lambda t: t[0]), vq, x_lo)
        s = solve(vq, xg, x_lo)
        if not s.report.converged:
            raise RuntimeError(f"forward solve did not converge at V_D={vq}")
        local.append((vq, s.X.copy()))
        return index(mod, s.X) - target

    bv = brentq(f, v_lo, v_hi, xtol=xtol, rtol=4 * np.finfo(float).eps)
    v_from[0] = v_lo
    xg = _guess(sorted(local, key=lambda t: t[0]), bv, x_lo)
    s = solve(bv, xg, x_lo)
    if not s.report.converged:
        raise RuntimeError("final polish solve did not converge")
    res.bv, res.X, res.bc = float(bv), s.X.copy(), bc_of(bv)
    res.cache = cache + [(bv, s.X.copy())]
    return res


def bv_gradient_c(mod: MOSDriftDiffusion, X: np.ndarray, bc, drain_nodes, index: FieldIndex):
    """(dBV_s/dC at nodes, dB_s/dV_D, B_s) at the breakdown state by the adjoint and the
    implicit function theorem."""
    b_s, gx = index(mod, X, grad=True)
    lam = adjoint_solve(mod, X, bc, gx)
    dn = np.asarray(drain_nodes, dtype=int)
    rows = np.r_[3 * dn, 3 * dn + 1, 3 * dn + 2]
    db_dv = float(lam[rows].sum())  # -lam^T dF/dV_D with dF/dV_D = -1 on the drain rows
    db_dc = dj_dc(mod, lam, bc)
    return -db_dc / db_dv, db_dv, b_s
