"""Discrete adjoint doping sensitivities (stage 5).

Independent implementation of the adjoint doping-sensitivity method of
Zhu & Andrei, AIP Adv. 9, 095301 (2019), doi:10.1063/1.5113764, sec. III,
on the verified forward solvers of stages 2-4
(docs/benchmarks/stage5_adjoint.md).

For F(X, C) = 0 and objective J(X): solve F_X^T lam = dJ/dX, then
dJ/dC_i = -lam^T dF/dC_i = lam_psi,i * V_i for interior nodes i (C enters the
Poisson rows only, as -V_i C_i). Contact nodes are not design variables.
Sensitivities are with respect to C = (N_D - N_A)/N0; ZA19 use D = N_A - N_D,
so gamma_D = -gamma_C.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import solve_banded

from widegap.fv import equilibrium as eq
from widegap.fv import reverse_bias as rb
from widegap.fv.drift_diffusion import DriftDiffusion1D, bernoulli, bernoulli_prime
from widegap.fv.equilibrium import control_volumes
from widegap.linalg import safe_spsolve


def banded_transpose(ab: np.ndarray) -> np.ndarray:
    """Transpose of a tridiagonal matrix in scipy (1, 1) banded storage."""
    at = np.zeros_like(ab)
    at[1] = ab[1]
    at[0, 1:] = ab[2, :-1]
    at[2, :-1] = ab[0, 1:]
    return at


def _interior(grad: np.ndarray) -> np.ndarray:
    grad = grad.copy()
    grad[0] = grad[-1] = 0.0
    return grad


# ----- A1: stage-2 equilibrium, J = int psi dx ------------------------------
def equilibrium_integral_psi(psi, x, delta):
    """Returns (J, dJ/dC) for J = trapezoid integral of psi."""
    vol = control_volumes(x)  # trapezoid weights equal control volumes
    lam = solve_banded((1, 1), banded_transpose(eq.jacobian_banded(psi, x, delta)), vol)
    return float(np.sum(vol * psi)), _interior(lam * vol)


# ----- A2: stage-3 drift-diffusion, J = mean element current ----------------
def dd_mean_current_grad_x(dd: DriftDiffusion1D, X: np.ndarray) -> np.ndarray:
    """d(mean_e (jn_e + jp_e)) / dX, from the analytic SG flux derivatives."""
    psi, n, p = dd.split(X)
    h, mu_n, mu_p = dd.h, dd.prm.mu_n, dd.prm.mu_p
    dl = np.diff(psi)
    bp, bm = bernoulli(dl), bernoulli(-dl)
    dbp, dbm = bernoulli_prime(dl), bernoulli_prime(-dl)
    g = 1.0 / h
    ne = len(h)
    w = 1.0 / ne
    grad = np.zeros_like(X)
    # electron flux jn = mu_n/h [B(D) n_b - B(-D) n_a]
    np.add.at(grad, 3 * np.arange(1, ne + 1) + 1, w * mu_n * g * bp)
    np.add.at(grad, 3 * np.arange(ne) + 1, -w * mu_n * g * bm)
    djn_dD = mu_n * g * (dbp * n[1:] + dbm * n[:-1])
    # hole flux jp = -mu_p/h [B(-D) p_b - B(D) p_a]
    np.add.at(grad, 3 * np.arange(1, ne + 1) + 2, -w * mu_p * g * bm)
    np.add.at(grad, 3 * np.arange(ne) + 2, w * mu_p * g * bp)
    djp_dD = mu_p * g * (dbm * p[1:] + dbp * p[:-1])
    dD = w * (djn_dD + djp_dD)
    np.add.at(grad, 3 * np.arange(1, ne + 1), dD)
    np.add.at(grad, 3 * np.arange(ne), -dD)
    return grad


def dd_mean_current(dd: DriftDiffusion1D, X: np.ndarray):
    """Returns (J, dJ/dC) for J = mean element total current at fixed bias."""
    psi, n, p = dd.split(X)
    jn, jp = dd.fluxes(psi, n, p)
    lam = safe_spsolve(dd.jacobian(X).T.tocsc(), dd_mean_current_grad_x(dd, X))
    return float(np.mean(jn + jp)), _interior(lam[0::3] * dd.vol)


# ----- A3/A4: stage-4 breakdown index and field-based BV --------------------
def breakdown_index(d: rb.Diode1D, psi, v_r, e_crit):
    """Returns (B, dB/dC, dB/dV_R, e_star) at reverse bias v_r."""
    phi_n, phi_p = 0.0, -float(v_r)
    e = d.element_field(psi)
    k = int(np.argmax(np.abs(e)))
    s = np.sign(e[k])
    h = d.x[k + 1] - d.x[k]
    gx = np.zeros_like(psi)
    gx[k + 1] = -s / (h * e_crit)  # E_k = -(psi_{k+1} - psi_k)/h
    gx[k] = s / (h * e_crit)
    ab = rb.jacobian_banded(psi, d.x, d.delta, phi_n, phi_p)
    lam = solve_banded((1, 1), banded_transpose(ab), gx)
    vol = control_volumes(d.x)
    # dF/dV_R: interior rows via phi_p = -V_R; anode Dirichlet row psi_0 - (psi_eq - V_R)
    f_v = np.zeros_like(psi)
    f_v[1:-1] = vol[1:-1] * d.delta * np.exp(phi_p - psi[1:-1])
    f_v[0] = 1.0
    b = float(abs(e[k]) / e_crit)
    return b, _interior(lam * vol), float(-lam @ f_v), k


def breakdown_voltage_sensitivity(d: rb.Diode1D, psi_bv, bv, e_crit):
    """dBV/dC = -(dB/dC) / (dB/dV_R) at B(BV) = 1 (ZA19 eq. 14)."""
    _, db_dc, db_dv, k = breakdown_index(d, psi_bv, bv, e_crit)
    return -db_dc / db_dv, k


# ----- polishing Newton steps for finite-difference forward solves ----------
def polish_equilibrium(psi, x, c, delta, steps=3):
    s0 = np.zeros_like(psi)
    for _ in range(steps):
        f = eq.residual(psi, x, c, delta, s0)
        psi = psi + solve_banded((1, 1), eq.jacobian_banded(psi, x, delta), -f)
    return psi


def polish_dd(dd: DriftDiffusion1D, X, bc, steps=3):
    for _ in range(steps):
        X = X + safe_spsolve(dd.jacobian(X), -dd.residual(X, bc))
    return X


def polish_reverse(d: rb.Diode1D, psi, v_r, steps=3):
    args = (d.x, d.c_cv, d.delta, 0.0, -float(v_r))
    for _ in range(steps):
        f = rb.residual(psi, *args)
        psi = psi + solve_banded(
            (1, 1), rb.jacobian_banded(psi, d.x, d.delta, 0.0, -float(v_r)), -f
        )
    return psi
