"""Incomplete ionization and equilibrium carrier statistics (stage 7).

All functions work in the scaled units of docs/equations.md section 1.2
(potential in V_T, densities in N0) unless a docstring says otherwise.

Ionized fractions (Boltzmann carriers, quasi-Fermi level 0):
    N_D^+ = N_D / (1 + g_D (n / N_c) e^{dE_D/kT}) = N_D / (1 + exp(psi + l_D))
    N_A^- = N_A / (1 + g_A (p / N_v) e^{dE_A/kT}) = N_A / (1 + exp(-psi + l_A))
with l_D = log(g_D * delta * e^{dE_D/kT} / nu_c), l_A = log(g_A * delta * e^{dE_A/kT} / nu_v),
delta = n_i / N0, nu_c = N_c / N0, nu_v = N_v / N0. Logs keep the arithmetic
safe for 4H-SiC (delta ~ 1e-25).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

K_B = 1.380649e-23  # J/K, exact (SI 2019)
Q = 1.602176634e-19  # C, exact (SI 2019)


@dataclass(frozen=True)
class Level:
    """One dopant level: fraction of the species on this site and its log-gamma l."""

    fraction: float
    log_gamma: float


def log_gamma_donor(g_d, de_d_ev, t_k, delta, nu_c):
    return float(np.log(g_d) + np.log(delta) + de_d_ev * Q / (K_B * t_k) - np.log(nu_c))


def log_gamma_acceptor(g_a, de_a_ev, t_k, delta, nu_v):
    return float(np.log(g_a) + np.log(delta) + de_a_ev * Q / (K_B * t_k) - np.log(nu_v))


def ionized_donors(psi, nd, levels):
    return sum(lv.fraction * nd / (1.0 + np.exp(psi + lv.log_gamma)) for lv in levels)


def ionized_acceptors(psi, na, levels):
    return sum(lv.fraction * na / (1.0 + np.exp(-psi + lv.log_gamma)) for lv in levels)


def charge(psi, nd, na, delta, donors, acceptors):
    """Scaled space charge rho = p - n + N_D^+ - N_A^-."""
    return (
        delta * (np.exp(-psi) - np.exp(psi))
        + ionized_donors(psi, nd, donors)
        - ionized_acceptors(psi, na, acceptors)
    )


def dcharge_dpsi(psi, nd, na, delta, donors, acceptors):
    out = -delta * (np.exp(-psi) + np.exp(psi))
    for lv in donors:
        s = 1.0 / (1.0 + np.exp(-(psi + lv.log_gamma)))  # e^z/(1+e^z)
        out = out - lv.fraction * nd * s * (1.0 - s)
    for lv in acceptors:
        t = 1.0 / (1.0 + np.exp(-(-psi + lv.log_gamma)))
        out = out - lv.fraction * na * t * (1.0 - t)
    return out


def charge_antiderivative(psi, nd, na, delta, donors, acceptors):
    """Exact antiderivative of ``charge`` in psi (for first-integral references)."""
    out = -delta * (np.exp(-psi) + np.exp(psi))
    for lv in donors:  # int 1/(1+e^{psi+l}) = psi - log(1 + e^{psi+l})
        out = out + lv.fraction * nd * (psi - np.logaddexp(0.0, psi + lv.log_gamma))
    for lv in acceptors:  # int 1/(1+e^{-psi+l}) = log(1 + e^{psi-l})
        out = out - lv.fraction * na * np.logaddexp(0.0, psi - lv.log_gamma)
    return out


def neutral_potential(nd, na, delta, donors, acceptors, lo=-400.0, hi=400.0):
    """Independent bulk charge-neutrality root (rho is strictly decreasing in psi)."""

    def f(p):
        return float(charge(p, nd, na, delta, donors, acceptors))

    return brentq(f, lo, hi, xtol=1e-14, rtol=4 * np.finfo(float).eps, maxiter=500)


def single_donor_electrons(nd, n1):
    """Closed form of n^2 + n n1 - N_D n1 = 0 (cancellation-free); any consistent units."""
    return 2.0 * nd * n1 / (n1 + np.sqrt(n1 * n1 + 4.0 * nd * n1))
