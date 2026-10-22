"""Exact 1D MOS-capacitor reference (stage 8.1, spec section 2).

Conventions (as in the finite-volume code, scaled units, eps_s = 1):

    -(eps psi')' = rho,  rho = p - n + C,  n = delta e^{psi},  p = delta e^{-psi}
    (equilibrium phi_n = phi_p = 0), complete ionisation, p-type C = -a.

Neutral bulk (rho = 0, psi = psi_b):  p0 = delta e^{-psi_b},  n0 = delta e^{psi_b},
p0 - n0 = a, p0 n0 = delta^2, i.e. p0 = a/2 + sqrt(a^2/4 + delta^2), n0 = delta^2/p0,
psi_b = -ln(p0/delta) = -asinh(a/(2 delta)).

Semiconductor charge.  With u = psi - psi_b (u = psi_s at the surface) the charge
density is rho(u) = p0 (e^{-u} - 1) - n0 (e^{u} - 1).  Let y point from the bulk
towards the surface (psi'(bulk) = 0).  From psi'' = -rho(u), multiply by psi' and
integrate from the bulk:

    (1/2) psi'^2 = -int_0^u rho du' = p0 (e^{-u} + u - 1) + n0 (e^{u} - u - 1) =: W(u) >= 0,

so psi'(surface) = sign(psi_s) sqrt(2 W(psi_s)).  The total semiconductor charge per
area is Q_s = int rho dy = -int psi'' dy = -psi'(surface):

    Q_s(psi_s) = -sign(psi_s) sqrt(2 [p0 (e^{-psi_s} + psi_s - 1) + n0 (e^{psi_s} - psi_s - 1)]).

(Positive psi_s: depletion/inversion, negative charge; negative psi_s: accumulation,
positive charge.)  This is the spec form with a ~ p0 and delta e^{psi_b} = n0.

Oxide.  No charge in the oxide (thickness t_ox, eps_ox = eps_r eps_s) and no interface
charge: D is continuous, eps_s psi'(0-) = eps_ox (psi_G - psi_surf)/t_ox with psi_G the
gate potential (psi_G = V_G - V_FB) and psi_surf = psi_b + psi_s.  Using
psi'(surface) = -Q_s from above:

    V_G - V_FB - psi_b = psi_s - (t_ox / eps_r) Q_s(psi_s)        (eps_s = 1).

The right side is strictly increasing in psi_s, so psi_s(V_G) is found by bracketing
(Brent), independent of the finite-volume code.  Overflow-safe: the n0 terms are
evaluated as exp(ln n0 + u), g(t) = e^t - t - 1 uses a Taylor series for |t| < 0.5 and
expm1 otherwise (so delta = 1e-24 works).
"""

from __future__ import annotations

import math

import numpy as np
from scipy.optimize import brentq


def bulk(a: float, delta: float):
    """(psi_b, p0, n0) of the neutral p-type bulk."""
    p0 = a / 2.0 + math.sqrt(a * a / 4.0 + delta * delta)
    n0 = delta * delta / p0
    return -math.log(p0 / delta), p0, n0


def _g(t: float) -> float:
    """e^t - t - 1 without cancellation."""
    if abs(t) < 0.5:
        term, s = t * t / 2.0, 0.0
        for k in range(3, 40):
            s += term
            term *= t / k
            if abs(term) < 1e-18 * abs(s):
                break
        return s + term
    return math.expm1(t) - t


def W(psi_s: float, a: float, delta: float) -> float:
    """First integral (1/2) psi'^2 at the surface."""
    _, p0, n0 = bulk(a, delta)
    u = psi_s
    pterm = p0 * _g(-u)  # e^{-u} + u - 1 = g(-u)
    if u > 30.0:
        nterm = math.exp(math.log(n0) + u) - n0 * (u + 1.0)
    else:
        nterm = n0 * _g(u)
    return pterm + nterm


def Qs(psi_s: float, a: float, delta: float) -> float:
    if psi_s == 0.0:
        return 0.0
    return -math.copysign(math.sqrt(2.0 * W(psi_s, a, delta)), psi_s)


def gate_voltage(psi_s: float, a: float, delta: float, t_ox: float, eps_r: float, vfb=0.0):
    """V_G for a given surface potential (Gauss across the oxide)."""
    psi_b = bulk(a, delta)[0]
    return vfb + psi_b + psi_s - (t_ox / eps_r) * Qs(psi_s, a, delta)


def psi_s_of_vg(vg: float, a: float, delta: float, t_ox: float, eps_r: float, vfb=0.0):
    """Surface potential (from the bulk) for a gate voltage, by bracketing."""
    f = lambda u: gate_voltage(u, a, delta, t_ox, eps_r, vfb) - vg  # noqa: E731
    lo, hi = -1.0, 1.0
    while f(lo) > 0:
        lo *= 2.0
    while f(hi) < 0:
        hi *= 2.0
    return brentq(f, lo, hi, xtol=1e-14, rtol=4 * np.finfo(float).eps, maxiter=500)
