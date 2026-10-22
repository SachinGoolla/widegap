"""Unit tests for the quasi-Fermi drift-diffusion solver (stage 7.1)."""

from __future__ import annotations

import numpy as np

from widegap.fv.drift_diffusion import DDParams, DriftDiffusion1D
from widegap.fv.equilibrium import cv_average_piecewise, neutral_potential
from widegap.fv.equilibrium import solve as solve_eq
from widegap.fv.qf_drift_diffusion import (
    QFDriftDiffusion1D,
    exact_sg_identity_currents,
    pack,
)


def _diode(delta, a, b, lo, hi, n_el, **prm):
    x = np.linspace(lo, hi, n_el + 1)
    c = cv_average_piecewise(x, [0.0], [-a, b])
    pl, pr = float(neutral_potential(-a, delta)), float(neutral_potential(b, delta))
    eq = solve_eq(x, c, delta, pl, pr, max_iter=400)
    assert eq.report.converged
    params = DDParams(delta=delta, **prm)
    return x, c, pl, pr, eq, params


def test_jacobian_matches_central_fd():
    rng = np.random.default_rng(3)
    for delta, rec in ((1e-3, True), (1e-24, True), (1e-24, False)):
        x = np.linspace(-3.0, 3.0, 25)
        c = cv_average_piecewise(x, [0.0], [-2.0, 1.0])
        prm = DDParams(delta=delta, mu_n=1.3, mu_p=0.7, tau_n=2.0, tau_p=3.0, recombination=rec)
        s = {
            k: (0.1 if delta > 1e-6 else 0.0) * rng.standard_normal(len(x))
            for k in ("psi", "n", "p")
        }
        dd = QFDriftDiffusion1D(x, c, prm, np.log(delta), s)
        X = pack(
            rng.uniform(-3, 3, len(x)),
            rng.uniform(-1, 1, len(x)),
            rng.uniform(-1, 1, len(x)),
        )
        bc = [tuple(X[:3]), tuple(X[-3:])]
        J = dd.jacobian(X)
        for _ in range(3):
            v = rng.standard_normal(len(X))
            v[[0, 1, 2, -3, -2, -1]] = 0.0
            e = 1e-6
            fd = (dd.residual(X + e * v, bc) - dd.residual(X - e * v, bc)) / (2 * e)
            jv = J @ v
            den = abs(J) @ np.abs(v)
            mask = den > 0
            assert np.max(np.abs(jv - fd)[mask] / den[mask]) < 1e-6


def test_equivalent_to_stage3_at_moderate_delta():
    delta = 1e-3
    x, c, pl, pr, eq, params = _diode(delta, 1.0, 1.0, -10.0, 10.0, 256, tau_n=2.0, tau_p=3.0)
    dd3 = DriftDiffusion1D(x, c, params)
    bc3 = [
        (pl + 2.0, delta * np.exp(pl), delta * np.exp(-pl)),
        (pr, delta * np.exp(pr), delta * np.exp(-pr)),
    ]
    X = pack(eq.psi, eq.n, eq.p)
    for v in (0.0, 1.0, 2.0):
        b = [
            (pl + v, delta * np.exp(pl), delta * np.exp(-pl)),
            (pr, delta * np.exp(pr), delta * np.exp(-pr)),
        ]
        sol3 = dd3.solve(X, b, tol_mode="jacobian")
        assert sol3.report.converged
        X = pack(sol3.psi, sol3.n, sol3.p)
    del bc3
    ld = np.log(delta)
    psi, n, p = sol3.psi, sol3.n, sol3.p
    dq = QFDriftDiffusion1D(x, c, params, ld)
    Xq = pack(psi, ld + psi - np.log(n), np.log(p) - ld + psi)
    bcq = [tuple(Xq[:3]), tuple(Xq[-3:])]
    F = dq.residual(Xq, bcq)
    S = dq.row_scale(Xq)
    assert np.max(np.abs(F[3:-3]) / S[3:-3]) < 1e-9
    jn, jp = dq.fluxes(Xq)
    jn3, jp3 = dd3.fluxes(psi, n, p)
    assert np.max(np.abs(jn - jn3)) <= 1e-9 * np.max(np.abs(jn3))
    assert np.max(np.abs(jp - jp3)) <= 1e-9 * np.max(np.abs(jp3))


def test_equilibrium_is_stage2_with_zero_quasi_fermi():
    delta = 1e-24
    x, c, pl, pr, eq, params = _diode(delta, 10.0, 1.0, -15.0, 30.0, 360, recombination=False)
    dd = QFDriftDiffusion1D(x, c, params, np.log(delta))
    z = np.zeros_like(x)
    sol = dd.solve(pack(eq.psi, z, z), [(pl, 0.0, 0.0), (pr, 0.0, 0.0)])
    assert sol.report.converged
    assert np.max(np.abs(sol.psi - eq.psi)) <= 1e-9
    assert max(np.max(np.abs(sol.phi_n)), np.max(np.abs(sol.phi_p))) <= 1e-12


def test_exact_identity_matches_solver_currents():
    delta = 1e-24
    x, c, pl, pr, eq, params = _diode(delta, 10.0, 1.0, -15.0, 30.0, 360, recombination=False)
    dd = QFDriftDiffusion1D(x, c, params, np.log(delta))
    z = np.zeros_like(x)
    X = pack(eq.psi, z, z)
    for v in range(0, 31, 3):
        sol = dd.solve(X, [(pl + v, v, v), (pr, 0.0, 0.0)])
        assert sol.report.converged
        X = sol.X
    jn_e, jp_e = exact_sg_identity_currents(
        x, sol.psi, 30.0, 0.0, 30.0, 0.0, np.log(delta), 1.0, 1.0
    )
    assert abs(sol.jn[0] - jn_e) <= 1e-8 * abs(jn_e)
    assert abs(sol.jp[-1] - jp_e) <= 1e-8 * abs(jp_e)
