"""Unit tests for stage 9b (impact ionization, avalanche breakdown, 1D)."""

from __future__ import annotations

import numpy as np

from widegap.fv.avalanche import (
    AvalancheQFDD1D,
    Chynoweth,
    avalanche_bv,
    injected_reverse_solutions,
    ionization_integrals,
    uniform_field_critical,
    uniform_field_integral,
)
from widegap.fv.drift_diffusion import DDParams
from widegap.fv.equilibrium import cv_average_piecewise
from widegap.fv.qf_drift_diffusion import QFDriftDiffusion1D, pack
from widegap.fv.reverse_bias import Diode1D

CH = Chynoweth(6.0, 30.0, 3.0, 39.0)


def _pn(h=0.2, delta=1e-3):
    x = np.linspace(-5.0, 60.0, int(round(65.0 / h)) + 1)
    return Diode1D(x, cv_average_piecewise(x, [0.0], [-100.0, 1.0]), delta)


def test_chynoweth_basics():
    e = np.array([0.0, 1e-3, 1.0, 5.0, -5.0, 50.0])
    a = CH.alpha_n(e)
    assert a[0] == 0.0 and np.all(np.diff(a[[1, 2, 3, 5]]) > 0) and a[3] == a[4]
    assert np.isclose(a[3], 6.0 * np.exp(-6.0))
    eps = 1e-6
    fd = (CH.alpha_p(np.array([4.0 + eps])) - CH.alpha_p(np.array([4.0 - eps]))) / (2 * eps)
    assert np.isclose(CH.dalpha_p(np.array([4.0]))[0], fd[0], rtol=1e-8)


def test_uniform_field_closed_form_and_quadrature():
    rng = np.random.default_rng(0)
    for an, ap in ((0.3, 0.1), (0.1, 0.3), (0.2, 0.2), (0.2, 0.2 + 1e-9)):
        x = np.cumsum(np.r_[0.0, rng.uniform(0.05, 0.4, 40)])
        w = x[-1]
        # a constant field with element-wise alphas (use a Chynoweth with equal exp factor)
        ch = Chynoweth(an, 0.0, ap, 0.0)  # b = 0: alpha = a for E != 0
        i_n, i_p, _ = ionization_integrals(x, np.ones(len(x) - 1), ch)
        ref = uniform_field_integral(an, ap, w)
        assert abs(i_n - ref) <= 1e-12 * ref
        assert abs(i_p - uniform_field_integral(ap, an, w)) <= 1e-12 * ref
    assert np.isclose(uniform_field_integral(0.25, 0.25, 4.0), 1.0, rtol=1e-14)
    # critical field: round trip, including alpha_n == alpha_p (alpha W = 1)
    for ch, w in ((CH, 3.0), (Chynoweth(1.0, 5.0, 1.0, 5.0), 2.0)):
        e = uniform_field_critical(ch, w)
        i = uniform_field_integral(ch.alpha_n(np.array([e]))[0], ch.alpha_p(np.array([e]))[0], w)
        assert abs(i - 1.0) < 1e-12
    e = uniform_field_critical(Chynoweth(1.0, 5.0, 1.0, 5.0), 2.0)
    assert np.isclose(np.exp(-5.0 / e) * 2.0, 1.0, rtol=1e-12)


def test_discrete_identity_electron_vs_hole_form():
    rng = np.random.default_rng(1)
    x = np.cumsum(np.r_[0.0, rng.uniform(0.02, 0.3, 120)])
    for scale in (0.2, 1.0, 3.0):
        e = scale * rng.uniform(0.0, 12.0, len(x) - 1) * (rng.random(len(x) - 1) > 0.2)
        i_n, i_p, a_w = ionization_integrals(x, e, CH)
        lhs, rhs = 1.0 - i_n, np.exp(-a_w) * (1.0 - i_p)
        assert abs(lhs - rhs) <= 1e-12 * max(1.0, abs(lhs), abs(rhs))


def test_both_forms_same_bv():
    d = _pn(0.2)
    rn = avalanche_bv(d, CH, "n")
    rp = avalanche_bv(d, CH, "p")
    assert 20.0 < rn.bv < 150.0
    assert abs(rn.bv - rp.bv) <= 1e-9 * rn.bv
    assert abs(rn.integral - 1.0) < 1e-10 and abs(rp.integral - 1.0) < 1e-10


def _random_model(delta, rec, strong=True):
    rng = np.random.default_rng(7)
    x = np.linspace(-3.0, 3.0, 25)
    c = cv_average_piecewise(x, [0.0], [-2.0, 1.0])
    prm = DDParams(delta=delta, mu_n=1.3, mu_p=0.7, tau_n=2.0, tau_p=3.0, recombination=rec)
    s = {k: (0.1 if delta > 1e-6 else 0.0) * rng.standard_normal(len(x)) for k in ("psi", "n", "p")}
    dd = AvalancheQFDD1D(x, c, prm, np.log(delta), Chynoweth(3.0, 4.0, 1.5, 5.0), s)
    X = pack(rng.uniform(-3, 3, len(x)), rng.uniform(-1, 1, len(x)), rng.uniform(-1, 1, len(x)))
    return dd, X, rng


def test_jacobian_matches_fd_per_row():
    for delta, rec in ((1e-3, True), (1e-24, True), (1e-24, False)):
        dd, X, rng = _random_model(delta, rec)
        bc = [tuple(X[:3]), tuple(X[-3:])]
        for gen_only in (False, True):
            jac = dd.generation_jacobian(X) if gen_only else dd.jacobian(X)

            def fun(y, gen_only=gen_only, dd=dd, bc=bc):
                return dd.generation_residual(y) if gen_only else dd.residual(y, bc)

            for _ in range(3):
                v = rng.standard_normal(len(X))
                v[[0, 1, 2, -3, -2, -1]] = 0.0
                e = 1e-6
                fd = (fun(X + e * v) - fun(X - e * v)) / (2 * e)
                den = abs(jac) @ np.abs(v)
                m = den > 0
                assert np.max(np.abs(jac @ v - fd)[m] / den[m]) < 1e-6


def test_generation_off_reproduces_base_and_is_active():
    dd, X, _ = _random_model(1e-3, True)
    bc = [tuple(X[:3]), tuple(X[-3:])]
    base = QFDriftDiffusion1D(
        dd.x, dd.c, dd.prm, dd.ld, {"psi": dd.s_psi, "n": dd.s_n, "p": dd.s_p}
    )
    off = AvalancheQFDD1D(dd.x, dd.c, dd.prm, dd.ld, dd.ch, gen_scale=0.0)
    off.s_psi, off.s_n, off.s_p = dd.s_psi, dd.s_n, dd.s_p
    assert np.array_equal(off.residual(X, bc), base.residual(X, bc))
    assert abs(off.jacobian(X) - base.jacobian(X)).max() == 0.0
    assert np.max(np.abs(dd.generation_residual(X))) > 0.0
    # node generation partitions the element generation exactly
    g_e = dd.element_generation(X)[0]
    assert np.isclose(dd.node_generation(X).sum(), np.sum(dd.h * g_e), rtol=1e-13)


def test_multiplication_matches_ionization_integral():
    d = _pn(0.2)
    for m_t in (2.0, 8.0):
        rb = avalanche_bv(d, CH, "n", level=1.0 - 1.0 / m_t)
        m0, s0, m1, s1 = injected_reverse_solutions(d, CH, rb.bv)
        i_off, i_on = -(s0.jn + s0.jp), -(s1.jn + s1.jp)
        assert np.all(i_on > 0) and i_on.max() / i_on.min() - 1.0 < 1e-6  # conserved flux
        assert m1.sign_violation(s1.X) < 1e-6
        i_n = m1.multiplication_integral(s1)[0]
        m_dd = i_on[-1] / i_off[-1]
        assert abs(m_dd * (1.0 - i_n) - 1.0) < 1e-2
        assert 0.8 * m_t < m_dd < 1.1 * m_t
