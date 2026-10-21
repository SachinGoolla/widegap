"""Unit tests for stage-6 optimization infrastructure (fast checks, not the gate)."""

import numpy as np

from widegap.fv import reverse_bias as rb
from widegap.optimize import PTDiodeDesign, smoothed_index


def test_parametrization_keeps_resistance_exact():
    des = PTDiodeDesign(h=0.1)
    for seed in range(3):
        th = np.random.default_rng(seed).standard_normal(des.n_modes) * 0.5
        assert abs(des.resistance(th) / des.R0 - 1.0) < 1e-12


def test_db_dtheta_matches_fd():
    des = PTDiodeDesign(h=0.1)
    th = np.array([0.3, -0.2, 0.1, 0.0, 0.05])
    xs = np.linspace(0.0, des.W, 7)
    _, db = des.b(th, xs)
    e = 1e-6
    for k in range(des.n_modes):
        u = np.zeros(des.n_modes)
        u[k] = e
        fd = (des.b(th + u, xs)[0] - des.b(th - u, xs)[0]) / (2 * e)
        assert np.allclose(db[k], fd, rtol=1e-7, atol=1e-9)


def test_smoothed_index_bounds_max():
    des = PTDiodeDesign(h=0.1)
    c, _ = des.c_cv(np.zeros(des.n_modes))
    d = rb.Diode1D(des.x, c, des.delta)
    psi, _ = d.solve(10.0, psi0=d.solve(5.0)[0])
    bs, _ = smoothed_index(d, psi, 15.0, 200.0)
    b = d.max_field(psi) / 15.0
    assert b <= bs <= b + np.log(len(des.x) - 1) / 200.0
