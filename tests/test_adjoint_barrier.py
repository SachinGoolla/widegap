"""Stage 24: nested soft-min barrier and its fixed-voltage off-state adjoint (tiny device)."""

import numpy as np
import pytest

from widegap.adjoint_barrier import BarrierAdjoint, SoftBarrier
from widegap.adjoint_sic import richardson_fd
from widegap.adjoint_sic_oxide import OxideAdjoint
from widegap.bodyredesign import barrier_min, body_design, column_index_matrix, source_columns
from widegap.sic import DMOSDesign, SiC4H, TracedRegistry

pytestmark = pytest.mark.filterwarnings(
    "ignore::RuntimeWarning", "ignore::scipy.sparse.linalg.MatrixRankWarning"
)

BASE = DMOSDesign(
    x_jfet=0.8, ch_len=0.8, src_w=1.0, bc_w=0.8, pad_w=0.4,
    d_src=0.3, d_body=0.8, w_drift=3.0, t_drain=0.3, t_ox=0.05,
    n_drift=5e16, n_jfet=1e17, n_body=3e17, n_src=1e19, n_drain=1e19,
    mesh_kw=dict(nx=(2, 1, 2, 2, 2, 1), ny_drain=1, ny_drift=5, ny_body=4, ny_surf=5, ny_ox=1,
                 q_surf=0.6),
)  # fmt: skip
TH = np.array([0.1, -0.05, 0.03, 0.02, 0.1])
V_TEST = 150.0


def test_softmin_bounds_and_psi_gradient():
    rng = np.random.default_rng(0)
    idx = np.arange(24).reshape(3, 8)
    psi = rng.normal(size=24) * 8.0
    sb = SoftBarrier(idx, 0.5)
    b, g, info = sb.value_grad(psi)
    hard = barrier_min(psi, idx)[0]
    assert info["hard"] == pytest.approx(hard)
    # nested soft-minima: within (ln n_y + ln n_c) tau of the hard barrier
    slack = 0.5 * (np.log(8) + np.log(3))
    assert hard - slack - 1e-9 <= b <= hard + slack + 1e-9
    fd = np.zeros_like(psi)
    for i in range(24):
        e = np.zeros(24)
        e[i] = 1e-6
        fd[i] = (sb.value(psi + e) - sb.value(psi - e)) / 2e-6
    assert np.allclose(g, fd, atol=1e-8)
    # tau -> 0 recovers the hard barrier
    assert SoftBarrier(idx, 1e-3).value(psi) == pytest.approx(hard, abs=0.02)


def test_adjoint_matches_fd_on_tiny_device():
    sc = SiC4H(TracedRegistry.from_json())
    ev = OxideAdjoint(sc, body_design(BASE, 1.2), 0, v_blk_v=V_TEST)
    ev.set_theta(np.zeros(5))
    idx = column_index_matrix(ev, source_columns(ev))
    ba = BarrierAdjoint(ev, idx, v_star_v=V_TEST, tau=0.5)
    o = ba.evaluate(TH)
    assert o["B"] > 0.0 and np.linalg.norm(o["grad"]) > 0.0
    fd = np.zeros(5)
    for k in range(5):
        u = np.zeros(5)
        u[k] = 1.0
        fd[k] = richardson_fd(lambda t: ba.evaluate(t, grad=False)["B"], TH, u, 1e-3)[0]
    assert np.linalg.norm(o["grad"] - fd) <= 1e-6 * np.linalg.norm(fd)
