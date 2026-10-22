"""Stage 11 study helpers: fast checks with tiny configurations."""

import numpy as np
import pytest

from widegap import pinn_study as ps

torch = pytest.importorskip("torch")


def test_fixed_configuration_matches_spec():
    c = ps.StudyConfig()
    assert (c.hidden, c.n_colloc, c.resample_every) == ((64,) * 4, 2048, 200)
    assert (c.adam_steps, c.adam_lr, c.lbfgs_max_iter) == (10000, 1e-3, 2000)
    b = ps.parametric_config(3)
    assert (b.adam_steps, b.lbfgs_max_iter, b.n_colloc, b.seed) == (20000, 4000, 2048, 3)


def test_fv_symmetric_junction_has_zero_J_and_one_sided_does_not():
    sol, _ = ps.fv_junction(512, 1.0, 1.0)
    assert abs(ps.trapezoid_J(sol.x, sol.psi)) < 1e-9
    sol2, _ = ps.fv_junction(512, 1.5, 1.0)
    assert ps.trapezoid_J(sol2.x, sol2.psi) > 1.0


def test_adjoint_chain_matches_library_and_richardson_fd():
    n, a = 256, 1.3
    r = ps.fv_J_and_grad(n, a)
    _, g_lib = ps.fv_dJ_dC_library(n, a)
    assert np.max(np.abs(r["dJdC"] - g_lib)) < 1e-12
    fd = ps.richardson_fd_dJda(n, a)
    assert abs(r["dJda"] - fd["R2"]) / abs(fd["R2"]) < 1e-8
    assert r["dJda_contact"] != 0.0  # the Dirichlet datum depends on a


def test_fv_second_order_and_interp_bound():
    xe = ps.eval_grid(4001)
    errs = []
    ref, _ = ps.fv_junction(4096)
    pr = np.interp(xe, ref.x, ref.psi)
    for n in (128, 256):
        s, _ = ps.fv_junction(n)
        errs.append(np.max(np.abs(np.interp(xe, s.x, s.psi) - pr)))
    assert errs[1] < errs[0] / 2.5
    assert ps.interp_error_bound(ref.x, ref.psi) < 1e-3


def test_residual_of_zero_network_is_the_algebraic_charge():
    net = ps.build_net((8,), 1)
    with torch.no_grad():
        for p in net.parameters():
            p.zero_()
    x = torch.tensor([[-30.0], [-1.0], [0.0], [2.0], [30.0]], dtype=torch.float64)
    r = ps.residual(net, x).detach().numpy().ravel()
    t = (x.numpy().ravel() + 40.0) / 80.0
    pl, pr = ps.contact_potentials()
    psi = pl * (1 - t) + pr * t
    c = np.array([-1.0, -1.0, 0.0, 1.0, 1.0])
    expected = ps.DELTA * (np.exp(-psi) - np.exp(psi)) + c  # psi'' = 0 for a zero network
    assert np.allclose(r, expected, rtol=1e-10, atol=1e-12)


def test_tiny_forward_training_has_exact_dirichlet_and_chunked_lbfgs():
    cfg = ps.StudyConfig(hidden=(8, 8), n_colloc=64, adam_steps=20, lbfgs_max_iter=25, seed=1)
    seen = []
    res = ps.train(
        cfg, checkpoint=lambda net: seen.append(1) or {}, adam_ckpt_every=10, lbfgs_chunk=10
    )
    assert res.finite and res.n_loss_evals > 20 and len(seen) >= 2
    p = ps.predict_psi(res.model, np.array([-40.0, 40.0]))
    pl, pr = ps.contact_potentials()
    assert abs(p[0] - pl) < 1e-12 and abs(p[1] - pr) < 1e-12
    assert 0.0 < res.adam_seconds <= res.train_seconds


def test_tiny_parametric_autograd_gradient_matches_finite_difference():
    cfg = ps.StudyConfig(hidden=(8, 8), n_colloc=64, adam_steps=10, lbfgs_max_iter=0, seed=2)
    res = ps.train(cfg, parametric=True)
    a = np.array([0.7, 1.4])
    j, g = ps.surrogate_J_and_grad(res.model, a)
    h = 1e-5
    jp, _ = ps.surrogate_J_and_grad(res.model, a + h)
    jm, _ = ps.surrogate_J_and_grad(res.model, a - h)
    assert np.allclose(g, (jp - jm) / (2 * h), rtol=1e-5, atol=1e-6)
    for av in (0.5, 2.0):
        pl, pr = ps.contact_potentials(av)
        p = ps.predict_psi(res.model, np.array([-40.0, 40.0]), av)
        assert abs(p[0] - pl) < 1e-10 and abs(p[1] - pr) < 1e-10
    # coarse quadrature is consistent with a finer one for a smooth network
    j8, _ = ps.surrogate_J_and_grad(res.model, a, n_quad=4001)
    assert np.allclose(j, j8, atol=1e-3)
