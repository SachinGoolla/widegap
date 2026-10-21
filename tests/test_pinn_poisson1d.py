"""PINN unit and smoke tests. The gating benchmark (N1) runs via benchmarks/poisson1d/run.py."""

import math

import numpy as np
import pytest

torch = pytest.importorskip("torch")

from widegap.pinn.poisson1d import PINNConfig, residual, train  # noqa: E402
from widegap.problems.poisson1d import evaluate, p1_smooth, p2_variable_eps  # noqa: E402
from widegap.verification import rel_l2_error  # noqa: E402

pytestmark = pytest.mark.pinn


class _ExactBubble(torch.nn.Module):
    """N(x) such that the ansatz reproduces u = sin(pi x) + x exactly."""

    def forward(self, x):
        return torch.sin(math.pi * x) / (x * (1.0 - x))


@pytest.mark.parametrize("make", [p1_smooth, p2_variable_eps])
def test_residual_vanishes_on_exact_solution(make):
    x = torch.linspace(0.05, 0.95, 50, dtype=torch.float64).reshape(-1, 1)
    r = residual(_ExactBubble(), make(), x)
    assert torch.max(torch.abs(r)).item() < 1e-10


def test_ansatz_satisfies_dirichlet_exactly():
    prob = p2_variable_eps()
    res = train(prob, PINNConfig(adam_steps=1, lbfgs_max_iter=0))
    u = res.predict(prob, np.array([0.0, 1.0]))
    assert u[0] == prob.g0 and u[1] == prob.g1


def test_training_is_deterministic_for_fixed_seed():
    prob = p1_smooth()
    cfg = PINNConfig(adam_steps=50, lbfgs_max_iter=0, seed=3)
    a, b = train(prob, cfg), train(prob, cfg)
    assert a.loss_history == b.loss_history


def test_smoke_short_training_p1():
    """CI smoke test: a short run must beat a loose 5e-2 bound.

    Not a benchmark result; the pre-registered N1 criterion uses the full
    configuration in PINNConfig defaults.
    """
    prob = p1_smooth()
    res = train(prob, PINNConfig(adam_steps=300, lbfgs_max_iter=100, seed=0))
    x = np.linspace(0.0, 1.0, 2001)
    assert rel_l2_error(x, res.predict(prob, x), evaluate(prob.u_exact, x)) < 5e-2
