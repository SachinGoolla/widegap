"""Physics-informed neural network for 1D Poisson (docs/benchmarks/poisson_1d.md).

The configuration defaults are the pre-registered benchmark configuration.
Dirichlet data are imposed exactly through the ansatz
u(x) = g0 (1 - x) + g1 x + x (1 - x) N(x), so the loss is the PDE residual
alone. This module deliberately imports nothing from ``wbgtcad.fv``.
"""

from __future__ import annotations

import time
from dataclasses import asdict, dataclass, field

import numpy as np
import torch

from wbgtcad.problems.poisson1d import PoissonProblem1D


@dataclass(frozen=True)
class PINNConfig:
    hidden: tuple[int, ...] = (32, 32, 32)
    n_colloc: int = 256
    resample_every: int = 100
    adam_steps: int = 5000
    adam_lr: float = 1e-3
    lbfgs_max_iter: int = 1000
    lbfgs_colloc: int = 512
    seed: int = 0


@dataclass
class PINNResult:
    config: dict
    model: torch.nn.Module = field(repr=False)
    loss_history: list[float]
    final_loss: float
    train_seconds: float

    def predict(self, problem: PoissonProblem1D, x: np.ndarray) -> np.ndarray:
        xt = torch.as_tensor(np.asarray(x, dtype=np.float64)).reshape(-1, 1)
        with torch.no_grad():
            return _ansatz(self.model, problem, xt).reshape(-1).numpy()


def _mlp(hidden: tuple[int, ...]) -> torch.nn.Sequential:
    layers: list[torch.nn.Module] = []
    width_in = 1
    for w in hidden:
        layers += [torch.nn.Linear(width_in, w), torch.nn.Tanh()]
        width_in = w
    layers.append(torch.nn.Linear(width_in, 1))
    return torch.nn.Sequential(*layers).double()


def _ansatz(net, problem: PoissonProblem1D, x: torch.Tensor) -> torch.Tensor:
    return problem.g0 * (1.0 - x) + problem.g1 * x + x * (1.0 - x) * net(x)


def residual(net, problem: PoissonProblem1D, x: torch.Tensor) -> torch.Tensor:
    """Strong-form residual -(eps u')' - f at collocation points x (shape (n, 1))."""
    x = x.detach().requires_grad_(True)
    u = _ansatz(net, problem, x)
    (du,) = torch.autograd.grad(u.sum(), x, create_graph=True)
    flux = problem.eps(x, torch) * du
    (dflux,) = torch.autograd.grad(flux.sum(), x, create_graph=True)
    return -dflux - problem.f(x, torch)


def train(problem: PoissonProblem1D, cfg: PINNConfig = PINNConfig()) -> PINNResult:
    torch.manual_seed(cfg.seed)
    gen = torch.Generator().manual_seed(cfg.seed)
    net = _mlp(cfg.hidden)

    def draw(n):
        return torch.rand(n, 1, generator=gen, dtype=torch.float64)

    history: list[float] = []
    t0 = time.perf_counter()

    opt = torch.optim.Adam(net.parameters(), lr=cfg.adam_lr)
    xc = draw(cfg.n_colloc)
    for step in range(cfg.adam_steps):
        if step and step % cfg.resample_every == 0:
            xc = draw(cfg.n_colloc)
        opt.zero_grad()
        loss = residual(net, problem, xc).pow(2).mean()
        loss.backward()
        opt.step()
        history.append(loss.item())

    if cfg.lbfgs_max_iter > 0:
        xl = draw(cfg.lbfgs_colloc)
        lbfgs = torch.optim.LBFGS(
            net.parameters(),
            max_iter=cfg.lbfgs_max_iter,
            line_search_fn="strong_wolfe",
            tolerance_grad=1e-12,
            tolerance_change=1e-15,
            history_size=50,
        )

        def closure():
            lbfgs.zero_grad()
            loss = residual(net, problem, xl).pow(2).mean()
            loss.backward()
            history.append(loss.item())
            return loss

        lbfgs.step(closure)

    final = residual(net, problem, draw(cfg.lbfgs_colloc)).pow(2).mean().item()
    return PINNResult(
        config=asdict(cfg),
        model=net,
        loss_history=history,
        final_loss=final,
        train_seconds=time.perf_counter() - t0,
    )
