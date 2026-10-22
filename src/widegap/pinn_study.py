"""Stage 11: PINN vs FV/adjoint comparison study (docs/benchmarks/stage11_pinn_question.md).

Q-A: forward PINN for the stage-2 symmetric equilibrium pn junction
    psi'' + delta (e^-psi - e^psi) + C(x) = 0 on [-40, 40], delta = 1e-8,
with a hard Dirichlet ansatz and output scaling for the ~37 V_T potential swing:

    psi(x) = psi_L (1 - t) + psi_R t + t (1 - t) S N(x / 40),   t = (x + 40) / 80,
    S = |psi_R - psi_L|   (so the network output N is O(1)).

Q-B: parametric PINN, same ansatz with the doping magnitude ``a`` as a second
network input; J(a) = int psi dx by a fixed trapezoid rule, dJ/da by autograd.

Doping parametrisation (SPEC ISSUE, see benchmarks/stage11/run.py): for the
symmetric junction C = -a (x<0), +a (x>0) the potential is odd, so J = 0 and
dJ/da = 0 identically. Q-B therefore scales the *n side only*:
C = -1 (x<0), +a (x>0); a = 1 is the stage-2 junction.

The FV reference helpers use ``widegap.fv.equilibrium`` and the stage-5 adjoint
``widegap.adjoint``. ``torch`` is imported lazily so the FV helpers (and their
peak-RSS measurement) do not pay for it.

Chain rule for dJ/da. The adjoint gives dJ/dC_i for interior nodes (C enters
the Poisson rows only). With C_i(a) = cv_average of [-1, a] the doping
derivative is dC_i/da = cv_average of [0, 1] (1 for x>0, 1/2 at the junction
node, 0 for x<0). The contact value psi_R(a) = asinh(a / 2 delta) is a Dirichlet
datum that also depends on a: the Dirichlet row F_N = psi_N - psi_R(a) gives an
extra term + lam_N psi_R'(a), psi_R' = (1 / 2 delta) / sqrt(1 + (a / 2 delta)^2).
Total:  dJ/da = sum_i (dJ/dC_i)(dC_i/da) + lam_N psi_R'(a).
"""

from __future__ import annotations

import time
from collections.abc import Callable
from dataclasses import asdict, dataclass

import numpy as np
from scipy.linalg import solve_banded

from widegap import adjoint
from widegap.fv import equilibrium as eq

DELTA = 1e-8
HALF_WIDTH = 40.0
A_RANGE = (0.5, 2.0)
N_QUAD = 1001  # trapezoid points for J(a) of the surrogate


# ============================================================ FV reference
def junction_doping(x: np.ndarray, a_n: float = 1.0, a_p: float = 1.0) -> np.ndarray:
    """CV-averaged abrupt junction at x = 0: C = -a_p (x<0), +a_n (x>0)."""
    return eq.cv_average_piecewise(x, [0.0], [-a_p, a_n])


def contact_potentials(a_n: float = 1.0, a_p: float = 1.0) -> tuple[float, float]:
    return float(eq.neutral_potential(-a_p, DELTA)), float(eq.neutral_potential(a_n, DELTA))


def fv_junction(n: int, a_n: float = 1.0, a_p: float = 1.0, psi0=None):
    """Stage-2 FV equilibrium solution on a uniform ``n``-interval mesh."""
    x = np.linspace(-HALF_WIDTH, HALF_WIDTH, n + 1)
    c = junction_doping(x, a_n, a_p)
    pl, pr = contact_potentials(a_n, a_p)
    # Newton tolerance on the flux-form residual scales with 1/h (round-off ~ eps*psi/h)
    tol = 1e-12 * max(1.0, n / 1024.0)
    sol = eq.solve(x, c, DELTA, pl, pr, psi0=psi0, tol=tol)
    return sol, c


def trapezoid_J(x: np.ndarray, psi: np.ndarray) -> float:
    return float(np.sum(eq.control_volumes(x) * psi))


def dpsi_right_da(a: float) -> float:
    z = a / (2.0 * DELTA)
    return float((1.0 / (2.0 * DELTA)) / np.sqrt(1.0 + z * z))


def fv_J_and_grad(n: int, a: float, a_p: float = 1.0) -> dict:
    """One FV forward solve + one adjoint solve: J(a) and dJ/da (single transposed solve)."""
    sol, _ = fv_junction(n, a, a_p)
    x, psi = sol.x, sol.psi
    vol = eq.control_volumes(x)
    lam = solve_banded((1, 1), adjoint.banded_transpose(eq.jacobian_banded(psi, x, DELTA)), vol)
    dj_dc = lam * vol
    dj_dc[0] = dj_dc[-1] = 0.0  # contact nodes are not doping design variables
    dc_da = eq.cv_average_piecewise(x, [0.0], [0.0, 1.0])
    d_interior = float(dj_dc @ dc_da)
    d_contact = float(lam[-1] * dpsi_right_da(a))
    return {
        "J": float(np.sum(vol * psi)),
        "dJda": d_interior + d_contact,
        "dJda_interior": d_interior,
        "dJda_contact": d_contact,
        "newton_converged": sol.report.converged,
        "x": x,
        "psi": psi,
        "dJdC": dj_dc,
    }


def fv_dJ_dC_library(n: int, a: float):
    """dJ/dC_i exactly as returned by ``widegap.adjoint.equilibrium_integral_psi``."""
    sol, _ = fv_junction(n, a)
    return adjoint.equilibrium_integral_psi(sol.psi, sol.x, DELTA)


def fv_J_polished(n: int, a: float, psi0=None, steps: int = 3) -> float:
    """Forward-only J(a) with polishing Newton steps (round-off-floor FD evaluations)."""
    sol, c = fv_junction(n, a, psi0=psi0)
    psi = adjoint.polish_equilibrium(sol.psi, sol.x, c, DELTA, steps=steps)
    return trapezoid_J(sol.x, psi)


def richardson_fd_dJda(n: int, a: float, rel_h: float = 0.02, levels: int = 3) -> dict:
    """Richardson-extrapolated central FD of the discrete J(a).

    D(h) central differences at h, h/2, h/4 (h = rel_h * a); R1 eliminates h^2,
    R2 eliminates h^4 as well. Returns D, R1, R2 and the R2-R1 gap as error proxy.
    """
    h = rel_h * a
    d = []
    for k in range(levels):
        hk = h / 2**k
        d.append((fv_J_polished(n, a + hk) - fv_J_polished(n, a - hk)) / (2.0 * hk))
    r1 = [(4.0 * d[k + 1] - d[k]) / 3.0 for k in range(levels - 1)]
    r2 = (16.0 * r1[1] - r1[0]) / 15.0 if levels >= 3 else r1[-1]
    return {"D": d, "R1": r1, "R2": float(r2), "error_proxy": float(abs(r2 - r1[-1]))}


def interp_to_grid(x: np.ndarray, psi: np.ndarray, xe: np.ndarray) -> np.ndarray:
    return np.interp(xe, x, psi)


def interp_error_bound(x: np.ndarray, psi: np.ndarray) -> float:
    """Linear-interpolation error bound h^2/8 max|psi''| from nodal second differences."""
    h = x[1] - x[0]
    d2 = np.abs(psi[2:] - 2.0 * psi[1:-1] + psi[:-2]) / h**2
    return float(h**2 / 8.0 * d2.max())


def eval_grid(n: int = 4001) -> np.ndarray:
    return np.linspace(-HALF_WIDTH, HALF_WIDTH, n)


# ============================================================ PINN
@dataclass(frozen=True)
class StudyConfig:
    """Fixed by the stage-11 spec (Q-A). Q-B uses ``parametric_config`` (budget x 2)."""

    hidden: tuple[int, ...] = (64, 64, 64, 64)
    n_colloc: int = 2048
    resample_every: int = 200
    adam_steps: int = 10000
    adam_lr: float = 1e-3
    lbfgs_max_iter: int = 2000
    seed: int = 0


def parametric_config(seed: int = 0) -> StudyConfig:
    """Q-B: same configuration, optimisation budget x 2 (collocation count unchanged)."""
    return StudyConfig(adam_steps=20000, lbfgs_max_iter=4000, seed=seed)


@dataclass
class TrainResult:
    config: dict
    model: object
    final_loss: float
    train_seconds: float  # excludes checkpoint evaluations
    adam_seconds: float  # the Adam phase alone (for load normalisation)
    n_loss_evals: int
    loss_history: list
    checkpoints: list  # [{"phase","step","train_seconds", ...callback output}]
    finite: bool


def _torch():
    import torch

    return torch


def build_net(hidden: tuple[int, ...], in_dim: int):
    torch = _torch()
    layers: list = []
    w_in = in_dim
    for w in hidden:
        layers += [torch.nn.Linear(w_in, w), torch.nn.Tanh()]
        w_in = w
    layers.append(torch.nn.Linear(w_in, 1))
    return torch.nn.Sequential(*layers).double()


def _bc_torch(a):
    """psi_L, psi_R (torch) for n-side doping ``a`` (tensor) and p-side 1."""
    torch = _torch()
    pl = -torch.asinh(torch.as_tensor(1.0 / (2.0 * DELTA), dtype=torch.float64))
    pr = torch.asinh(a / (2.0 * DELTA))
    return pl, pr


def psi_ansatz(net, x, a=None):
    """psi(x[, a]); x shape (n, 1). ``a=None`` is the fixed a = 1 junction (Q-A)."""
    torch = _torch()
    t = (x + HALF_WIDTH) / (2.0 * HALF_WIDTH)
    xi = x / HALF_WIDTH
    if a is None:
        a_t = torch.ones_like(x)
        inp = xi
    else:
        a_t = a
        inp = torch.cat(
            [xi, (a - 0.5 * (A_RANGE[0] + A_RANGE[1])) / (0.5 * (A_RANGE[1] - A_RANGE[0]))], 1
        )
    pl, pr = _bc_torch(a_t)
    scale = (pr - pl).abs()
    return pl * (1.0 - t) + pr * t + t * (1.0 - t) * scale * net(inp)


def residual(net, x, a=None):
    """psi'' + delta (e^-psi - e^psi) + C at collocation points (the strong-form PDE)."""
    torch = _torch()
    x = x.detach().requires_grad_(True)
    psi = psi_ansatz(net, x, a)
    (d1,) = torch.autograd.grad(psi.sum(), x, create_graph=True)
    (d2,) = torch.autograd.grad(d1.sum(), x, create_graph=True)
    a_t = torch.ones_like(x) if a is None else a
    c = torch.where(x < 0, -torch.ones_like(x), a_t) * (x != 0)
    return d2 + DELTA * (torch.exp(-psi) - torch.exp(psi)) + c


def train(
    cfg: StudyConfig | None = None,
    parametric: bool = False,
    checkpoint: Callable[[object], dict] | None = None,
    adam_ckpt_every: int = 250,
    lbfgs_chunk: int = 100,
) -> TrainResult:
    """Adam then L-BFGS on the PDE residual (the BCs hold exactly through the ansatz).

    ``checkpoint(net)`` may return a dict of diagnostics; its wall time is excluded
    from ``train_seconds``. L-BFGS runs in chunks of ``lbfgs_chunk`` iterations
    (optimizer state persists across ``step`` calls), total <= ``lbfgs_max_iter``.
    """
    torch = _torch()
    cfg = cfg or (parametric_config() if parametric else StudyConfig())
    torch.manual_seed(cfg.seed)
    gen = torch.Generator().manual_seed(cfg.seed)
    net = build_net(cfg.hidden, 2 if parametric else 1)

    def draw():
        x = (torch.rand(cfg.n_colloc, 1, generator=gen, dtype=torch.float64) * 2 - 1) * HALF_WIDTH
        if not parametric:
            return x, None
        a = A_RANGE[0] + (A_RANGE[1] - A_RANGE[0]) * torch.rand(
            cfg.n_colloc, 1, generator=gen, dtype=torch.float64
        )
        return x, a

    history: list[float] = []
    ckpts: list[dict] = []
    clock = {"t": 0.0, "mark": time.perf_counter()}

    def pause():
        clock["t"] += time.perf_counter() - clock["mark"]

    def resume():
        clock["mark"] = time.perf_counter()

    def snap(phase: str, step: int):
        pause()
        rec = {"phase": phase, "step": step, "train_seconds": clock["t"]}
        if checkpoint is not None:
            rec.update(checkpoint(net))
        ckpts.append(rec)
        resume()

    xc, ac = draw()
    opt = torch.optim.Adam(net.parameters(), lr=cfg.adam_lr)
    for step in range(cfg.adam_steps):
        if step and step % cfg.resample_every == 0:
            xc, ac = draw()
        opt.zero_grad()
        loss = residual(net, xc, ac).pow(2).mean()
        loss.backward()
        opt.step()
        history.append(loss.item())
        if checkpoint is not None and (step + 1) % adam_ckpt_every == 0:
            snap("adam", step + 1)

    adam_seconds = clock["t"] + time.perf_counter() - clock["mark"]
    if cfg.lbfgs_max_iter > 0:
        xl, al = draw()  # fixed collocation set for L-BFGS (same count as Adam)
        lbfgs = torch.optim.LBFGS(
            net.parameters(),
            max_iter=lbfgs_chunk,
            line_search_fn="strong_wolfe",
            tolerance_grad=1e-12,
            tolerance_change=1e-15,
            history_size=50,
        )

        def closure():
            lbfgs.zero_grad()
            loss = residual(net, xl, al).pow(2).mean()
            loss.backward()
            history.append(loss.item())
            return loss

        done = 0
        p0 = lbfgs._params[0]
        while done < cfg.lbfgs_max_iter:
            before = lbfgs.state[p0].get("n_iter", 0)
            lbfgs.step(closure)
            used = lbfgs.state[p0]["n_iter"] - before
            done += used
            if checkpoint is not None:
                snap("lbfgs", done)
            if used < lbfgs_chunk:  # converged / line search stalled
                break

    xf, af = draw()
    final = residual(net, xf, af).pow(2).mean().item()
    pause()
    finite = bool(
        np.isfinite(final) and all(np.isfinite(p.detach().numpy()).all() for p in net.parameters())
    )
    return TrainResult(
        config=asdict(cfg),
        model=net,
        final_loss=final,
        train_seconds=clock["t"],
        adam_seconds=adam_seconds,
        n_loss_evals=len(history),
        loss_history=history,
        checkpoints=ckpts,
        finite=finite,
    )


def predict_psi(net, x: np.ndarray, a: float | None = None) -> np.ndarray:
    torch = _torch()
    xt = torch.as_tensor(np.asarray(x, dtype=np.float64)).reshape(-1, 1)
    at = None if a is None else torch.full_like(xt, float(a))
    with torch.no_grad():
        return psi_ansatz(net, xt, at).reshape(-1).numpy()


def surrogate_J_and_grad(net, a_values: np.ndarray, n_quad: int = N_QUAD, chunk: int = 25):
    """J(a) = trapezoid int psi dx and dJ/da by autograd, for a batch of a values."""
    torch = _torch()
    a_values = np.atleast_1d(np.asarray(a_values, dtype=np.float64))
    x = torch.linspace(-HALF_WIDTH, HALF_WIDTH, n_quad, dtype=torch.float64).reshape(1, -1, 1)
    w = torch.full((n_quad,), 2 * HALF_WIDTH / (n_quad - 1), dtype=torch.float64)
    w[0] *= 0.5
    w[-1] *= 0.5
    js, gs = [], []
    for i in range(0, len(a_values), chunk):
        a = torch.as_tensor(a_values[i : i + chunk]).reshape(-1, 1, 1).requires_grad_(True)
        xs = x.expand(a.shape[0], -1, -1)
        psi = psi_ansatz(net, xs.reshape(-1, 1), a.expand(-1, n_quad, -1).reshape(-1, 1))
        j = (psi.reshape(-1, n_quad) * w).sum(1)
        (g,) = torch.autograd.grad(j.sum(), a)
        js.append(j.detach().numpy())
        gs.append(g.reshape(-1).numpy())
    return np.concatenate(js), np.concatenate(gs)
