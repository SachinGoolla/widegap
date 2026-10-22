# Stage 11: do PINNs add value next to the verified FV/adjoint stack? (research question)

Status: **pre-registered** (2026-10-22). No stage-11 code exists. This is a
**comparison study**. The gates check the *methodology* (fair, complete and
reproducible comparison), not whether PINNs "win". Either outcome is a
valid result.

## Tasks

- **Q-A (forward, nonlinear).** The stage-2 symmetric pn junction
  ($\delta=10^{-8}$, $[-40,40]$). Reference: the stage-2 FV solution on a
  fine mesh, converged to ≤ 1e-6 relative to the next finer mesh.
  - PINN: the stage-1 architecture family (tanh MLP, float64, hard
    Dirichlet ansatz) with output scaling for the ~37 $V_T$ potential swing.
  - The configuration is **fixed in this document**: 4×64 network, 2048
    collocation points resampled every 200 steps, Adam (1e-3, 10 000
    steps) then L-BFGS (≤ 2000 iterations), seeds 0–4.
- **Q-B (parametric surrogate and sensitivity).** The same junction with
  doping magnitude $a\in[0.5,2]$ as a network input (a parametric PINN,
  same budget × 2). Objective $J(a)=\int\psi\,dx$.
  - Compare $dJ/da$ from automatic differentiation through the trained
    surrogate against the FV adjoint (stage 5) at 9 points in $a$.
  - Compare the **amortized cost**: the time to evaluate $J$ and $dJ/da$
    at $M$ new points, including the PINN's one-off training cost, for
    $M\in\{1,10,10^2,10^3,10^4\}$.

## Acceptance criteria (methodology)

| ID | Criterion | Gating |
|----|-----------|--------|
| M1 | FV references converged: the change between the two finest meshes is ≤ 1e-6 relative for $\psi$ and $J$ | yes |
| M2 | all 5 seeds trained and reported for each task (no selection); failed runs reported as failures | yes |
| M3 | every PINN error is measured on an independent evaluation grid (≥ 4001 points) against the FV reference interpolated with error ≤ 1e-7 | yes |
| M4 | timing: wall time and peak RSS for FV and PINN on the same machine, with thread counts recorded; PINN training time is included in the amortized cost | yes |
| M5 | the FV adjoint $dJ/da$ is verified against Richardson FD to ≤ 1e-8 relative | yes |
| R1 | Q-A: PINN relative L² error and max error (all seeds) vs FV; cost-to-accuracy, i.e. the time each method needs to reach a given error | reported |
| R2 | Q-B: surrogate $J$ and $dJ/da$ errors vs FV (all seeds, all 9 points); the break-even $M$, if any, at which the PINN surrogate becomes cheaper than repeated FV solves at matched accuracy | reported |
| R3 | conclusion stated with evidence, including "no advantage found" if that is the result | reported |

Stage 11 passes if M1–M5 hold.
