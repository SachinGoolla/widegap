# Stage 3b: re-registered convergence tolerance (2026-10-20)

Status: **pre-registered** before the code that implements it.

**Disclosure.** This replaces the stage-3 convergence tolerance *after* run
1 was seen (`docs/benchmarks/stage3_result.md`). Only the convergence test
changes. Discretization, cases and every other criterion (D1–D5, D7, D8) are
unchanged and re-evaluated in a fresh run.

**New tolerance.** For each block $b$,
$\lVert F_b\rVert_\infty\le\tau_b=10^{2}\,\epsilon_{mach}\,S_b$, where

$$
S_b=\max_{i\in b}\Big(\sum_j |J_{ij}\,X_j| + V_i\,|\text{source}_i|\Big),
$$

with $J$ the analytic Jacobian at the current iterate and $X$ the current
iterate. This is the standard first-order round-off bound for evaluating
row $i$, and it includes the $|\psi|/h$ scale that the run-1 definition
omitted. As before, the final update must also satisfy
$\lVert\Delta\psi\rVert_\infty\le10^{-8}$.

Gate unchanged: D1, D2, D3a–D3d, D5 and D6 must hold.
