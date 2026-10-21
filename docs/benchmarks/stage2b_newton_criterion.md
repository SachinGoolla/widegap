# Stage 2b: re-registered Newton criterion S2′ (2026-10-20)

Status: **pre-registered** before the code that evaluates it.

**Disclosure.** S2′ replaces S2 *after* the stage-2 run was seen
(`docs/benchmarks/stage2_result.md`). A rerun is therefore not a blind test
of the convergence-rate part. To compensate, S2′ adds a check that run 1 did
not contain: the analytic Jacobian is verified directly against finite
differences, and correctness of the Jacobian is what makes Newton quadratic.
S1 and S3–S5 are unchanged and are re-evaluated.

| ID | Criterion | Gating |
|----|-----------|--------|
| S2′a | Analytic Jacobian: for 5 random directions $v$ (seeded) at 3 states per case (the initial guess, a midway iterate, the solution) on E1 ($N=256$) and E2 ($N=1024$), $\lVert J v-[F(\psi+\epsilon v)-F(\psi-\epsilon v)]/(2\epsilon)\rVert_\infty \le 10^{-6}\,\lVert J v\rVert_\infty$ with $\epsilon=10^{-6}$. Dirichlet rows excluded | yes |
| S2′b | Newton converges to $\lVert F\rVert_\infty\le10^{-12}$ within 40 iterations on E1 (all $N$ of S1), E2 and E3 | yes |
| S2′c | E2: every unlimited step whose *output* residual exceeds the round-off floor $\phi=10^3\,\epsilon_{mach}\,s$ has $r_{k+1}/r_k^2\le10^3$, and there are at least 2 such steps. Here $s=\max_i(|F^{flux}_{i+1/2}|+|F^{flux}_{i-1/2}|+|V_i\,\rho_i|)$ at the solution | yes |
| S2′d | Same as S2′c for E1 and E3 | reported |

Stage 2 passes if S1, S2′a, S2′b, S2′c, S3, S4 and S5 hold.
