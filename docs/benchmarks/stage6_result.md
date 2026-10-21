# Stage 6 results

## Run 1, pre-registered criteria (2026-10-20): gate FAILED (S6a only)

Raw data: `results/stage6_run1/`. A first invocation crashed before any
criterion was evaluated. Some designs in the linear-grading grid break down
already at $V_R=0$ ($B_s\ge1$, i.e. $BV_s\le0$), and the runner did not
handle that case. The fix records them explicitly (4 of 41); they cannot be
the grid maximum.

| Criterion | Result | Value |
|---|---|---|
| S6a gradient vs FD, initial design | pass | 1.9e-9 to 5.4e-9 relative |
| S6a gradient vs FD, final design | **fail** | 4.1e-4 to 7.0e-3 relative (absolute 0.8e-8 to 2.5e-8) |
| S6b monotone ascent | pass | 80.626 → 89.755 over 23 accepted iterations |
| S6c constraints at every iterate | pass | $\max|R/R_0-1|=3.3\times10^{-16}$; $b\in[0.62,15.3]$ |
| S6d stopping rule | pass | "gradient" after 23 iterations (24 evaluations) |
| S6e independent exact-BV re-verification | pass | exact BV 80.66 → 89.80 (h=0.025, +9.13); 79.53 → 88.73 (h=0.0125, +9.20) |
| S6f nested linear family | pass | 89.755 ≥ best linear grading 87.951 |

**Diagnosis of S6a.** At the final design the pre-registered stopping rule
*guarantees* a vanishing gradient
($\lVert\nabla_\theta BV_s\rVert=3.1\times10^{-5}$), so the directional
derivatives are only $4\times10^{-6}$ to $2\times10^{-5}$. The absolute
adjoint–FD differences, 0.8e-8 to 2.5e-8, sit exactly at the Brent noise
floor ($10^{-12}/\epsilon=10^{-8}$) budgeted in the pre-registration. That
budget assumed an O(1) signal at both points, which cannot hold at a
stationary point, so the relative measure there is ill-posed. The criterion
is at fault, not the gradient. The gate stays failed. The re-registration is
in `docs/benchmarks/stage6b_gradient_check.md`.

## Run 2, criteria 6b (2026-10-20): gate PASSED

Raw data: `results/stage6b/` (commit `dd3bd40`, clean tree). The run is
bit-identical to run 1: the same $\theta$, the same ascent path, and the same
exact-BV re-verification.

| Criterion | Result | Value |
|---|---|---|
| S6a′-init | pass | ≤ 5.4e-9 relative |
| S6a′-final | pass | \|adjoint − FD\| = 0.8e-8 to 2.5e-8 (≤ 1e-6); \|FD\| = 3.6e-6 to 2.3e-5 (≤ 1e-4), so the design is stationary |
| S6b–S6f | pass | as in run 1 |

**Result** (`figures/stage6_optimized_profile.png`): at a fixed drift
resistance proxy, the optimizer lowers the drift doping near the junction
($b\approx0.62$) and raises it toward the n⁺ layer ($b\approx15$). This
flattens the field, so its integral (the voltage) at the same peak field is
larger. The exact field-based BV rises from 80.66 to 89.80 $V_T$ (+11.3%)
on $h=0.025$, and from 79.53 to 88.73 (+11.6%) on $h=0.0125$. It is the
textbook mechanism behind graded drift layers, found here by verified
adjoint-based optimization.

**Not claimed:**
- a fabrication-ready profile;
- real (avalanche) breakdown;
- that the resistance proxy equals the on-resistance under depletion or
  injection;
- that this is a global optimum (it is a stationary point within a
  5-mode Legendre family that beats the best linear grading);
- anything about 2D or SiC.
