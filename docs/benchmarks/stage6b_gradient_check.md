# Stage 6b: re-registered gradient check at the optimum (2026-10-20)

**Disclosure.** Written after run 1 (`docs/benchmarks/stage6_result.md`).
Only S6a changes. S6b–S6f and the optimizer, problem and parametrization
are unchanged and are re-evaluated in a fresh run.

| ID | Criterion | Gating |
|----|-----------|--------|
| S6a′-init | initial design: relative error ≤ 1e-5 (unchanged) | yes |
| S6a′-final | final design: absolute error $|\text{adjoint}-\text{FD}|\le10^{-6}$ (V_T per unit θ, 100× the budgeted noise floor) along the same 3 directions, **and** FD confirms stationarity, $|\text{FD}|\le10^{-4}$ along each direction | yes |

Stage 6 passes if S6a′-init, S6a′-final and S6b–S6f hold.
