# Stage 3 results

## Run 1 (2026-10-20): gate FAILED (D6 only)

Raw data: `results/stage3_run1/`. Its manifest records `git_dirty: true`
only because a concurrently written, unimported file
(`benchmarks/stage3/devsim_crosscheck.py`) was untracked.

| Criterion | Result | Value |
|---|---|---|
| D1 MMS orders (ψ, n, p) | pass | all 2.00 (last three per field) |
| D2 analytic Jacobian vs central FD | pass | 6.8e-10 to 1.7e-9 (≤ 1e-6) |
| D3a equilibrium vs stage 2 | pass | max\|Δψ\| 2.7e-15; max\|j\| 1.1e-14 |
| D3b current constancy (\|V\| ≥ 3) | pass | ≤ 4.3e-7 (≤ 1e-6) |
| D3c d ln j/dV, V = 4…8 | pass | 1.007, 0.998, 0.994, 0.990 |
| D3d j / j_short, V = 2…8 | pass | 1.000 … 0.994 |
| D5 positivity | pass | all runs |
| D6 Newton convergence | **fail** | D1 at N=512 did not meet the tolerance in 60 iterations. Every other solve converged in ≤ 4 iterations |
| D3e reverse bias (reported) | n/a | j / (−j_short) = 1.00001 at V = −2 and −5 |
| D4 SRH (reported) | n/a | d ln j/dV = 1.07, 0.86, 0.82 for V = 1–2, 2–3, 3–4 |

**Diagnosis of D6.** For D1 at N=512 the ψ-block residual stalls at
2.0e-13 against a tolerance of 1.7e-13. The amended tolerance used
$S_\psi=\max_i\sum|\text{flux terms}|\approx8$, but the round-off of a flux
$-(\psi_{i+1}-\psi_i)/h$ scales with $|\psi|/h\approx10^3$, not with the flux
value. The criterion's scale is wrong. The solver is not: the solution
errors at N=512 continue the exact second-order sequence.

**Note on D4 wording.** The pre-registration calls $d\ln j/dV$ "ideality",
but it is the *inverse* ideality factor $1/n$. D3c's window [0.97, 1.03] is
symmetric about 1, so D3c is unaffected. For D4 the expected trend is toward
0.5, that is $n\to2$. The observed 0.86 and 0.82 ($n\approx1.2$) are
consistent with a growing depletion-region recombination share.
