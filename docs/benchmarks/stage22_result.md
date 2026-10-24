# Stage 22 result: gate PASSED (Q1 to Q3); P22a CONFIRMED (27.7 % R_on gradient error), P22b REFUTED (same optimum, 37 against 10 iterations) (2026-10-23)

Raw data: `results/stage22/` (`summary.json`, `results.json`, `manifest.json`). Figures:
`figures/stage22_gradient_error.png`, `figures/stage22_slsqp_trajectories.png`,
`figures/stage22_final_theta.png`. Pre-registration: `stage22_adjoint_simplifications.md`. Context:
`stage21_result.md` (the E-run and the problem), `stage20_result.md`. All numbers below are quoted from the
JSON files above; anything not recorded there is marked "derived". Device: stage-20 design at
d_body = 1.4 um; figures of merit R_on,ch (h/2), soft-min BV (h), E_ox at 400 V (h). Design points:
theta = 0 and two seeded random points (seeds 0 and 1, `uniform(-0.1, 0.1, 5)` clipped to the 13d box).
Variants: E exact; A drops F_mu mu_X (R_on only); B drops the bulk-mobility doping dependence (R_on only);
C uses complete ionization in dF/dtheta (all three figures of merit); D = A + B + C.

## Verdict

`criteria`: **Q1 true, Q2 true, Q3 true**; `gate_passed: true`, `dev: false`. **Stage 22 passes.**

| ID | Criterion | Measured | Result |
|---|---|---|---|
| Q1 | E matches Richardson FD (eps 1e-3, 5e-4) to <= 1e-5 relative norm, every FoM, every design point | `max_rel_norm_error` 2.54e-7 (tolerance 1e-5); e.g. theta = 0 R_on,ch 8.79e-9, seed 0 BV 3.43e-10 | pass |
| Q2 | variant gradients finite, deterministic, from the same converged state as E | `finite` true, `bitwise_identical_reruns` true, `same_converged_states` true; 2242 bitwise E re-build checks (2197 in the D SLSQP run over 507 assemblies) | pass |
| Q3 | every solve converged; no "On entry to" BLAS lines | `on_entry_to_lines` 0, `job_errors` [], `failed_picard` 0, `failed_evals_on/off` 0, `aborted_slsqp` [], `dense_all_resolved` true | pass |

Non-finite SuperLU calls: `nonfinite_calls_total` 118, all inside rejected (failed) Newton trial steps
(`nonfinite_in_failed_newton` 118; adjoint 0, converged Newton 0).

## R1: gradient error of each variant against E

Relative norm error = ||g_variant - g_E|| / ||g_E|| (`rel_norm_error`), cosine to E (`cosine`). A and B
apply to R_on only, so they have no BV or E_ox entry (B is also exactly zero in the BV and E_ox columns of
`R1`; not shown).

| FoM | Variant | theta = 0 error (cosine) | seed 0 error (cosine) | seed 1 error (cosine) |
|---|---|---|---|---|
| R_on,ch | A | 0.04328 (0.999582) | 0.04657 (0.999527) | 0.04022 (0.999612) |
| R_on,ch | B | 0.18962 (0.999901) | 0.19104 (0.999882) | 0.18553 (0.999900) |
| R_on,ch | C | 0.11840 (0.999971) | 0.12180 (0.999971) | 0.11360 (0.999939) |
| R_on,ch | **D** | **0.27681 (0.999504)** | **0.27920 (0.999411)** | **0.27030 (0.999582)** |
| BV (soft-min) | C | 2.221e-6 (1.0) | 1.873e-6 (1.0) | 2.403e-6 (1.0) |
| BV (soft-min) | D | 2.221e-6 (1.0) | 1.873e-6 (1.0) | 2.403e-6 (1.0) |
| E_ox | C | 3.833e-6 (1.0) | 3.678e-6 (1.0) | 4.101e-6 (1.0) |
| E_ox | D | 3.833e-6 (1.0) | 3.678e-6 (1.0) | 4.101e-6 (1.0) |

(Values rounded to 4 significant digits from `R1`; the A and B rows at the random points are as recorded.)

### Per-component error of the R_on,ch gradient (error / ||g_E||; components 1 to 4 are the drift Legendre coefficients, component 5 the JFET level)

| Variant | theta = 0 | seed 0 | seed 1 |
|---|---|---|---|
| A | -0.0015, 0.0006, -0.0003, 0.0001, **0.0432** | -0.0015, 0.0006, -0.0003, 0.0001, 0.0465 | -0.0015, 0.0007, -0.0004, 0.0001, 0.0402 |
| B | -0.1170, 0.0443, -0.0145, 0.0042, -0.1417 | -0.1181, 0.0443, -0.0145, 0.0043, -0.1427 | -0.1154, 0.0445, -0.0144, 0.0037, -0.1374 |
| C | -0.0646, 0.0234, -0.0086, 0.0018, -0.0960 | -0.0661, 0.0238, -0.0087, 0.0021, -0.0991 | -0.0627, 0.0224, -0.0077, 0.0010, -0.0917 |
| D | -0.1837, 0.0685, -0.0235, 0.0061, **-0.1940** | -0.1861, 0.0690, -0.0237, 0.0065, -0.1947 | -0.1801, 0.0677, -0.0226, 0.0049, -0.1885 |

For D at theta = 0 the largest errors are on the JFET level (-19.4 % of ||g_E||) and the first drift
coefficient (-18.4 %). For BV and E_ox, D at theta = 0 has component errors of order 1e-6 (BV: -1.4e-6,
-1.3e-6, -9.6e-7, -5.7e-7, 0; E_ox: 2.9e-6, 2.2e-6, 1.2e-6, 2.3e-7, 0).

Stage 17 reported about 5 % (8 % on the JFET variable) for dropping F_mu mu_X alone; variant A here gives
4.3 % (theta = 0), with 4.3 % on the JFET component.

## R2: P22a and P22b

### P22a: CONFIRMED

Prediction: for at least one FoM, D differs from E by >= 1 % in relative norm at theta = 0.
`P22a.D_vs_E_rel_norm_theta0`: **R 0.2768 (27.7 %)**, BV 2.22e-6, E_ox 3.83e-6; `threshold` 0.01;
`prediction_holds: true`. The R_on,ch gradient is corrupted by 27.7 % in norm, but its direction is
close (cosine 0.9995);
the BV and E_ox gradients are not materially changed (~1e-6).

### P22b: REFUTED

Prediction: an SLSQP run with D gradients ends at a design differing from the E run in at least one of
(i) dense BV >= 1 V, (ii) a constraint violated by > 1e-6, (iii) KKT ratio (E gradients) >= 10 %,
(iv) exit status != 0. `P22b.conditions`: **all four false**; `prediction_holds: false`.

| Quantity | E run (stage 21 main run, reused) | D run |
|---|---|---|
| exit status | 0 | 0 ("Optimization terminated successfully") |
| iterations (`nit`) | 10 | **37** |
| function evaluations | 10 (derived: length of the recorded evaluation list) | **169** (`n_evals`) |
| wall time | not recorded for the reused run (stage 21 main run 270.4 s) | 1171.1 s |
| dense true-min BV, h (V) | 608.0977 | 608.0550 |
| BV_obj (V) | 608.4223 | 608.3970 |
| R_on,ch (ohm cm2) | 2.4707265053683e-3 | 2.4707265053539e-3 |
| R relative violation | -8.7e-14 | -5.9e-12 |
| E_ox relative violation | -0.010023 | -0.009434 |
| classification: min margin value | 2.6e-16 | 2.8e-17 (violation 0.0) |
| KKT ratio with E gradients | 1.46e-9 (baseline) | **0.008975 (0.90 %)** |
| theta final | (-0.15692, -0.02399, 0.85007, 0.02399, 0.008567) | (-0.16148, -0.03826, 0.85462, 0.03826, 0.008567) |

Differences: dense BV D - E = **-0.0427 V** (limit 1 V); `worst_violation_D` 0.0; `theta_diff_norm`
0.0212; KKT ratio 0.9 % (limit 10 %); both exit status 0. Both runs are feasible and end at the same
active set as far as recorded (R_on,ch active). The dense-BV argmin is at x = 1.8417 um in both. The D run
needs 37 iterations and 169 evaluations against 10 and 10 for the E run (3.7 and 16.9 times, derived).
Trajectories: `figures/stage22_slsqp_trajectories.png`; final theta comparison:
`figures/stage22_final_theta.png`.

The E-run source is stage 21 (`E_run_source` `results/stage21/results.json`, reused, not re-run); the D run
used the same settings (`max_iter` 60, margin m = 0.02, 35 classification constraints).

## Interpretation (not a measured result)

On this problem the simplifications corrupt only the R_on,ch gradient (the BV and E_ox gradients are
affected at the 1e-6 level), and R_on,ch enters SLSQP as a constraint whose value is evaluated exactly at
every iterate. A wrong constraint gradient therefore cannot move the feasible set; it degrades the local
model SLSQP builds of that constraint, which costs convergence speed (37 against 10 iterations, 169
against 10 evaluations) and leaves a small residual stationarity error (KKT ratio 0.9 % with exact
gradients), but the optimum is essentially the same (0.04 V in dense BV). The reduced but nonzero risk that
remains: a problem in which the corrupted gradient belongs to the objective, or to a constraint whose
exactness is not restored by value re-evaluation, could converge to a different or non-stationary design.
Stage 22 did not test such a problem. One device, one start, one coefficient set.

## Caveats

- One device geometry (d_body = 1.4 um), one SLSQP start (theta = 0), two random gradient-comparison points.
- Variants C and D are implemented as in the spec (complete-ionization derivative in dF/dtheta); they are
  not a reproduction of any published code.
- The E-run's wall time and per-evaluation count are not recorded as `n_evals` in the reused stage-21
  record; 10 evaluations is derived from the evaluation list (10 entries, 10 iterations).
- The D run's evaluation count (169) includes SLSQP line-search evaluations.
- Gradient comparison is at three points only; the 27.7 % error is at theta = 0 (27.0 % to 27.9 % across
  the three points).

## Manifest

Code commit `c943a728b2761f92082fb3e64ccbe7d5d12cc285` (`c943a72`), `git_dirty: false`,
2026-10-23T23:58:52 UTC, `benchmarks/stage22/run.py --workers 12`, runtime 1222.60 s. `widegap 0.6.0`
(installed metadata 0.6.0). Python 3.12.3, NumPy 2.5.4, SciPy 1.18.1, torch 2.14.1+cpu (1 thread),
AMD Ryzen 7 5700U, 16 logical CPUs, Linux 7.0.0-38-generic. `seeds` [0, 1]; `dev: false`.

## Not claimed

- That Brungi & Andrei's results are wrong. Their device, model and objective differ, and their
  derivation is to be re-read from the PDF; the motivation here rests on a fetched summary.
- That these simplifications are common practice (the survey is too thin; "at least one recent study").
- That the simplifications are harmless in general: only that they did not change the optimum of this
  problem.
- Anything about problems where the corrupted gradient is in the objective.
