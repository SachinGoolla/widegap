# Stage 11 result: gate PASSED (M1-M5); no PINN advantage found for this problem (2026-10-22)

Raw data: `results/stage11/` (`summary.json`, `r1.json`, `r2.json`, `fv.json`,
`calibration.json`, `manifest.json`, `pinn_A_seed{0..4}.json`,
`pinn_B_seed{0..4}.json`). Figures: `figures/stage11_qa_solution.png`,
`stage11_qa_cost_accuracy.png`, `stage11_qb_surrogate_errors.png`,
`stage11_qb_amortized_cost.png`. Spec: `stage11_pinn_question.md`
(pre-registered 2026-10-22, with one disclosed amendment, below).
Implemented by a Sonnet subagent from the pre-registered spec, then reviewed.
Label: **new (comparison study)**. The gates check the *methodology*
(fair, complete, reproducible comparison), not whether PINNs "win".

## Verdict

**M1-M5 all pass** (`gate_passed: true`, `failures: []`,
`pre_registered_configuration: true`). The comparison found **no advantage of
the PINN over the FV/adjoint stack on this problem** (R3).

| ID | Criterion | Measured | Threshold | Result |
|---|---|---|---|---|
| M1 | FV references converged, change between the two finest meshes (131072 vs 262144 cells) | Q-A: psi max 3.6e-10, psi L2 9.7e-11 (relative). Q-B, worst of 9 points: psi max 8.7e-10, J 7.8e-10 (rel. to \|J\|; 3.9e-11 rel. to int\|psi\|), dJ/da 1.1e-9 | <= 1e-6 | pass |
| M2 | all 5 seeds trained and reported per task, no selection | seeds 0-4 for Q-A and Q-B, all `finite: true`, `failed: null`, `failures: []` | all seeds, failures reported | pass |
| M3 | independent evaluation grid and reference interpolation error | 4001 evaluation points; interpolation error bound 1.2e-8 (Q-A), at most 2.3e-8 over the nine Q-B points (empirical estimate at most 3.1e-8) | >= 4001 points, <= 1e-7 | pass |
| M4 | wall time and peak RSS for FV and PINN, thread counts recorded | recorded for every run (see Cost and provenance); PINN training time included in the amortized cost | recorded | pass |
| M5 | FV adjoint dJ/da vs Richardson FD (mesh 1024, 3 levels, rel_h 0.02) | worst relative error 3.7e-13 vs the R2 extrapolation (8.5e-10 vs R1) | <= 1e-8 | pass |

## Setup (fixed in the pre-registration, no tuning)

- **PINN:** tanh MLP, 4 x 64, float64, hard Dirichlet ansatz with output
  scaling, 2048 collocation points resampled every 200 steps. Q-A: Adam
  (1e-3, 10 000 steps) then L-BFGS (<= 2000 iterations). Q-B (parametric
  PINN, doping scale $a$ as an input): same budget x 2 (20 000 Adam steps,
  <= 4000 L-BFGS iterations). Seeds 0-4.
- **Reference:** stage-2 FV/Newton solution, 262 144 cells, `ref_seconds_two_meshes`
  0.577 s; Q-B reference $J$ and $dJ/da$ from the FV adjoint (stage 5).

### Disclosed amendment (made before the official run)

**Q-B was degenerate as originally written.** With symmetric doping
$C=\mp a$ the solution is odd, so $J=\int\psi\,dx\equiv0$ and $dJ/da\equiv0$
for all $a$ (confirmed numerically: `symmetric_degeneracy` in `fv.json` has
$J=0.0$ exactly at $a=0.5,\,1.25,\,2$). Q-B was restated: only the n side is
scaled ($C=-1$ for $x<0$, $C=+a$ for $x>0$; $a=1$ is the stage-2 junction),
the nine evaluation points are `linspace(0.5, 2, 9)` (which avoids $a=1$,
where $J=0$), $J$ errors are additionally normalized by $\int|\psi|dx$
because relative $J$ errors near $a\approx1$ are ill-conditioned, and the FV
$dJ/da$ includes the Dirichlet-datum term $\lambda_N\psi_R'(a)$. The PINN
configuration was not changed. Run with `--jobs 1`.

## R1 (reported): Q-A forward problem, stage-2 junction

Per seed, on the independent 4001-point grid (`r1.json`):

| Seed | rel. L2 error | max error | max error / max\|psi\| | train time (s) | final residual MSE | peak RSS (MiB) |
|---|---|---|---|---|---|---|
| 0 | 6.11e-3 | 0.483 | 2.62e-2 | 405.7 | 7.40e-3 | 446.2 |
| 1 | 1.87e-3 | 0.149 | 8.11e-3 | 455.3 | 1.96e-3 | 450.5 |
| 2 | 3.03e-3 | 0.240 | 1.30e-2 | 437.4 | 4.96e-3 | 450.5 |
| 3 | **7.71e-4** (best) | 0.0624 | 3.39e-3 | 447.5 | 1.53e-3 | 452.9 |
| 4 | **1.71e-2** (worst) | 1.344 | 7.30e-2 | 477.3 | 4.29e-2 | 450.1 |

The seed-to-seed spread is more than a factor of 20 in rel. L2 error. The
median training time is 447.5 s (`summary.json` quotes 448 s).

**FV cost to reach the same error** (`fv.json` `ladder_A`, wall time per
solve): N=128 reaches rel. L2 6.23e-4 (max error 0.048) in **1.46e-3 s**,
which is at or below the best PINN error (7.71e-4), against a median PINN
training time of **448 s**. Cost-to-accuracy from `r1.json`:

| Target rel. L2 | FV: N, time | PINN time to first reach it (checkpoints, s), seeds 0-4 |
|---|---|---|
| 1e-1 | 16, 6.0e-4 s | 9.9, 9.6, 9.9, 10.2, 10.2 |
| 1e-2 | 64, 1.13e-3 s | 9.9, 9.6, 9.9, 10.2, 10.2 |
| 1e-3 | 128, 1.46e-3 s | 114.2, 63.9, 215.1, 56.0, 138.0 |
| 1e-4 | 512, 2.92e-3 s | not reached (seed 3 only, 400.3 s) |
| 1e-5 | 1024, 2.67e-3 s | not reached by any seed |

The PINN times are the first logged checkpoint (every 250 Adam steps) at
which the error was at or below the target, so the first two rows are
bounded below by the checkpoint spacing (the first checkpoint, step 250,
is at about 10 s). The FV time at N=1024 is below that at N=512, which is
timing noise at the millisecond level, not a trend. Load factors were 1.0
for all seeds.

## R2 (reported): Q-B parametric surrogate and sensitivity

Nine points $a\in\{0.5, 0.6875, \ldots, 2.0\}$, all 5 seeds
(`r2.json`). Maximum and median over the nine points of the relative error
against the FV reference:

| Seed | J err, max | J err, median | J err / int\|psi\|, max | dJ/da err, max | dJ/da err, median | PINN train (s) | PINN s/query |
|---|---|---|---|---|---|---|---|
| 0 | 5.65e-1 | 6.0e-2 | 6.32e-3 | 1.75e-1 | 9.9e-2 | 995.7 | 2.73e-3 |
| 1 | 4.26e-1 | 8.5e-2 | 6.27e-3 | 5.33e-1 | 1.19e-1 | 938.8 | 2.67e-3 |
| 2 | 4.23e-1 | 2.7e-2 | 4.45e-3 | 2.48e-1 | 7.4e-2 | 943.8 | 3.25e-3 |
| 3 | 9.08e-2 | 2.3e-2 | 4.78e-3 | 2.13e-1 | 5.4e-2 | 938.2 | 2.62e-3 |
| 4 | 6.96e-1 | 1.0e-1 | 4.12e-3 | 2.93e-1 | 1.8e-2 | 865.6 | 2.44e-3 |

Across seeds the worst-point $J$ error is 9.1e-2 to 6.96e-1 and the
worst-point $dJ/da$ error is 1.75e-1 to 5.33e-1 (matching `summary.json`).
Normalized by $\int|\psi|dx$, the $J$ error is at most 6.3e-3, so the large
relative $J$ errors are mostly the ill-conditioning near $J\approx0$
($a\approx1$) that the amendment anticipated. The $dJ/da$ errors, which do
not suffer from this, are large (worst-point errors 1.75e-1 to 5.33e-1). The per-point $\psi$ relative L2 errors are in `r2.json`
(`psi_rel_l2`, range 1.6e-4 to 1.7e-2 over all seeds and points). Peak RSS
for Q-B training was 644 to 672 MiB per seed.

**FV at matched accuracy.** The smallest FV ladder rung whose worst-point
errors do not exceed the PINN's was N=16 for seed 1 and N=32 for seeds 0, 2,
3 and 4 (FV needs N=16 for $J$ and N=32 for $dJ/da$ in seeds 0, 2, 4).
FV cost per evaluation (forward + adjoint) is 1.07e-3 s (N=16) to
1.09e-3 s (N=32). The FV ladder itself (`fv.json` `ladder_B`) reaches
$J$ / $dJ/da$ worst-point errors of 1.8e-1 / 3.8e-1 (N=16, 1.07e-3 s),
3.0e-2 / 6.0e-2 (N=32, 1.09e-3 s), 4.4e-3 / 6.3e-3 (N=64, 1.25e-3 s) and
1.1e-3 / 1.6e-3 (N=128, 1.40e-3 s).

**Amortized cost** (seconds to evaluate $J$ and $dJ/da$ at $M$ new points,
PINN training time included; range over the five seeds, `costs_by_M`;
the "scaled" and "raw" values are identical because all load factors are 1.0):

| M | FV at matched N (s) | PINN incl. training (s) |
|---|---|---|
| 1 | 1.07e-3 to 1.09e-3 | 865.6 to 995.7 |
| 10 | 1.07e-2 to 1.09e-2 | 865.7 to 995.7 |
| 100 | 0.107 to 0.109 | 865.9 to 995.9 |
| 1000 | 1.07 to 1.09 | 868.1 to 998.4 |
| 10 000 | 10.7 to 10.9 | 890.0 to 1023.0 |

### Key finding: no break-even M

`break_even_M` is `null` and `break_even_exists` is `false` for **all five
seeds**. The reason is structural, not a matter of training cost being
large: the PINN's marginal cost per query (2.44e-3 to 3.25e-3 s, median
2.67e-3 s; this includes the autograd derivative) is **not smaller than** the
FV cost per query at matched accuracy (1.07e-3 to 1.09e-3 s at N=16-32). A
surrogate pays off only if its per-query cost is lower than the solver's,
so the one-off training cost (866 to 996 s) is never amortized, at any $M$.
At $M=10^4$ the FV total is about 10.8 s against about 890 to 1023 s for the
PINN. The per-query slope for the 1D FV solve is simply very small, so this
conclusion is specific to a problem where the solver is already
sub-millisecond.

## R3 (conclusion)

**No advantage found for this problem.** For the stage-2 junction (Q-A) the
FV solver reaches the best PINN error in about 1.5 ms where the PINN needs a
median 448 s of training and the error varies by more than 20x across
seeds. For the parametric task (Q-B) the PINN surrogate is both less
accurate (dJ/da errors 1.75e-1 to 5.33e-1 at the worst point) and, per
query, no faster than FV at matched accuracy, so there is no break-even $M$.
The PINN is a legitimate, working solver; it is simply not competitive on
cost or accuracy here.

**Scope limits.** One 1D problem; one fixed, untuned PINN configuration
(4 x 64 tanh, the pre-registered budget); CPU only (AMD Ryzen 7 5700U);
PINN training with 16 torch threads, one job at a time (`--jobs 1`). The
conclusion does not generalize beyond this setting.

## Cost and provenance

- Command `benchmarks/stage11/run.py --jobs 1`, 2026-10-22T19:12:58 UTC.
  Total runtime 6968 s (`total_runtime_seconds` 6968.1). Manifest git commit
  `f19ccce`, `git_dirty: false`.
- Python 3.12.3, widegap 0.2.0, NumPy 2.5.4, SciPy 1.18.1, torch 2.14.1+cpu,
  AMD Ryzen 7 5700U, 16 logical CPUs, Linux 7.0.0-38. Seeds 0-4.
- FV timing at the 262 144-cell reference mesh: forward 0.317 s, forward +
  adjoint 0.392 s. FV peak RSS (isolated subprocess): 76.2 MiB at N=1024,
  108.7 MiB at N=262 144. PINN peak RSS: 446 to 453 MiB (Q-A), 644 to 672 MiB
  (Q-B). Adam step time (calibration): 0.0372 s (Q-A), 0.0361 s (Q-B).
- Thread-count note: the per-run PINN records and `summary.json`
  (`pinn_threads_per_job: 16`) state 16 torch threads, while the top-level
  `manifest.json` and the FV `threads` block record `torch_threads: 8` (the
  parent process). Both are correct: the parent keeps torch's default (8),
  and each PINN subprocess is launched with `--threads 16`, which calls
  `torch.set_num_threads(16)` (`benchmarks/stage11/run.py`). Each PINN JSON
  records 16 from inside its own process. The FV solves use NumPy/SciPy, so
  the torch thread setting does not affect them.

## Not claimed

- That PINNs are useless in general. This is one comparison on one problem.
- Anything about tuned PINNs or other architectures, losses, samplers or
  training schedules (no tuning was done, by design).
- Anything about GPU execution; only CPU was measured.
- Anything about higher-dimensional, parametric-with-many-inputs, or
  expensive-solver problems (for example the 2D device stages), where the
  per-query FV cost is orders of magnitude higher and a surrogate could
  amortize. Not tested.
- Measured-device agreement (not applicable; this is a numerical comparison).
