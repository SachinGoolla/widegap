# Stage 21 result: gate PASSED (S1 to S6); first guarded, punch-through-safe optimized design: dense true-min BV 544.00 to 608.10 V on h (2026-10-23)

Raw data: `results/stage21/` (`summary.json`, `results.json`, `manifest.json`). Figures:
`figures/stage21_pareto.png`, `figures/stage21_per_cut_bv.png`, `figures/stage21_pt_guard.png`.
Pre-registration: `stage21_safe_geometry_opt.md`. Context: `stage20_result.md` (safe body depth),
`stage17b_result.md` (the problem on the 0.8 um geometry). All numbers below are quoted from the JSON files
above; anything not recorded there is marked "derived". Mesh levels: h = L0 (off-state), h/2 = L1,
h/4 = L2. Device: stage-20 design at d_body = 1.4 um (`d_body_um` 1.4, `ny_body` 25), start theta = 0, the
13d/17b variables and bounds, SLSQP (`max_iter` 60), 35 classification constraints (margin m = 0.02).

## Verdict

`criteria`: **S1 to S6 all true**; `gate_passed: true`, `dev: false`. **Stage 21 passes.**

| ID | Criterion | Measured | Result |
|---|---|---|---|
| S1 | gradient checks at theta = 0 and at the final design (17b J1/J2 thresholds) | J1 rel. error R_on,ch 8.79e-9, BV_obj 1.47e-8, E_ox 3.42e-10; J2 abs. (BV/BV0, R/R0, E/E0) 2.43e-11, 1.50e-10, 1.46e-10 against tolerance 1e-6; constraint gradient (J8, 35 constraints) rel. Frobenius error 3.29e-12 (initial), 4.09e-12 (final) | pass |
| S2 | feasibility (R_on,ch, E_ox, classification) | R/R0 - 1 = -8.69e-14; E_ox/E0 - 1 = -0.010023; bound violation 5.55e-16; `n_iterates_crossing` 0, `n_final_designs_crossing` 0, `clamp_violations_at_designs` 0; min constraint value at the final design 2.60e-16 (J3 true, J4b true) | pass |
| S3 | dense true-min BV improves on h and h/2, h/2 gain >= 50 % of h gain | +64.094 V (h), +64.609 V (h/2); ratio 1.0080 | pass |
| S4 | P-guard: V_pt,es >= 1.2 BV at every iterate, final, every Pareto design (h); final also on h/2 | 81 iterates checked, 0 violating; min iterate ratio 1.4419; final 1.4493 (h), 1.4667 (h/2); Pareto 1.4369, 1.4581, 1.4499; `ceiling_hits` 0, `iterate_vpt_errors` 0 | pass |
| S5 | field-line BV within 3 % of dense-cut BV, h and h/2 | rel. diff 0.82 % (h), 1.14 % (h/2) | pass |
| S6 | SLSQP exit 0; all solves converged; no "On entry to" BLAS lines | `slsqp_status` 0 ("Optimization terminated successfully"), `nit` 10; Pareto statuses 0, 0, 0; `on_entry_to_lines` 0; `failed_picard`, `failed_evals_on/off` 0; `job_errors` [] | pass |

Non-finite SuperLU calls: `nonfinite_calls_total` 133, all inside rejected (failed) Newton trial steps
(`nonfinite_in_failed_newton` 133; adjoint 0, converged Newton 0).

### S4: the ambiguity in "BV" and how it was handled

The spec writes the guard as V_pt,es(theta) >= 1.2 BV_obj(theta) in the problem text and as "BV" in
criterion S4, without fixing which BV (soft-min objective or dense true minimum) applies at each design. The
choice is documented in the `benchmarks/stage21/run.py` docstring; `results.json` has no separate
interpretations field for stage 21, so the record is that docstring and the `BV_ref` field of `S4.final_h`.
Handling:

- **SLSQP iterates:** BV_obj (the soft-min objective on h), the only BV available at an iterate (the dense
  BV is not computed per iterate).
- **Final design and Pareto designs:** max(BV_obj, dense true-min BV), the larger (stricter) reference. For
  the final design this is BV_obj = 608.42 V (`BV_ref` 608.4223, dense 608.0977); for the Pareto designs it
  is the dense value (735.07, 744.15, 745.90 V).
- The ratio against the dense BV alone is also recorded for the initial and final designs
  (`V_pt_over_BV_dense`: 1.5724 initial, 1.4501 final). Every ratio is far above 1.2 under either reading,
  so the verdict does not depend on the choice.

## R1: final versus initial

| Quantity | Initial (theta = 0) | Final | JSON field |
|---|---|---|---|
| dense true-min BV, h (V) | 544.0033 | 608.0977 | `bv_dense_h` |
| dense true-min BV, h/2 (V) | 543.1461 | 607.7553 | `bv_dense_h2` |
| soft-min BV_obj, h (V) | 544.0414 | 608.4223 | `bv_obj` |
| field-line BV, h (V) | 541.3099 | 603.0898 | `bv_fieldline_h` |
| field-line BV, h/2 (V) | 538.1088 | 600.8126 | `bv_fieldline_h2` |
| R_on,ch, h/2 (ohm cm2) | 2.470727e-3 | 2.470727e-3 | `R_on_ch_h2` |
| R_on,ch, h/4 (ohm cm2) | 2.480286e-3 | 2.480212e-3 | `R_on_ch_h4` |
| accumulation bracket L' (R_Lp), h/2 (ohm cm2) | 2.549350e-3 | 2.548228e-3 | `accumulation_bracket_Lp_U.*_L1` |
| accumulation bracket L' (R_Lp), h/4 (ohm cm2) | 2.592581e-3 | 2.591223e-3 | `accumulation_bracket_Lp_U.*_L2` |
| V_T, electrostatic (V) | 10.517977 | 10.517970 | `V_T_es_V` |
| E_ox at 400 V, h (MV/cm) | 3.8729 | 3.8341 | `E_ox_400V.*.E_ox_L0` |
| E_ox at 400 V, h/2 (MV/cm) | 3.9080 | 3.8675 | `E_ox_400V.*.E_ox_L1` |
| V_pt,es, h (V) | 855.41 | 881.77 | `V_pt_es.*_h` |
| V_pt,es, h/2 (V) | 866.22 | 892.37 | `V_pt_es.*_h2` |
| V_pt,es / dense BV, h | 1.5724 | 1.4501 | `V_pt_over_BV_dense` |

Dense BV gains: +64.094 V (h; +11.8 %, derived) and +64.609 V (h/2). Field-line gains: +61.78 V (h) and
+62.70 V (h/2) (derived). The R_Lp values are the L' ("excess-charge accumulation mobility") side of the
stage-18b bracket; the U side equals the R_on,ch above (`R_U`). The h/2 on-state mesh uncertainty (about
1.3 %, 16 H5) and the 18b B3 narrow failure apply to the bracket and R_on values.

- **theta final** = (-0.15692, -0.02399, 0.85007, 0.02399, 0.008567) (`theta_final`; initial all 0).
- **Active constraints** (`active_constraints`): one classification constraint, level 1 at x = 0.7333 um
  (value 2.6e-16). The R_on,ch constraint is active (`active_R` true, R/R0 - 1 = -8.7e-14); the E_ox
  constraint is not (`active_E_ox` false, E_ox -1.0 %).
- **Limiting cut of BV_obj** (`cut_limiting_x_um`): x = 2.249 um initially, x = 1.6 um finally. The dense
  minimum of the final design sits at x = 1.8417 um (`dense_x_argmin_um` in the stage-22 record of this run). See
  `figures/stage21_per_cut_bv.png`.
- **Cost of the main run:** 10 iterations, 10 evaluations; `slsqp_main` 270.4 s; total runtime 1538.5 s
  (`runtime_seconds`).

### V_pt,es margin at every iterate (`figures/stage21_pt_guard.png`)

81 iterates were checked (main run 10, Pareto runs 16, 22, 33). Minimum ratio V_pt,es / BV over all
iterates: **1.4419** (Pareto 1.05 run). Main run: 1.4493 to 1.5186 (derived from `iterates`). Pareto runs'
minima: 1.4419 (1.05), 1.4638 (1.1), 1.4469 (1.2). No violation of the 1.2 guard, no ceiling (1200 V) hit,
no V_pt,es evaluation error. The guard was not an SLSQP constraint; it was never active.

### Pareto set (`figures/stage21_pareto.png`)

Designs at fixed R_on,ch / R_0 (SLSQP, same settings, started from theta = 0; all exit status 0).

| R/R_0 | BV_obj (V) | dense BV, h (V) | R_on,ch (ohm cm2) | E_ox 400 V (MV/cm) | V_pt,es (V) | V_pt,es / BV | nit |
|---|---|---|---|---|---|---|---|
| 1.00 (main) | 608.42 | 608.10 | 2.470727e-3 | 3.8341 | 881.77 | 1.4493 | 10 |
| 1.05 | 731.66 | 735.07 | 2.594263e-3 | 3.3318 | 1056.24 | 1.4369 | 16 |
| 1.10 | 741.23 | 744.15 | 2.717799e-3 | 3.2557 | 1085.04 | 1.4581 | 22 |
| 1.20 | 743.22 | 745.90 | 2.965e-3 (derived: R_0 x 1.2) | 3.2649 | 1081.50 | 1.4499 | 33 |

The R/R_0 = 1.0 row is the main run; the 1.2 row's R_on,ch is quoted from `R_over_R0` 1.2000 and R_0.
The step from the main run to R/R_0 = 1.05 gains 127.0 V dense BV (derived); beyond 1.05 the gain flattens
(+9.1 V to 1.10, +1.7 V to 1.20; derived).

## R2: comparison with stage 17b (different geometry)

**Different geometry:** stage 17b used the 0.8 um body, which punches through at about 90 to 130 V
(stage 19); stage 21 uses d_body = 1.4 um. The designs are not interchangeable and the BV values are not
a like-for-like comparison of optimizers. Values from `R2` in `results.json`:

| | 17b, 0.8 um | 21, 1.4 um |
|---|---|---|
| initial dense BV, h (V) | 516.87 | 544.00 |
| final dense BV, h (V) | 572.46 | 608.10 |
| gain (V, derived) | +55.59 | +64.09 |
| BV_obj initial to final (V) | 525.08 to 582.89 | 544.04 to 608.42 |
| theta final | (-0.14186, 0.01346, 0.83501, -0.01346, 0.008567) | (-0.15692, -0.02399, 0.85007, 0.02399, 0.008567) |

The stage-17 final (without the classification constraint) had dense BV 607.61 V on h (`R2`, key `17`).
The two theta vectors are similar in the first, third and fifth components; the second and fourth have
opposite sign (small magnitudes, 0.013 and 0.024). Whether the similarity says anything about the design
problem was not tested.

## Caveats

- Carried over from the spec: stage-16 channel-model caveats; h/2 on-state mesh uncertainty of about 1.3 %
  (16 H5); the 18b B3 narrow failure (the bracket values are reported, not gated).
- Punch-through safety above 400 V rests on the electrostatic 10 V_T barrier criterion (stage 20); the DD
  runs there reached only 400 V for d_body >= 1.2 um.
- The soft-min BV_obj and the dense BV differ by 0.3 V at the final design (608.42 vs 608.10); the guard
  at iterates uses BV_obj.
- BV (zero-current avalanche) and V_pt,es (electrostatic barrier) are evaluated separately; there is no
  coupled current-carrying breakdown simulation.
- Local optimum from one start (theta = 0) in a 5-parameter family; the 400 V E_ox reference is a choice.
- The R_on,ch constraint is active at the optimum; its h/2 value is not the converged h/4 value
  (2.4707e-3 vs 2.4802e-3).

## Manifest

Code commit `63b92b137bb4bfb28d0e1dd6713a47839163123e` (`63b92b1`), `git_dirty: false`,
2026-10-23T23:33:09 UTC, `benchmarks/stage21/run.py --workers 12`, runtime 1538.47 s. `widegap 0.6.0`
(installed metadata 0.6.0). Python 3.12.3, NumPy 2.5.4, SciPy 1.18.1, torch 2.14.1+cpu (1 thread),
AMD Ryzen 7 5700U, 16 logical CPUs, Linux 7.0.0-38-generic. `seeds` [0, 1, 2]; `dev: false`.

## Not claimed

- That the design is punch-through-free above 400 V in drift-diffusion (the claim rests on the
  electrostatic barrier criterion).
- A 2D current-carrying breakdown simulation, fabrication feasibility, or agreement with measured devices.
- That the BV gain is a property of the optimizer rather than of the geometry: the 17b and 21 numbers are
  on different geometries (R2).
- Global optimality; this is one SLSQP start in a 5-parameter family.
