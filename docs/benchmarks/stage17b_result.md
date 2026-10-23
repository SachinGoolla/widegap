# Stage 17b result: gate PASSED (J1-J3, J4b, J5-J8 all hold) (2026-10-23)

Raw data: `results/stage17b/` (`summary.json`, `results.json`, `manifest.json`). Figures:
`figures/stage17b_theta.png`, `figures/stage17b_pareto.png`, `figures/stage17b_per_cut_bv.png`.
Pre-registration: `stage17b_classification_constraint.md` (written after the stage-17 dev run, before the
official stage-17 result). Comparison data: `results/stage17/` ([stage 17 result](stage17_result.md)) and
`results/stage13d/`. Label: **new** (stage 17 with the surface-classification guard enforced as 35 SLSQP
inequality constraints, margin m = 0.02). All numbers below are quoted from the JSON files above; anything
not recorded there is marked "not recorded" or "derived" (computed from recorded values).

## Verdict

The stage pass rule is "J1, J2, J3, J4b, J5, J6, J7 and J8 hold". All eight hold
(`gate_passed: true`, `dev: false`), so **stage 17b passes**.

The pass is a statement about the constrained problem: the optimum stays inside the region where the
frozen surface classification is valid. **It does not show that the design improves on the 13d optimum**
(see "Comparison to 13d" below: the 17b dense true-min BV is 5.5 V *below* the 13d final on both meshes).

## Criteria

| ID | Criterion | Measured | Threshold | Result |
|---|---|---|---|---|
| J1 | exact adjoint dR_on,ch/dtheta vs Richardson FD, theta = 0, h/2 | rel. error (vector norm) 5.19e-9; BV (D2) worst rel. 7.83e-9; E_ox worst rel. 1.43e-12 | R_on <= 1e-5; BV, E_ox <= 1e-4 | pass |
| J2 | gradients at the final design, absolute, in units of the initial values | BV/BV0 2.66e-13; R/R0 1.76e-10; E_ox/E0 4.05e-9 | <= 1e-6 max(1, norm) | pass |
| J3 | feasibility | R/R0 - 1 = -1.10e-13; E_ox/E0 - 1 = -0.01871 (inactive); bound violation 5.55e-16; min classification constraint value -6.56e-16 | <= 1e-8 | pass |
| J4b | no gated surface column changes sign vs theta = 0 on h/2 or h/4, at every SLSQP iterate, final and Pareto designs; clamp margins as J4 | 93 iterates checked, 0 crossing; 0 final/Pareto designs crossing; 0 clamp-margin violations; 0 line-search trial evaluations crossing (reported); `min_gated_ns_cm2` 2.836e12; `min_gated_eeff_raw_MV_cm` 0.586 | 0 crossings | pass |
| J5 | dense true-min BV (49 cuts) improves on h and h/2; h/2 gain >= 50 % of h gain | h: 516.87 to 572.46 V (+55.59); h/2: 518.09 to 574.91 V (+56.81); ratio 1.0219 (derived) | ratio >= 0.5 | pass |
| J6 | independent h/4 on-state re-check | R_on,ch 2.446694e-3 to 2.446599e-3 ohm cm2, final/initial - 1 = -3.88e-5 | <= +1e-3 | pass |
| J7 | SLSQP exit 0; solves converged; non-finite call accounting | exit 0 ("Optimization terminated successfully"), 11 iterations; Pareto exit codes {1.05: 0, 1.1: 0, 1.2: 0}; 0 failed Picard; 0 failed evaluations (on and off state); 115 non-finite calls, all in rejected Newton trial steps (adjoint 0, converged Newton 0); 0 "On entry to" lines; job_errors empty | see stage 17 J7 | pass |
| J8 | classification-constraint gradients vs FD | rel. error (Frobenius) 3.43e-12 at theta = 0 and 5.07e-12 at the final design; max abs error 2.23e-12 and 2.91e-12; 35 constraints | <= 1e-8 | pass |

The clean J4b is by construction (the constraint enforces it), so it is a consistency check of the
constraint implementation and the SLSQP path, not independent evidence of model validity.

## Headline results

Mesh policy as in stage 17: on-state (R_on,ch, adjoint, constraint) on h/2; off-state (soft-min BV,
E_ox) on h. The classification constraints are placed on both h/2 (level 0) and h/4 (level 1).

| Quantity | Initial | Stage-17b final | Stage-17 final | Stage-13d final (channel model) |
|---|---|---|---|---|
| Dense true-min BV, h (49 cuts) | 516.871 V | **572.464 V** (+55.59) | 607.608 V (+90.74) | 577.940 V (+61.07) |
| Dense true-min BV, h/2 | 518.094 V | **574.906 V** (+56.81) | 609.708 V (+91.61) | 580.413 V (+62.32) |
| Soft-min BV (h) | 525.077 V | 582.891 V | 617.127 V | 588.286 V |
| R_on,ch (h/2, ohm cm2) | 2.436882e-3 | 2.436882e-3 (R/R0 - 1 = -1.1e-13, active) | 2.436882e-3 (6.5e-13) | 2.433282e-3 (-1.48e-3) |
| R_on,ch (h/4, ohm cm2) | 2.446694e-3 | 2.446599e-3 | 2.446383e-3 | 2.442934e-3 |
| E_ox at 400 V (h, MV/cm) | 5.2525 | 5.1542 (E/E0 - 1 = -1.87e-2) | 5.2525 | 5.2525 |

Initial and 17b-final values are from `summary.json` (`bv_obj_*`, `R4`, `J5`, `J6`) and `results.json`
(`final`); stage-17 and 13d columns are from `results/stage17/summary.json` (`R1`). The 17b soft-min BV on
h/2 and the 17b E_ox on h/2 and h/4 are not recorded. Dense argmin cut: x = 1.7 um on h and 1.7708 um on
h/2 (as for the initial design, stage 17 and 13d).

**theta (5 variables):**

| | theta_1 | theta_2 | theta_3 | theta_4 | theta_5 (JFET level) |
|---|---|---|---|---|---|
| stage-17b final | -0.14186 | 0.01346 | 0.83501 | -0.01346 | 0.008567 |
| stage-17 final | -0.28282 | 0.02068 | 0.81863 | -0.17336 | 0.14666 |
| 13d final | -0.15874 | 0.03536 | 0.85189 | -0.03536 | 0.06798 |

Figure: `figures/stage17b_theta.png`. The 17b JFET level theta_5 = 0.008567 sits on the classification
constraint and is 0.058 of the stage-17 value (derived: 0.008567 / 0.14666). The stage-17 final (0.147)
and the 13d final (0.068) are both larger. The 17b theta_4 equals -theta_2 to the printed digits (as in the
13d final); the reason is not analysed here.

### Active constraints

- `R4.active_R`: **true** (R_on,ch = R_on,ch(0), R/R0 - 1 = -1.1e-13).
- `R4.active_E_ox`: **false**. E_ox = 5.1542 MV/cm, E/E0 - 1 = -0.01871; the constraint has slack.
- Classification constraints: exactly one active, `{"level": 1, "x_um": 0.7333333333333335, "value": -6.557254739192331e-16}`, i.e.
  the gated surface column at x = 0.7333 um (the rounded JFET/body corner column of the stage-17 finding)
  **on the h/4 mesh** (level 1 = h/4; level 0 = h/2, per `ClassificationConstraint(evs=[h/2, h/4])` in
  `benchmarks/stage17b/run.py`). The same column on h/2 has value 0.04294 (slack); the next-nearest
  constraints are h/2 x = 1.4022 um (0.2021), h/4 x = 1.4022 um (0.2102), h/4 x = 0.7667 um (0.4104),
  h/4 x = 1.3800 um (0.5010).
- At theta = 0 the smallest constraint value is 0.008359 (`initial.min_constraint_value`), i.e. already close
  to the boundary on h/4. So the constraint is binding because the h/4 corner column is the more
  restrictive one, not because the optimizer was pushed against a loose h/2 limit.

## R4 (reported): cost of the classification constraint, and the stage-17 final reclassified

**BV cost.** Dense true-min BV, 17b vs the stage-17 final: h 572.464 vs 607.608 V, a cost of
`bv_cost_of_constraint_h` = **-35.14 V** (-5.8 %, derived); h/2 574.906 vs 609.708 V, -34.80 V (derived). The
17b gain over the initial design is 55.59 V (h) against 90.74 V (stage 17): about 61 % of the stage-17 gain
survives the constraint (derived). So roughly 39 % of the stage-17 improvement was an artifact of
evaluating the binding R_on constraint outside the validity of the frozen classification.

**Reclassified R_on,ch of the stage-17 final** (channel map rebuilt at the final theta; `R4`):

| Mesh | Frozen map | Reclassified | Change |
|---|---|---|---|
| h/2 | 2.436882e-3 | 2.337085e-3 (20 channel columns, 180 edges) | -4.095 % (`reclass_rel_change_h2` -0.04095) |
| h/4 | 2.446383e-3 | 2.402331e-3 (39 channel columns, 740 edges) | -1.80 % (derived) |

The h/2 value reproduces the 4.1 % figure of the 17b spec (dev run), now on the official run. The
reclassified R_on,ch is *lower*: the column that becomes n-type switches from the inversion mobility to
the bulk mobility, so the frozen-map R_on was pessimistic for the stage-17 final. This means the stage-17
constraint "R_on,ch <= R_on,ch(0)" was conservative at its own final design (the reclassified model would
have allowed more room), yet the model is not self-consistent there; the 17b result is for the
consistent regime. Whether the *reclassified* problem would give a higher BV than stage 17 is not
tested (a dynamic reclassification is not differentiable in the present formulation, spec "Not claimed").
The h/2 and h/4 changes differ (-4.1 % vs -1.8 %), and that difference is itself not explained here.

## Pareto set (R_on,ch levels R/R0 = 1.05, 1.10, 1.20; all with the classification constraints)

| R/R0 (channel, h/2) | Soft-min BV (h) | True-min BV (h) | E_ox / E_ox,0 (derived) | theta_5 | SLSQP status | Stage-17 true-min BV (h) |
|---|---|---|---|---|---|---|
| 1.00 (the final) | 582.89 | 572.46 | 0.9813 | 0.00857 | 0 | 607.61 |
| 1.05 | 717.89 | 715.49 | 0.8029 | 0.00857 | 0 | 720.81 |
| 1.10 | 730.87 | 732.24 | 0.7620 | 0.00857 | 0 | 732.44 |
| 1.20 | 732.99 | 735.42 | 0.6883 | -0.14808 | 0 | 735.42 |

Figure: `figures/stage17b_pareto.png`. All three Pareto runs exit 0, and there were 0 crossings at every
iterate and 0 flagged evaluations (`max_sign_changes` 0). In stage 17 the R/R0 = 1.2 run hit the iteration
limit (exit 9); here it converged (39 iterations). The cost of the constraint shrinks with R_on headroom:
-35.1 V at R/R0 = 1.00 (h), -5.3 V at 1.05, -0.2 V at 1.10, 0 at 1.20 (the 1.2 design has theta_5 = -0.148,
away from the constraint, and the true-min BV agrees with stage 17 to the printed digits). Differences
are derived from recorded values. The stage-17 Pareto points at 1.05 and 1.10 themselves had J4 flags on
most evaluations, so they are not a clean reference.

## Comparison to 13d (bulk-mobility optimum)

The 13d optimum (optimized under the bulk-mobility R_on) re-evaluated under the channel model has dense
true-min BV **577.940 V on h** (580.413 V on h/2) at R/R0 - 1 = -1.48e-3 and E_ox unchanged
(`results/stage17/summary.json`, `R1`).

- **17b is not an improvement over 13d. It is lower**: 572.464 vs 577.940 V on h (-5.48 V, -0.95 %) and
  574.906 vs 580.413 V on h/2 (-5.51 V), derived from recorded values. In R_on,ch the 13d final is also
  lower (2.433282e-3 vs 2.436882e-3 at h/2, 0.148 % lower, derived).
- Classification validity of the 13d final (from `results/stage17/results.json`, the stage-17 re-evaluation jobs): under the frozen stage-16 classification the 13d final flags a surface-sign change in 3 of 4 evaluations on h/2 and 1 of 1 on h/4. Its channel-model R_on is therefore itself evaluated partly outside the classification's validity, so the 13d-vs-17b comparison in the channel model is not like for like. What the record supports: within the region where the model is valid, 17b is the best design found (a local optimum). The 13d design lies outside that region on h/4, and its channel-model figures inherit the same artifact stage 17 had.
- Reading: the 13d final dominates the 17b final on dense true-min BV and on R_on,ch in the channel
  model. That is not a contradiction of optimality: the 13d final was not a feasible point of the 17b
  problem (its theta_5 = 0.068 is far above the 17b limit 0.0086). Whether the 13d final satisfies the 17b
  classification constraints is **not recorded** (the constraint values were not evaluated at the 13d
  final in 17b); its J4 status in stage 17 is in `results/stage17/results.json`. If the 13d final violates
  the m = 0.02 margin, that margin (a model-validity restriction, not physical) is what separates the two
  designs; if it satisfies it, the 17b SLSQP point would be a worse local optimum of a problem that contains
  the 13d final. Which case holds is not determined here.
- Meaningfulness against uncertainty. The 5.5 V difference is a *deficit* of 17b, so no gain over 13d is
  claimed. The R_on,ch difference between the two designs (0.148 %) is far smaller than the h/2 on-state
  mesh uncertainty of about 1.3 % (stage 16: +3.15 % then +1.27 %; here the h/2 to h/4 change is
  +0.39 % at the 17b final and +0.40 % at the initial design), so the R_on,ch ordering of 17b vs 13d is
  not resolved by this model. The open accumulation-layer mobility (below) changes R_on by an unquantified
  amount on all designs; since theta_5 differs between the designs, it would likely affect them differently,
  but how is not computed. So the data support only the statement that, under the channel model with the
  frozen-classification validity constraint, no BV improvement over the bulk-mobility 13d design was found
  at equal R_on,ch.
- The substantive stage-17b conclusion: the large stage-17 improvement over 13d (+29.7 V on h) was
  largely (about 35 V of the 91 V gain over the initial design) a consequence of operating past the
  classification crossing; the consistent-regime optimum gives back about 35 V and ends within 6 V of the
  13d final.

## Caveats

- **The margin m = 0.02 is not physical.** It is a model-validity restriction that keeps the frozen
  classification valid. The 17b optimum depends on the choice of m and on the resolution of the corner
  column (the h/4 column binds, the h/2 one does not).
- **The compensated corner region is unmodelled.** No model of the compensated surface region at the
  rounded JFET/body corner is available (no dynamic reclassification, no smooth inversion/accumulation
  blend, which would need a sourced model). The stage-17 final reclassified R_on,ch (-4.1 % on h/2,
  -1.8 % on h/4) shows the size of the model ambiguity there.
- **All stage-16 caveats carry over**: Hall factor 1, NO-annealed parameters, no N_it dependence in the
  Coulomb term, 300 K; clamps (N_S 1e11, E_eff 1e-2 MV/cm) are model choices; the stage-16 gate was not
  passed (H5), so the channel R_on is not shown mesh-converged (about 1.3 % uncertainty at h/2).
- **No accumulation-layer mobility; stage 18 is pending.** The n-type surface under the gate uses the bulk
  model. The exploratory bracket (20-50 cm2/Vs) gave +2.0 % to +6.2 % in R_on; the sourced bracket
  (stage 18, `stage18_accumulation_bracket.md`, Zhao 2007 lower end 17.5 cm2/Vs) has not run for the 17b
  design. Until it does, the R_on,ch constraint is optimistic by an unquantified amount, and this bears
  directly on the 17b vs 13d comparison.
- Stage-12c zero-current off-state; vertical-cut BV, not a full 2D avalanche solution; E_ox and the W
  measures are not mesh-converged (15d); local optimum of one SLSQP run from one start in a 5-variable
  family; the 400 V blocking reference is a choice; no p-body doping variable.
- The pre-registered J4b counts line-search trial points that cross and reports them; the count is 0.

## Manifest

Code commit `69eb1c367a76dfa891abf56d1fecdd0bbe0e4219`, `git_dirty: false`, 2026-10-23T07:00:34 UTC,
`benchmarks/stage17b/run.py --workers 12`, runtime 1043.8 s (`phase2_parallel` 912.9 s, `slsqp_main`
100.7 s). Python 3.12.3, NumPy 2.5.4, SciPy 1.18.1, torch 2.14.1+cpu (1 thread), AMD Ryzen 7 5700U, 16 logical
CPUs, Linux 7.0.0-38-generic. `seeds: [0, 1, 2]`; `dev: false`. The manifest records `widegap 0.4.0`
(and `widegap_installed_metadata` 0.4.0), i.e. the stage-17 manifest-version defect (0.2.0 recorded) is fixed
here. The stage-17 results used for R4 are read from `results/stage17/results.json`.

## Next

- **Stage 18** (pre-registered, not run): sourced bracket for the accumulation-layer mobility; it should be
  applied to the 13d, stage-17 and 17b designs for a like-for-like comparison.

## Not claimed

- That the 17b design is better than 13d: it is 5.5 V lower in dense true-min BV on both meshes.
- That m = 0.02 is physical, or that the 17b optimum is a global optimum.
- The behavior of the compensated corner region, or the BV of a consistently reclassified optimum.
- Measured-device agreement, mesh-converged channel R_on (stage 16 H5), or any accumulation-layer mobility effect.
