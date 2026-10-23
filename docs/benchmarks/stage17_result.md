# Stage 17 result: gate NOT passed (J4 failed; J1-J3 and J5-J7 pass) (2026-10-23)

Raw data: `results/stage17/` (`summary.json`, `results.json`, `manifest.json`). Figures:
`figures/stage17_theta.png`, `figures/stage17_pareto.png`, `figures/stage17_per_cut_bv.png`,
`figures/stage17_jfet_sensitivity.png`. Pre-registration: `stage17_channel_aware_opt.md`. Follow-up
pre-registration (written after the dev run, before this official result): `stage17b_classification_constraint.md`.
Label: **new** (5-variable re-optimization with the stage-16 channel-mobility R_on and its exact
adjoint). All numbers below are quoted from the JSON files above; anything not recorded there is
marked "not recorded".

## Verdict

The stage pass rule is "J1-J7 hold". **J4 failed, so stage 17 does not pass** (`gate_passed: false`,
`dev: false`). J4 is a smoothness guard, not a constraint; it failed on 113 of 173 evaluated
designs.

**The optimum numbers in this document are not trustworthy as design results.** The binding R_on,ch
constraint was evaluated partly outside the validity of the frozen surface classification (J4, see
below), so "R_on,ch <= R_on,ch(0)" at the stage-17 final design is a statement about the frozen-map
model, not about the model the spec intended. The numbers are reported as recorded, as a statement of
what SLSQP converged to under that model. Stage 17b (pre-registered, not run) enforces the
classification guard as a constraint and is the follow-up.

## Criteria

| ID | Criterion | Measured | Threshold | Result |
|---|---|---|---|---|
| J1 | exact adjoint dR_on,ch/dtheta vs Richardson FD, 5 unit directions, theta = 0, h/2 | rel. error (vector norm) 5.19e-9; FD spread max 2.46e-7; BV (D2) worst rel. 7.83e-9; E_ox (D2) worst rel. 1.43e-12 | R_on <= 1e-5; BV, E_ox <= 1e-4 | pass |
| J2 | gradients at the final design, absolute, in units of the initial values | BV/BV0 3.05e-10; R/R0 3.37e-10; E_ox/E0 7.15e-13 | <= 1e-6 max(1, norm) | pass |
| J3 | feasibility | R/R0 - 1 = 6.53e-13; E_ox/E0 - 1 = 4.55e-12; bound violation -4.66e-3 (negative: inside bounds) | <= 1e-8 | pass |
| J4 | smoothness guards at every evaluated theta | 173 evaluations, **113 violations**; `max_sign_changes` 1 (a surface column changed the sign of N_D - N_A); `min_gated_ns_cm2` 2.836e12; `min_gated_eeff_raw_MV_cm` 0.588; `max_ungated_near_ns_clamp` 0; `max_surface_NA_rel_change` 0.0 | 0 violations | **fail** |
| J5 | dense true-min BV (49 cuts) improves on h and h/2; h/2 gain >= 50 % of h gain | h: 516.87 to 607.61 V (+90.74); h/2: 518.09 to 609.71 V (+91.61); ratio 1.0097 | ratio >= 0.5 | pass |
| J6 | independent h/4 on-state re-check | R_on,ch 2.446694e-3 to 2.446383e-3 ohm cm2, final/initial - 1 = -1.27e-4 | <= +1e-3 | pass |
| J7 | SLSQP exit 0; solves converged; non-finite call accounting | exit 0 ("Optimization terminated successfully"), 15 iterations; 0 failed Picard; 0 failed evaluations; 119 non-finite calls, all inside rejected Newton trial steps; 0 "On entry to" lines; 0 stderr bytes | see section J7 | pass |

Only the clamp margins (N_S, E_eff) and the fixed N_A surface were clean; the single failing part of
J4 is the sign-change guard.

### J4 detail

- **What the guard checks.** The channel model classifies each surface column once, at theta = 0, as
  p-type (inversion, Noguchi mobility) or n-type (bulk mobility). The adjoint is exact for that frozen
  map. J4 asks that no surface column changes the sign of N_D - N_A at any evaluated theta. The
  theta-dependent drift/JFET doping can compensate a p-type surface column at the rounded JFET/body
  corner.
- **Count.** 173 evaluations in the J4 population (every channel on-state evaluation of every job:
  FD points, SLSQP, Pareto, h/4, 13d final); 113 flagged, each with `sign_changes` = 1 (the maximum
  recorded is 1). Per job (`results.json`, `diagnostics`, `guards`): main SLSQP 15 of 17; Pareto
  R/R0 = 1.05, 26 of 26; 1.10, 31 of 33; 1.20, 12 of 43; the stage-17 final design itself is flagged
  (`final.guard.sign_changes` 1, `ok: false`), at h/2 and at h/4 (`ron_pt|final_L1`, `final_L2`). The J1
  FD points at theta = 0 are clean (0 of 4 in each direction); the J2 FD points at the final design
  are all flagged (4 of 4 in each of the 5 jobs).
- **Not recorded in the JSON, quoted from the 17b pre-registration (dev run):** the flagged column is
  the one at the rounded JFET/body corner (x about 0.733 um on h/2); it starts p-type at theta = 0
  with N_D/N_A = 0.93 (so the margin to the flip is 7 % of N_A); it becomes n-type at
  theta_J of about 0.074. The official `results.json` records only the sign-change count per
  evaluation, not the column index, its N_D/N_A or the theta_J crossing. The stage-17 final theta_J
  is 0.1467 (past the crossing), consistent with `sign_changes` = 1.
- **Reclassified R_on at the final design.** Not recorded in `results/stage17/`. The 17b
  pre-registration states, from the dev run, that re-classifying the column at the stage-17 final
  lowers R_on,ch by 4.1 %. This is a dev-run number and has not been re-measured on the official run.
- **Consequence.** The binding constraint R_on,ch <= R_on,ch(0) holds only for the frozen
  classification. Re-classifying the column gives an R_on,ch that differs (4.1 % lower per the
  dev run). The reported BV gain is therefore the gain of a problem whose binding constraint was
  evaluated in a regime the frozen model does not represent. What the BV would be under a consistently
  classified constraint is not quantified here (that is the 17b question).
- **A discrepancy for the Pareto point at R/R0 = 1.2:** its final theta_J is -0.148 (below the
  crossing at about 0.074), yet 12 of its 43 evaluations are flagged. The 17b text says "every Pareto
  design pass that point"; that holds for the 1.05 and 1.10 final designs (theta_J 0.416, 0.400) and
  the stage-17 final, not for the 1.2 final design (the 1.2 run passed through the crossing on the
  way).

## Headline results (quoted, not trustworthy as optima; see Verdict)

Mesh policy as pre-registered: on-state (R_on,ch, adjoint, constraint) on h/2 (4293 nodes);
off-state (soft-min BV, E_ox) on h (1107 nodes). Channel columns 21, gated 12, channel edges 190,
d_ch = 0.02 um, Picard tolerance 1e-10 (maximum 4 Picard iterations).

| Quantity | Initial | Stage-17 final | Stage-13d final (under the channel model) |
|---|---|---|---|
| Soft-min BV (h) | 525.077 V | 617.127 V | 588.286 V |
| Dense true-min BV, h (49 cuts) | 516.871 V | 607.608 V | 577.940 V |
| Dense true-min BV, h/2 | 518.094 V | 609.708 V | 580.413 V |
| True-min BV gain over initial, h / h/2 | | +90.74 V / +91.61 V | +61.07 V / +62.32 V |
| R_on,ch (h/2, ohm cm2) | 2.436882e-3 | 2.436882e-3 (R/R0 - 1 = 6.5e-13, constraint active) | 2.433282e-3 (R/R0 - 1 = -1.48e-3) |
| R_on,ch (h/4, ohm cm2) | 2.446694e-3 | 2.446383e-3 | 2.442934e-3 |
| E_ox at 400 V (h, MV/cm) | 5.2525 | 5.2525 (E/E0 - 1 = 4.5e-12, constraint active) | 5.2525 |

The 13d soft-min value is `bv_obj_13d_final_on_level0`; the 13d R_on,ch is `R_ch_13d_final_h2` and
`R_ch_13d_final_h4` (13d was not run at h/2 in stage 16; here it is evaluated under the channel
model). The 13d h/2 true-min gain is the difference of the recorded values (580.413 - 518.094) and is
a derived number, not a recorded field. The initial and final soft-min BV on h/2 are not recorded.

**theta (5 variables, `R1`):**

| | theta_1 | theta_2 | theta_3 | theta_4 | theta_5 (JFET level) |
|---|---|---|---|---|---|
| stage-17 final | -0.28282 | 0.02068 | 0.81863 | -0.17336 | 0.14666 |
| 13d final | -0.15874 | 0.03536 | 0.85189 | -0.03536 | 0.06798 |
| difference | -0.12408 | -0.01467 | -0.03325 | -0.13800 | +0.07868 |

The JFET level theta_5 is 2.16 times the 13d value and is past the classification crossing (about
0.074, see J4), while the 13d final (0.068) is just below it.

- Stage-17 vs 13d final: +29.67 V on h and +29.30 V on h/2 in dense true-min BV (derived from the
  recorded values), at the same E_ox and with the 13d final 0.148 % below the R_on,ch limit. These
  gains are for the frozen-classification model.
- The binding limiting cut of the final design is x = 1.6 um (`limiting_cut_x_um`); the dense minimum
  is at x = 1.7 um on h and 1.7708 um on h/2.

### Pareto set (R_on,ch levels R/R0 = 1.05, 1.10, 1.20)

| R/R0 (channel, h/2) | Soft-min BV (h) | True-min BV (h) | E_ox / E_ox,0 | SLSQP status | theta_5 |
|---|---|---|---|---|---|
| 1.00 (the final) | 617.13 | 607.61 | 1.0000 | 0 | 0.147 |
| 1.05 | 722.91 | 720.81 | 1.0000 | 0 | 0.416 |
| 1.10 | 731.17 | 732.44 | 0.9591 | 0 | 0.400 |
| 1.20 | 732.99 | 735.42 | 0.6883 | **9 (iteration limit reached)** | -0.148 |

The R/R0 = 1.2 run ended on the iteration limit (40 iterations), not on convergence; the Pareto exit
codes are reported, not gated (as in 13d). 13d Pareto designs (R_on fixed under the bulk model at
1.05, 1.10, 1.20 times) re-evaluated under the channel model have channel R/R0 = 1.0017, 1.0060,
1.0159 and true-min BV 624.4, 651.9, 684.0 V, i.e. they lie at much smaller channel R_on than the
stage-17 Pareto levels and are not a same-constraint comparison. The stage-17 Pareto points are even
less reliable than the final design: 26 of 26 (1.05) and 31 of 33 (1.10) evaluations are flagged by J4.

## R1 (reported) beyond the above

See `figures/stage17_theta.png` (theta) and `figures/stage17_pareto.png` (Pareto vs 13d). Fields at
400 V are in R3. Under the channel model the 13d final has R_on,ch 0.148 % below the initial design
at h/2 and 0.154 % below at h/4 (derived from 2.442934e-3 vs 2.446694e-3), i.e. it is not worse than
the initial design in R_on,ch on either mesh.

## R2 (reported): lagged-gradient error and the JFET-level sensitivity

**Lagged-gradient error** (gradient without F_mu mu_X, relative to the exact adjoint, vector norm):
5.07 % at the initial design and 5.25 % at the final design (`R2_lagged_error`). Per component the
error is 0.26-0.28 % for the four drift coefficients (theta_1-4) and 7.2 % (initial) and 8.1 % (final)
for the JFET level theta_5. So the lagged shortcut would have mis-stated the quantity that matters
most for the optimum by 7-8 %, which is why the exact adjoint was built; the exact gradient matches
its own FD to 5.2e-9 (J1).

**JFET-level sensitivity.** Channel-model dR_on/dtheta_5 relative to the bulk model:

| | theta_1 | theta_2 | theta_3 | theta_4 | theta_5 |
|---|---|---|---|---|---|
| channel/bulk, initial | 0.987 | 0.986 | 0.983 | 0.981 | **2.174** |
| channel/bulk, final | 0.988 | 0.987 | 0.984 | 0.983 | **2.528** |

The four drift coefficients stay within 2 % of the bulk sensitivity (as in the exploratory
diagnostic); theta_5 is 2.17 (initial) and 2.53 (final) times the bulk value. The decomposition of the
lambda-weighted dF/dtheta_5 (`R2`, `theta5`; the two parts sum to `exact_total`, which agrees with
`theta5_check`) at the initial design:

| Part | channel exact | bulk |
|---|---|---|
| edge (mobility) part | 1.392e-5 | 1.523e-5 |
| node (ionization, Poisson) part | -1.082e-4 | -5.86e-5 |

The final design has the same pattern (edge 1.339e-5 vs 1.484e-5; node -1.068e-4 vs -5.18e-5). **The
excess is in the node (ionization/Poisson) part, not in the mobility edges**: the channel model's
mobility-edge contribution is slightly smaller than the bulk one, while the node part is about 1.85
(initial) and 2.06 (final) times larger. By depth band (y below the surface, node part, initial,
channel vs bulk, in 1e-6): 0-0.02 um -7.9 vs -0.6; 0.02-0.05 um -11.4 vs -1.5; 0.05-0.1 um -19.8 vs
-4.2; 0.1-0.3 um -28.0 vs -13.2; 0.3-0.8 um -39.6 vs -37.5; below 0.8 um -1.6 vs -1.6. About 96 % of
the node-part excess (derived from these numbers) lies within 0.3 um of the surface, so the excess sensitivity is the near-surface
JFET region (the Poisson/ionization response of the surface n density that sets N_S, through
mu(N_S)), not the mobility edges. The recorded top edge contributions (initial) are lateral edges at the surface
(y = 0, x = 0.50-0.70 um); the edge part is what the bulk model also has. See
`figures/stage17_jfet_sensitivity.png`. Mechanism interpretation beyond "node part near the JFET surface" was not
tested here.

## R3 (reported): window fields, E_ox, field-line BV

W1-W3 and E_ox at V_blk = 400 V (`R3`, MV/cm; W measures are not mesh-resolvable per 15d, so only
same-mesh comparisons are meaningful):

| Design | Mesh | E_ox | W1 | W2 | W3 |
|---|---|---|---|---|---|
| initial | h | 5.2525 | 2.3877 | 2.6514 | 2.1141 |
| 13d final | h | 5.2525 | 2.2901 | 2.5240 | 2.1125 |
| stage-17 final | h | 5.2525 | 2.1957 | 2.4129 | 2.1106 |
| initial | h/2 | 5.2867 | 2.4466 | 2.7274 | 2.1789 |
| 13d final | h/2 | 5.2819 | 2.3454 | 2.5964 | 2.1754 |
| stage-17 final | h/2 | 5.2894 | 2.2536 | 2.4885 | 2.1764 |
| initial | h/4 | 5.3160 | 2.5342 | 2.7999 | 2.2439 |
| stage-17 final | h/4 | 5.3202 | 2.3405 | 2.5607 | 2.2419 |

Stage-17 vs initial, same mesh (h / h/2): E_ox +4.5e-12 / +5.2e-4 (relative); W1 -8.04 % / -7.89 %; W2
-9.00 % / -8.76 %; W3 -0.17 % / -0.12 %. Stage-17 vs 13d final: W1 -4.12 % / -3.91 %; W2 -4.40 % /
-4.15 %; W3 -0.09 % / +0.04 %. E_ox on h/2 and h/4 rises slightly (+0.05 % h/2 vs initial) although
it is constrained on h only; E_ox is known not to be mesh-converged (15d).

**Field-line BV** (15b method, V) against the dense cut minimum:

| Design | Mesh | Field-line | Dense cut minimum |
|---|---|---|---|
| initial | h | 517.16 | 516.87 |
| initial | h/2 | 516.04 | 518.09 |
| 13d final | h | 577.46 | 577.94 |
| 13d final | h/2 | 577.78 | 580.41 |
| stage-17 final | h | 606.86 | 607.61 |
| stage-17 final | h/2 | 607.13 | 609.71 |

The field-line and cut BV agree to -0.12 % (h) and -0.42 % (h/2) for the stage-17 final; the
stage-17 gain survives (field-line +89.70 V on h, +91.09 V on h/2, derived). All field-line
integrations are `fieldline_converged`; all dense solves are `dense_resolved_all` true.

## J7 detail: non-finite call accounting

From `summary.json` (`values.J7`):

- `slsqp_status` 0, `nit` 15; `failed_picard` 0; `failed_evals_on_state` 0; `failed_evals_off_state` 0;
  `job_errors` empty; `max_picard_iterations` 4.
- `nonfinite_calls_total_over_jobs` **119**: `nonfinite_in_adjoint_solves` 0,
  `nonfinite_in_converged_newton` 0, `nonfinite_in_failed_newton` 119,
  `nonfinite_outside_newton_trial_steps` 0; `failed_newton_calls_rejected_trials` 119.
- Failed Newton *attempts* (rejected cold-march trial steps, all recovered by step halving):
  `on_state_failed_attempts` 111, `off_state_failed_attempts` 8, `off_state_warm_fallbacks` 2.
- `on_entry_to_lines` 0 and `stderr_bytes` 0: no BLAS "On entry to" line, no SuperLU illegal call
  (the stage-16 `safe_spsolve` wrapper is active; non-finite matrices are never passed to SuperLU).
- Pareto exit codes `{"1.05": 0, "1.1": 0, "1.2": 9}` (reported, not gated).

## Caveats

- All stage-16 caveats carry over: Hall factor 1, NO-annealed parameters, no N_it dependence in the
  Coulomb term, 300 K; clamps (N_S 1e11, E_eff 1e-2 MV/cm) are model choices; the stage-16 gate was not
  passed (H5), so the channel-model R_on is not shown mesh-converged. The stage-16 mesh change
  (+3.15 % then +1.27 %) means the constraint is evaluated on h/2 with about 1.3 % mesh uncertainty
  (the h2 to h4 change in R_on,ch here: +0.40 % initial and +0.39 % final).
- **No accumulation-layer mobility.** The n-type surface under the gate (JFET accumulation layer)
  uses the bulk model. The exploratory bracket (fixed 20-50 cm2/Vs) gave +2.0 % to +6.2 % in R_on; the
  sourced bracket is stage 18 (pre-registered, `stage18_accumulation_bracket.md`, Zhao 2007 lower
  end 17.5 cm2/Vs, now in the registry as entry 109). It is pending; until it runs, the R_on,ch
  constraint is optimistic by an unquantified amount.
- **Frozen surface classification (J4).** The result is for the theta = 0 classification; see above.
- h/2 on-state mesh uncertainty of about 1.3 % (stage 16); the on-state is evaluated on h/2 and the
  off-state on h; E_ox and the W measures are not mesh-converged (15d).
- Stage-12c zero-current off-state; vertical-cut and field-line BV, not a full 2D avalanche solution.
- Local optimum of one SLSQP run from one start in a 5-variable family; the 400 V blocking reference
  is a choice. The p-body doping variable (needs a V_T constraint and a punch-through guard) is not
  included.

## Manifest

Code commit `a2425f039250341c5a44518e11d381cb618d50d2`, `git_dirty: false`, 2026-10-23T06:32:08 UTC,
`benchmarks/stage17/run.py --workers 12`, runtime 1211.6 s (summary; `phase2_parallel` 1033.8 s,
`slsqp_main` 146.1 s). Python 3.12.3, NumPy 2.5.4, SciPy 1.18.1, torch 2.14.1+cpu (1 thread), AMD
Ryzen 7 5700U, 16 logical CPUs, Linux 7.0.0-38-generic. `seeds: [0, 1, 2]`; `dev: false`. The manifest
records `widegap 0.2.0` as the package version (the release is v0.4.0; not reconciled here). The results were committed in `b377668`.

## Next

- **Stage 17b** (pre-registered `stage17b_classification_constraint.md`, not run): the same problem
  with SLSQP inequality constraints keeping every gated, p-type-at-theta=0 surface column within a
  2 % margin of its classification on h/2 and h/4 (criteria J4b and J8; R4 reports the BV cost and
  the reclassified R_on of this stage-17 final). The margin m = 0.02 is a model-validity restriction,
  not a physical rule.
- **Stage 18** (pre-registered `stage18_accumulation_bracket.md`, not run): sourced bracket for the
  accumulation-layer mobility.

## Not claimed

- That stage 17 passed: J4 failed.
- That the stage-17 optimum is a valid optimum of the intended model: the binding R_on,ch constraint
  was evaluated partly with a frozen surface classification the final design has crossed.
- The reclassified R_on of the final design (4.1 % is a dev-run figure quoted from the 17b spec).
- Measured-device agreement, global optimality, or mesh-converged channel R_on (stage 16 H5).
- Any accumulation-layer mobility effect.
