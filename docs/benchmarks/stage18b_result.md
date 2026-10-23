# Stage 18b result: gate NOT passed (B1, B2, B4 hold; B3 mesh fails narrowly, 2.098 % vs 2 %) (2026-10-23)

Raw data: `results/stage18b/` (`summary.json`, `results.json`, `manifest.json`). Figures:
`figures/stage18b_bracket.png`, `figures/stage18b_sheet_charge.png`. Pre-registration:
`stage18b_accumulation_excess.md`. Comparison data: `results/stage18/` ([stage 18 result](stage18_result.md)).
All numbers below are quoted from the JSON files above; anything not recorded there is marked "not recorded"
or "derived" (computed from recorded values). Level 1 = h/2, level 2 = h/4. The L' end applies
mu_acc = 17.5 cm2/Vs (Zhao 2007) only to the excess (gate-induced) electron density on the JFET-surface
edges; the U end is the stage-16 default (bulk Hall mobility). "L" below is the stage-18 end
(mu_acc on all edge electrons in the layer).

## Verdict

The stage pass rule is "B1-B4 hold". B1, B2 and B4 hold; **B3 (mesh) fails narrowly**: R_on(L') of 12b
changes 0.020980596422307762 (2.098 %) from h/2 to h/4 against a 2 % limit (`gate_passed: false`,
`dev: false`). **Stage 18b does not pass.** R1 to R3 are reported as pre-registered.

The bracket numbers below are therefore not mesh-converged for the L' end and must be quoted with a mesh
uncertainty of about 2 % (of R_on(L'), 12b) as stated by B3. The U end reproduces stage 16 exactly (B1).

## Criteria

| ID | Criterion | Measured | Threshold | Result |
|---|---|---|---|---|
| B1 | mu_acc := mu_bulk reproduces stage-16 channel R_on (12b, h/2) | `R_on_B1_this_run` 0.003215047237922234 = `R_on_stage16` 0.003215047237922234; `rel_diff` 0.0 | <= 1e-12 | pass |
| B2 | every Picard and Newton solve converges; no "On entry to" BLAS lines | `failed: []`, `on_entry_to_lines` 0, `stderr_bytes` 0 | none | pass |
| B3 | R_on(L') h/2 to h/4 (12b) | 0.003297876387751674 to 0.00336706780129375 ohm cm2, `rel_diff` 0.020980596422307762 (2.098 %) | <= 2 % | **fail** |
| B4 | R_on(L') d_ch 20 vs 40 nm (12b, h/2) | 0.003297876387751674 to 0.0032994687592602006 ohm cm2, `rel_diff` 0.00048284754226712856 (0.048 %) | <= 1 % | pass |

## Comparison with stage 18

| Criterion | Stage 18 (L) | Stage 18b (L') |
|---|---|---|
| mesh h/2 to h/4 (12b), limit 2 % | A4: 2.71 % (`0.027141360022861605`), fail | B3: 2.10 % (`0.020980596422307762`), fail |
| d_ch 20 vs 40 nm (12b, h/2), limit 1 % | A5: 1.63 % (`0.01626244007477995`), fail | B4: 0.048 % (`0.00048284754226712856`), pass |

Interpretation (not a gate; no run isolates it):

- **d_ch.** Applying mu_acc only to the excess charge removes the d_ch dependence (1.63 % to 0.048 %). This is
  evidence for the stage-18 diagnosis that the d_ch dependence came from the model definition (mu_s applied
  to the equilibrium electrons in the layer). It is consistent with, not a proof of, that diagnosis.
- **Mesh.** The mesh dependence fell (2.71 % to 2.10 %) but did not reach the limit. A part of it is not
  accumulation-specific: R_on(U) of 12b itself moves from 0.003215047237922234 to 0.0032558290366153505
  ohm cm2 from h/2 to h/4, i.e. +1.27 % (derived, 0.012685), which matches the stage-16 record of +1.27 %
  (h/2 to h/4) for the 12b channel-model R_on. The U-model mesh change therefore accounts for roughly 1.3 of
  the 2.10 percentage points of the L' mesh change (derived). The remaining part (about 0.8 percentage
  points, derived by subtraction) is not explained by any run: it could be the number of accumulation edges
  (110 to 460) entering the excess-charge definition, a mesh dependence of the excess charge itself, or
  something else. This is unexplained.
- Because R_on(U) is itself not mesh-converged (stage 16, H5 failed), a mesh criterion on R_on(L') alone
  cannot be expected to pass below the U-model mesh change plus whatever the accumulation part adds.

## R1: bracket [R_on(U), R_on(L')] (ohm cm2)

`complete_results_file` for the 17b design: **true** (`designs.s17b_final`: `present: true`,
`complete_results_file: true`, `dev: false`, source `results/stage17b/results.json`). The stage-17 final is
from `results/stage17/results.json` (`present: true`). "L (st. 18)" is `R_on_L18_ohm_cm2` (stage 18 L, mu_s
on all edge electrons in the layer); "acc. edges" is `n_acc_edges`.

| Design | Level | R_on(U) | R_on(L') | L'/U - 1 | L (st. 18) | L/U - 1 (st. 18) | acc. edges |
|---|---|---|---|---|---|---|---|
| 12b | 1 (h/2) | 3.215047e-3 | 3.297876e-3 | +2.58 % | 3.338163e-3 | +3.83 % | 110 |
| 12b | 2 (h/4) | 3.255829e-3 | 3.367068e-3 | +3.42 % | 3.428765e-3 | +5.31 % | 460 |
| 13d initial | 1 | 2.436882e-3 | 2.515245e-3 | +3.22 % | 2.539322e-3 | +4.20 % | 100 |
| 13d initial | 2 | 2.446694e-3 | 2.558754e-3 | +4.58 % | 2.598883e-3 | +6.22 % | 420 |
| 13d final | 1 | 2.433282e-3 | 2.502895e-3 | +2.86 % | 2.523156e-3 | +3.69 % | 100 |
| 13d final | 2 | 2.398423e-3 | 2.560475e-3 | +6.76 % | 2.608255e-3 | +8.75 % | 440 |
| stage-17 final | 1 | 2.337085e-3 | 2.506143e-3 | +7.23 % | 2.550439e-3 | +9.13 % | 110 |
| stage-17 final | 2 | 2.402331e-3 | 2.550166e-3 | +6.15 % | 2.597068e-3 | +8.11 % | 440 |
| 17b final | 1 | 2.436882e-3 | 2.514058e-3 | +3.17 % | 2.537635e-3 | +4.13 % | 100 |
| 17b final | 2 | 2.446599e-3 | 2.557315e-3 | +4.53 % | 2.597122e-3 | +6.15 % | 420 |

Figure: `figures/stage18b_bracket.png`. The `bulk` column of the stage-18 table is also recorded in
`results/stage18b/summary.json` (`R_on_bulk_ohm_cm2`); it is unchanged and not repeated here. The L/U - 1
values are the recorded `L18_over_U_minus_1`.

Observations (derived from the table):

- **L'/U - 1 ranges from +2.58 % to +7.23 % over all rows** (level 1: +2.58 % to +7.23 %; level 2: +3.42 %
  to +6.76 %). The stage-18 L/U - 1 range was +3.69 % to +9.13 %. L' is below L in every row.
- **The L' bracket is mesh-sensitive and not monotonic across designs.** L'/U - 1 from h/2 to h/4: 12b
  +2.58 to +3.42 %, 13d initial +3.22 to +4.58 %, 13d final +2.86 to +6.76 %, stage-17 final +7.23 to
  +6.15 %, 17b +3.17 to +4.53 %. As in stage 18, the 13d-final level-2 value (+6.76 %) is out of line with
  the 13d-initial (+4.58 %) and 17b (+4.53 %) rows; the cause is not analysed here.
- Ordering of designs under L' is not resolved: the level-2 spread of R_on(L') (2.5502e-3 to 2.5605e-3 across
  the four designs other than 12b) is comparable to the mesh sensitivity above.
- The U values equal those of stage 18 (B1/A2 code path unchanged).

## R2: gate-induced sheet charge vs C_ox(V_GS - V_FB)/q (h/2)

Evaluated for 12b and the 13d initial design at level 1 (h/2) only; other designs and h/4 are not recorded.
Q_ex = integral of (n_bar - n_fb)^+ dy per column. Two flat-band references: V_FB,local (per column, from the
local doping) and V_FB,body (a single value). Both runs converged (`picard_converged: true`, 4 Picard
iterations, `final_rel` 8.96e-12 for 12b and 6.33e-12 for 13d initial). `C_ox` = 6.906266493984e-08 F/cm2
(`t_ox_um` 0.05), V_GS = 15 V, `V_FB_body_V` = -2.611021527117672.

| Design | Columns | Q_ex range (cm-2) | V_FB,local (V) range | ref (local) (cm-2) | ratio Q_ex / ref (local) | ref (body) (cm-2) | ratio Q_ex / ref (body) |
|---|---|---|---|---|---|---|---|
| 12b | 12 | 5.716e12 to 5.853e12 | 0.3124 (constant) | 6.3312e12 | 0.9029 to 0.9245 | 7.5913e12 | 0.7530 to 0.7711 |
| 13d initial | 11 | 5.782e12 to 5.854e12 | 0.2987 to 0.3124 | 6.3312e12 to 6.3371e12 | 0.9125 to 0.9246 | 7.5913e12 | 0.7617 to 0.7711 |

- The excess charge is 0.90 to 0.92 of the local-V_FB reference and 0.75 to 0.77 of the body-V_FB
  reference. Neither ratio is 1; the local definition is closer. The deficit is plausibly the surface
  potential drop / non-ideal charge control (the capacitor formula assumes all of V_GS - V_FB falls across
  the oxide); this is interpretation, not tested.
- The ratio is flat across the interior columns (differences in the 4th digit) and drops only in the last
  columns at the p-body edge (12b: 0.9029 at x = 0.733 um; 13d initial: 0.9125 at x = 0.667 um against 0.9246 at x = 0); this is the
  expected edge effect and is not further analysed.
- Conclusion: the "excess charge" definition behaves as a gate-capacitor charge to within about 8 to 10 %
  (local V_FB) or 23 to 25 % (body V_FB). The consistency check is reported, not gated.

Figure: `figures/stage18b_sheet_charge.png`.

## R3: JFET-level sensitivity dR_on/dtheta_J (FD, 13d initial design, h/2)

FD at eps = 1e-3 and 5e-4 (`spread` 2.3e-11 for U and 1.4e-11 for L'), 5 solves each, all converged.

| Model | R_on at theta = 0 | dR_on/dtheta_J (ohm cm2) | abs(dR/dtheta_J) / R_on (derived) |
|---|---|---|---|
| U (channel model, stage 16) | 2.436882e-3 | -9.426121e-5 | 0.0387 |
| L' (channel model + mu_acc on excess charge) | 2.515245e-3 | -2.250651e-4 | 0.0895 |

`R3.Lp_over_U` = **2.3877** (`2.387674625770438`). Stage 18 reported L/U = 2.9913 for the same quantity
(`R2.ratios` in `results/stage18`); the excess-only definition lowers the factor from about 3.0 to about
2.4. As in stage 18, the stage-16 channel model alone (U) already supplies the stage-17 factor over the bulk
model; L' raises the sensitivity by a further factor of about 2.4. The L' number inherits the B3 mesh
uncertainty; evaluated at one design (13d initial, h/2) and theta = 0 only.

## Bottom line

Within the stage-16 channel model, a JFET-surface accumulation layer at Zhao's un-nitrided mobility
(17.5 cm2/Vs) applied to the gate-induced charge only adds **+2.6 % to +7.2 % to R_on** (L'/U - 1, range over
all designs and both mesh levels; +3.4 % for 12b at h/4). This is a pessimistic, sourced bound (Zhao's oxide
is not nitrided, not a best estimate), not a prediction for an NO-nitrided device. The mesh uncertainty is
about 2 % of R_on(L') (B3, 2.098 % for 12b); part of it (about 1.3 percentage points) comes from the U-model
itself. Stage 18b did not pass B3, so these numbers are not mesh-converged.

## Caveats

- Zhao's oxide is not nitrided and its gate overdrive (about 3-7 V above flat band) is below our 15 V; L' is
  a pessimistic bound.
- No field dependence of mu_acc is modelled; the "excess" definition depends on the flat-band density n_fb
  and agrees with C_ox(V_GS - V_FB)/q only to the ratios in R2.
- L' is mesh-sensitive (B3 failed) and R_on(U) is itself not mesh-converged (stage 16, H5 failed).
- All stage-16 caveats carry over (Hall factor 1, NO-annealed parameters, no N_it dependence, 300 K).
- R2 and R3 are evaluated at h/2 and, for R3, at one design.

## Manifest

Code commit `f09529bf2a61718193e5efb2f047a174cc3e3bea`, `git_dirty: false`, 2026-10-23T18:05:57 UTC,
`benchmarks/stage18b/run.py --workers 12`, runtime 931.5 s (`runtime_seconds` 931.4639711809577). Python
3.12.3, NumPy 2.5.4, SciPy 1.18.1, torch 2.14.1+cpu (1 thread), AMD Ryzen 7 5700U, 16 logical CPUs, Linux
7.0.0-38-generic. `widegap 0.4.0` (installed metadata 0.4.0). `seeds: []`; `dev: false`.

## Not claimed

- Any best-estimate accumulation mobility for NO-nitrided oxide; any field, N_S or temperature dependence.
- That the L' bracket is mesh-converged (B3 failed).
- That the stage-18 diagnosis is proven: the d_ch improvement supports it, the remaining mesh dependence is
  not explained.
- Any ranking of the designs under L' (spread comparable to the mesh sensitivity).
