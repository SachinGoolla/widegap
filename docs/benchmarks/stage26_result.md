# Stage 26 result: gate PASSED (E1 to E3, C1 to C3 hold); Criterion A, P26b and P26c all FAILED (2026-10-24)

Raw data: `results/stage26/` (`summary.json`, `results.json`, `manifest.json`). Figures:
`figures/stage26_a_cross_simulator.png`, `figures/stage26_b_mobility.png`, `figures/stage26_c_bracket.png`,
`figures/stage26_c_oat.png`, `figures/stage26_c_wj_trend.png`. Pre-registration (including its dated Note):
`stage26_saha_comparison.md`. Reference device: Saha and Cooper, IEEE TED 54, 2786 (2007) ("S07").
All numbers are quoted from the JSON files above (`summary.json` unless another file is named); anything not
recorded there is marked "derived". Values in the JSON are in Ohm cm2; they are shown here in mOhm cm2.

## Verdict

`summary.json`: `criteria` C1 true, C2 true, C3 true; `gate_passed: true`, `dev: false`. E1 to E3 are the
26a extension gates; they are **not** recorded in `summary.json` or `results.json`. They are covered by
`tests/test_stage26_extensions.py` (E1 default bit-identity incl. the stage-21 initial on/off state, E2 CSL and
retrograde analytic doping, E3 series-term formulas); the nine tests pass (re-run for this write-up, 41.65 s).
**Stage 26 passes its methodology gate. That is all the gate says.**

**The three scientific outcomes the stage was designed to report all failed:**

- **Criterion A (cross-simulator vs S07 Table I): FAILED** (`R1.a.criterion_A_pass: false`).
- **P26b (mobility model vs S07 Fig. 12): FAILED** (`R1.b.P26b_pass: false`).
- **P26c (measured R_on,sp of 6.95 mOhm cm2): FAILED** (`R1.c.P26c.P26c_pass: false`); bracket labelled WEAK.

**Bottom line: the absolute R_on,sp of this model is NOT validated.** It underpredicts the one fabricated
device checked by 25 % (S07-fit mobility arm) to 46 % (Noguchi arm). The cross-simulator miss (about 26 %)
shows that part of the gap is geometry or contact representation, not physics alone.

| ID | Criterion | Measured | Result |
|---|---|---|---|
| E1 to E3 | 26a extensions (defaults bit-identical; analytic doping; series-term formulas) | `tests/test_stage26_extensions.py`, 9 passed; not in the JSON | pass (test record) |
| C1 | every solve converged; no "On entry to" lines | `C1`: `job_errors` [], `on_entry_to_lines` 0, `off_state_ok` true | pass |
| C2 | central R_on changes <= 3 % from h/2 to h/4 (both arms) | `C2.h4_vs_h2_total_rel`: Noguchi 0.2586 %, S07-fit 0.2353 % (PDE-only 0.563 %, 0.386 %) | pass |
| C3 | Picard convergence | `C3`: `max_iterations` 4, `max_final_rel` 9.81e-11 | pass |

## (a) Cross-simulator comparison and Criterion A

Setup (spec section 2a): S07 simulation conditions, constant mu = 20 cm2/Vs on channel edges, L_ch = 0.5 um,
pitch 8 um, LS = 2 um, W_J = 1 um, uniform body N_A bracketed, no substrate. Source: `R1.a`, `results.json` `a`.

| N_A (cm-3) | R_pde (mOhm cm2) | R_total, rho_c 5e-4 | R_total, rho_c 1e-3 |
|---|---|---|---|
| 2.5e17 | 1.256 | 2.276 | 3.297 |
| 5e17 (central) | 1.869 | 2.889 | 3.910 |
| 1e18 | 2.462 | 3.482 | 4.503 |
| **S07 Table I** | | **3.89** | **5.38** |

Criterion A (Table I within our bracket AND central within +-20 %):

| rho_c | Table I | our bracket | in bracket | central | central_rel_dev | within 20 % |
|---|---|---|---|---|---|---|
| 5e-4 | 3.89 | 2.276 to 3.482 | false | 2.889 | -25.7 % (-0.2573) | false |
| 1e-3 | 5.38 | 3.297 to 4.503 | false | 3.910 | -27.3 % (-0.2733) | false |

`criterion_A_pass: false`. Table I is above the whole N_A bracket for both contact resistivities. h/4 check of
the central PDE part: `a.h4_vs_h2_pde_rel` 2.65 %.

**Contact-geometry diagnostic.** The two Table I rows with the same mobility differ only by the contact term,
so (derived) 5.38 - 3.89 = 1.49 mOhm cm2 for a rho_c change of 5e-4 Ohm cm2, i.e. an implied
`table_I_implied_xmax_over_Lc` = **2.98** (`summary.json`). Our model uses `model_xmax_over_Lc` = **2.041**
(x_max 4 um over L_c 1.96 um). Subtracting the implied contact term (derived: 2.98 x 5e-4 = 1.49) from 3.89
leaves 2.40 mOhm cm2, which lies inside our R_pde bracket (1.256 to 2.462) though above the central 1.869.
This arithmetic assumes S07's total is additive in the same way as ours.

*Interpretation (not tested):* a ratio of 2.98 would arise if the simulated S07 design had a smaller effective
source-contact length for the current path, for instance segmented contacts, or a different current-crowding
geometry than our single contact of length LS = 2 um. The S07 text was not used to confirm this. The pitch
of the non-optimized Table I rows is itself inferred in the spec (section 0).

## (b) Mobility model vs the S07 Fig. 12 fit, and P26b

Gate-field to surface-field conversion (spec section 2b): E_s = (eps_ox / eps_SiC) E_gate with eps_ox = 3.9
and eps_SiC = 10.2 (`results.json` `registry`), so E_s = 0.3824 E_gate; then N_S = eps_SiC E_s / q - N_dep and
E_eff = (q / eps_SiC)(N_dep + N_S / 3). Ratio = Noguchi mu / Fig. 12 fit. Source: `R1.b.ends`.

| N_A (cm-3) | E_gate (MV/cm) | N_S (cm-2) | E_eff (MV/cm) | mu Noguchi | mu fit | ratio | within 1.5 |
|---|---|---|---|---|---|---|---|
| 2.5e17 | 1 | -7.6e11 | excluded (N_S <= 0) | | 13.6 | | |
| 2.5e17 | 2 | 1.39e12 | 0.600 | 15.32 | 11.5 | 1.332 | yes |
| 2.5e17 | 3 | 3.55e12 | 0.728 | 20.21 | 9.0 | 2.246 | no |
| 2.5e17 | 4 | 5.70e12 | 0.855 | 21.54 | 7.0 | 3.077 | no |
| 5e17 | 1 | -2.0e12 | excluded (N_S <= 0) | | 13.6 | | |
| 5e17 | 2 | 1.58e11 | 0.746 | 3.716 | 11.5 | 0.323 | no |
| 5e17 | 3 | 2.31e12 | 0.873 | 13.51 | 9.0 | 1.501 | no (1.5014) |
| 5e17 | 4 | 4.47e12 | 1.001 | 16.21 | 7.0 | 2.316 | no |

**Excluded points.** The 1 MV/cm point has N_S <= 0 for every N_A in the swept range
(`results.json` `b.excluded_points`: 11 entries, N_A 2.5e17 to 5e17, all at E_gate = 1.0). So only three
points (2, 3, 4 MV/cm) are included per N_A (`n_included` 3).

**P26b (within a factor 1.5 at every included point for at least one N_A in the range): FAILED.**
`n_na_passing` = 0 over the 11 swept N_A values (`b.na_plot`). At N_A = 2.5e17 the 3 and 4 MV/cm points are
2.25 and 3.08 times above the fit. At 5e17 the 3 MV/cm ratio is 1.5014, a hair above the threshold, which
does not change the outcome because the 2 and 4 MV/cm points fail by wide margins (0.32 and 2.32).

Reading: the measured constant-current mobility falls from 11.5 to 7 cm2/Vs between 2 and 4 MV/cm gate
field. After the gate-field conversion, E_eff rises only to 0.86 (N_A 2.5e17) or 1.00 MV/cm (5e17) over the
same range, so the Noguchi model, a function of E_eff and N_A, stays high (15 to 22 cm2/Vs) while the
measured mobility falls faster. Caveats stated in the spec still apply: the Fig. 12 values are read from a
figure, Hall (Noguchi) vs constant-current mobility, implanted vs epitaxial p-well, and N_A above the
Noguchi fit range (`counts.na_out_of_range` 3 at the central design).

## (c) Measured R_on,sp (6.95 mOhm cm2) and P26c

Central design: fabricated geometry (pitch 10 um, W_J 1 um, L_ch 0.3 um electrical, gate 3 um, source
contact 3 um), retrograde body 5e17 to 1e18, rho_c 1e-3, R_sub 0.35 mOhm cm2 (range midpoint, a choice made
in the spec Note), V_GS 20 V, V_DS 0.1 V. Source: `R1.c.central`.

| Arm | R_pde | Delta R_c | R_sub | R_total | ratio to measured | channel share | mean mu (gated) |
|---|---|---|---|---|---|---|---|
| Noguchi (main) | 1.735 | 1.689 | 0.35 | **3.774** | **0.543** (-45.7 %) | 0.492 | 20.3 |
| S07 Fig. 12 fit | 3.177 | 1.689 | 0.35 | **5.216** | **0.750** (-25.0 %) | 0.723 | 7.18 |
| measured (S07) | | | | 6.95 | 1 | | |

h/4 totals: 3.784 and 5.228 (C2 above). Picard converged in 4 iterations for both arms; 8 gated columns on
h/2, 15 on h/4 (`n_gated_columns_h2`, `_h4`).

**Bracket (Noguchi arm, 64 corners = 16 PDE x 4 series):** `bracket.min` 2.214, `max` 4.706 mOhm cm2,
`width_ratio` 2.126; `contains_measured: false`. Even the top corner (4.706) is 32 % (derived) below 6.95.
Because the width ratio is >= 2, the comparison is labelled **WEAK** (`weak_label: true`), per the spec rule.

**P26c (6.95 within the Noguchi bracket AND the central within +-20 %): FAILED.** `bracket_contains_measured`
false, `central_rel_dev` -0.457, `central_within_20pct` false. The pre-registered anchor was S07's own 13 %
miss (6.02 vs 6.95); our central miss is 46 % (Noguchi) or 25 % (S07-fit).

**Components.** At the central point the series terms are large: Delta R_c = 1.689 and R_sub = 0.35 mOhm cm2
together make up 54 % (derived, Noguchi arm: 2.039 / 3.774) of the Noguchi total. The `share_channel` values
(0.49 and 0.72) are the PDE channel share as recorded; S07's qualitative "almost equal components" is
consistent with the Noguchi arm only loosely and was not tested as a criterion.

**R_on vs W_J (pitch = 9 + W_J).** Source: `wj_trend`.

| W_J (um) | Noguchi R_total | S07-fit R_total | measured |
|---|---|---|---|
| 1 | 3.774 | 5.216 | 6.95 |
| 1.5 | 3.761 | 5.166 | 7.5 |
| 2 | 3.876 | 5.465 | 7.9 |

Increments from W_J = 1 (`increments_mohm`): Noguchi -0.013 (1 to 1.5) and +0.102 (1 to 2); S07-fit -0.049 and
+0.250; measured +0.55 and +0.95. The model trend is flat (even non-monotonic), while the measurement rises
by 0.95 mOhm cm2. The model does not reproduce the measured JFET-width dependence.

**Sensitivity (one-at-a-time, Noguchi arm; recorded in the official run, `--scan`).** Spans in mOhm cm2
(`oat_ranking`, `oat`): rho_c 0.845 (0.5e-3 to 1e-3), R_sub 0.700 (0 to 0.7), L_ch 0.483 (0.3 to 0.4), N_A,s
0.274, N_A,b 0.101, W_J,eff 0.099. Ranking: rho_c > R_sub > L_ch > N_A,s > N_A,b > W_J,eff. The two largest
are series terms with the least certain inputs (R_sub's upper bound is an unsourced guess). Scan rows (single
runs, `scan`): the literal L_ch mapping (`nominal_l_ch`) gives 3.068 (-0.706 vs central, 3 gated columns);
rounding radius r = 0.05 um gives 3.953 and r = 0.15 um gives 2.476; t_csl = 0.4 um gives 3.688; uniform
N_A = 5e17 gives 3.672. No variant reaches 6.95.

### Off-state of the central design

Source: `c.off_state`, `results.json` `c.off_state`.

- Dense avalanche BV (26 columns, level h, no termination): `bv_min` **1233.0 V** at x = 0.98 um, all columns
  resolved. S06's simulated value is 1230 V (`bv_s06_sim`): +0.24 % (derived). The measured BV of the
  fabricated device is 1050 V (S07), with a trench termination; not gated.
- Electrostatic punch-through (stage-20 method, barrier < 10 V_T): `V_pt_es` = **210.27 V**
  (`V_pt_ceiling: false`, worst column x = 1.365 um).

**Plain statement.** With the assumed body profile (the retrograde 5e17 to 1e18 profile borrowed from M04,
unverified for S07), our model predicts punch-through at about 210 V, whereas the fabricated device blocks
1050 V. Therefore the assumed body profile, the barrier criterion, or both are wrong for this device. This
does not by itself invalidate the stage 19 and 20 results (a different, documented study geometry), but it
means that **the punch-through criterion has not been validated against any measured device**; its first
comparison with one fails by a factor of about 5 (derived: 1050 / 210).

## Caveats

- One fabricated device; several inputs unsourced (S07 bias, temperature, substrate subtraction, body
  implant, gate material); ranges were fixed in the spec and nothing was tuned.
- The channel length is implemented as the electrical channel length (spec Note: chosen after a dev pre-run
  exposed the degeneracy of the literal mapping). Under the literal mapping the central R_on is even lower
  (3.068, scan row `nominal_l_ch`).
- The S07 Fig. 12 curve and Table I pitch are read or inferred, not tabulated.
- The 64-corner bracket reflects only the ranges fixed in the spec; the true uncertainty can be larger
  (the cross-simulator miss shows an unrepresented geometry effect).
- The 26a gate (E1 to E3) is evidenced by the test file, not by the result JSON.

## Manifest

Code commit `6904717b99da511cc6f9730d39a006ed78b89ec9` (`6904717`), `git_dirty: false`,
2026-10-24T05:30:53 UTC, `benchmarks/stage26/run.py --workers 12 --scan`, runtime 478.68 s
(`summary.json` `runtime_seconds`). `widegap 0.6.0`. Python 3.12.3, NumPy 2.5.4, SciPy 1.18.1, torch
2.14.1+cpu (1 thread), AMD Ryzen 7 5700U, 16 logical CPUs, Linux 7.0.0-38-generic. `seeds` [];
`dev: false`.

## Not claimed

- Validation of the model in general, or agreement with any device but this one.
- That the unsourced substrate value, or the assumed M04-borrowed body profile, is correct.
- That the S07 simulated contact geometry is segmented (an untested interpretation).
- A 2D punch-through voltage or a current-carrying breakdown.
