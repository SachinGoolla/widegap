# Stage 25 result: gate PASSED (H1 to H3 hold); P25a and P25c confirmed, P25b refuted (2026-10-24)

Raw data: `results/stage25/` (`summary.json`, `results.json`, `manifest.json`). Figure:
`figures/stage25_niwa_vs_hatakeyama.png`. Pre-registration: `stage25_generality_hatakeyama.md`.
Context: `stage23_result.md` (Niwa numbers, AB1 and AB3), `stage24_result.md`, `stage21_result.md`.
All numbers below are quoted from the JSON files above (`summary.json` unless another file is named);
anything not recorded there is marked "derived".

## Verdict

`criteria`: **H1 true, H2 true, H3 true**; `gate_passed: true`, `dev: false`.
**Stage 25 passes.** Predictions: **P25a confirmed, P25b refuted, P25c confirmed.**

| ID | Criterion | Measured | Result |
|---|---|---|---|
| H1 | code path reproduces the recorded Niwa values to <= 1e-6 relative | `results.json` `H1`: 13 checks (`n_checks` 13), all `rel_errors` 0.0 (tol 1e-6), including `AB1|initial|dense`, `AB3|13|dense`, `stage21|final|L0` | pass |
| H2 | every dense line resolved or flagged above 2x its first crossing; every solve converged; no "On entry to" lines | `results.json` `H2`: `job_errors` [], `on_entry_to_lines` 0; flagged lines (`n_unbracketed_lines_flagged_above_2x_first_crossing`): `dense|niwa|13|0.8|0` 40, `dense|hatakeyama|13|0.8|0` 40, `ab1|hatakeyama|13` 39, `ab1|hatakeyama|initial` 8, all others 0 | pass (flagged, not resolved; see caveats) |
| H3 | Hatakeyama dense minimum BV of the initial 0.8 um design differs by <= 3 % between h and h/2 | `results.json` `H3`: `bv_h` 536.717 V, `bv_h2` 536.802 V, `rel_change` 1.59e-4 (0.016 %) | pass |

## Tables

### AB1′: stage 13 planar -> dense BV gain, Niwa vs Hatakeyama

Source: `summary.json` `table` rows `AB1[niwa]` and `AB1[hatakeyama]`. Stage 13 initial and final designs,
sharp device, stage-12 off-state, mesh h. Claimed = planar-cut BV gain; real = dense true-minimum BV change;
artifact fraction = (claimed - real) / claimed (as recorded).

| Coefficient set | Claimed planar gain (V) | Real dense gain (V) | Artifact fraction | Planar BV (V) | Dense BV (V) |
|---|---|---|---|---|---|
| Niwa 2015 | 120.91 | 12.41 | 0.897 (89.7 %) | 512.27 -> 633.18 | 406.06 -> 418.47 |
| Hatakeyama 2004 | 126.70 | -98.25 | 1.775 | 521.92 -> 648.62 | 350.41 -> 252.15 |

The Niwa row matches stage 23 `AB1` (120.9 V claimed, 12.4 V real, 89.7 %). Under Hatakeyama the dense
minimum falls: the "optimized" stage-13 final (dense 252.15 V) has a lower dense BV than the initial
(350.41 V).

### AB3′: blocking voltage min(BV, V_pt) vs avalanche BV, 0.8 um body, 12c model, h

Avalanche BV: Niwa from `results/stage23/summary.json` `table` rows `AB3[...]` (`claimed`); Hatakeyama from
`summary.json` `table` rows `AB3[...][hatakeyama]` (`claimed`). V_pt and blocking are the stage-23 DD onset
values (`results/stage23/results.json` `AB3/<design>/dd_onset_V`), unchanged between coefficient sets as the
spec requires. In every design the blocking voltage equals V_pt, so all five are punch-through limited.

| Design | BV Niwa (V) | BV Hatakeyama (V) | V_pt = blocking (V) | Artifact Niwa | Artifact Hatakeyama |
|---|---|---|---|---|---|
| initial | 516.87 | 536.72 | 94.81 | 0.817 | 0.823 |
| 13 | 181.06 | 145.72 | 110.55 | 0.389 | 0.241 |
| 13c | 618.85 | 635.17 | 109.85 | 0.822 | 0.827 |
| 13d | 577.94 | 593.42 | 83.16 | 0.856 | 0.860 |
| 17b | 572.46 | 587.58 | 80.86 | 0.859 | 0.862 |

Artifact fraction = (BV - blocking) / BV (derived from the recorded BV and blocking values; the Hatakeyama
values equal `summary.json` `artifact_fraction`). Niwa values are the stage-23 `AB3` rows; their
`artifact_fraction` as recorded there: 0.817, 0.389, 0.822, 0.856, 0.859.

### AB6: stage 21 optimum on the 1.4 um geometry, dense true-min BV gain, Niwa vs Hatakeyama

Source: `summary.json` `table` rows `AB6[niwa]` and `AB6[hatakeyama]`; `results.json` `predictions` `P25c`.

| Coefficient set | Initial dense BV (V) | Final dense BV (V) | Gain (V) | Difference vs Niwa |
|---|---|---|---|---|
| Niwa 2015 | 544.00 | 608.10 | 64.09 | (reference) |
| Hatakeyama 2004 | 565.43 | 624.49 | 59.06 | -7.85 % (`P25c.rel_difference` -0.0785) |

## Predictions (thresholds unchanged from the spec)

| Prediction | Recorded value | Outcome |
|---|---|---|
| P25a: AB1′ artifact fraction >= 0.5 | `P25a.artifact_fraction` 1.775 (AB1 Hatakeyama) | **CONFIRMED** |
| P25b: every AB3′ design punch-through-limited, artifact >= 0.7 | all five `punch_through_limited` true; artifacts 0.823, 0.241, 0.827, 0.860, 0.862 (`P25b.by_design`) | **REFUTED** by design 13 (0.241 < 0.7) |
| P25c: AB6 Hatakeyama gain positive and within +/-30 % of Niwa | 59.06 V vs 64.09 V, -7.85 % (`P25c.rel_difference`, tolerance 0.3) | **CONFIRMED** |

- **P25a.** Most of the claimed planar-cut gain is absent from the dense minimum under Hatakeyama, and the
  dense BV falls by 98.25 V. The "optimized" stage-13 final is worse than its initial.
- **P25b.** Four of the five designs satisfy the artifact threshold (0.823 to 0.862). Design 13 does not: its
  Hatakeyama avalanche BV is 145.72 V, only about 35 V above its punch-through voltage of 110.55 V
  (derived), so the blocking voltage is set by punch-through and the avalanche-BV artifact is small. The
  prediction assumed the 0.8 um designs would keep a large avalanche-BV artifact; it does not hold for
  design 13 under Hatakeyama. This was pre-registered, so it stands as recorded.
- **P25c.** The stage-21 improvement (dense BV gain of about 60 V) survives the coefficient change, at
  -7.85 % of the Niwa gain.

## Key conclusion

The gain-inflation findings (claim C1a) do not depend on the choice between the Niwa 2015 and Hatakeyama 2004
ionization-coefficient sets, for the quantities tested here. The planar-cut gain of stage 13 is an artifact
under both sets (artifact 0.897 under Niwa, 1.775 under Hatakeyama). The avalanche-BV versus blocking-voltage
gap on the 0.8 um body is large under both sets for four of five designs (0.82 to 0.86). The stage-21 dense
gain on the 1.4 um body keeps its sign and is within 8 % of its Niwa value.

## Interpretation (not a measured result)

Under Hatakeyama the stage-13 planar-cut objective and the dense minimum move in opposite directions, so the
artifact is larger, not smaller. The direction of the effect is therefore the same under both sets, and the
size of the artifact depends on the coefficient set (derived: 0.897 vs 1.775 for AB1; 0.389 vs 0.241 for
design 13 in AB3).

## Caveats

- Two coefficient sets, one device family, one stage-13 pair and five 0.8 um designs (AB3), one 1.4 um
  optimum (AB6). Generality to other coefficient sets or devices is not tested.
- H2 passes because the high-field dense lines are flagged above 2x their first crossing, not resolved:
  40 lines for `dense|niwa|13` and 40 for `dense|hatakeyama|13`, 39 for `ab1|hatakeyama|13`, 8 for
  `ab1|hatakeyama|initial`. The stage-23 convention allows this. The flagged lines are the least resolved
  part of the record; their effect on the recorded BV values was not separately measured.
- V_pt is the stage-23 DD onset (J_n = 1 A/cm2) for both coefficient sets, as the spec pre-registered.
  Punch-through is taken as independent of the ionization coefficients; this was not re-run.
- The blocking voltage is min(dense BV, V_pt) from the electrostatic and DD criteria of stages 19 to 23, not
  a 2D current-carrying breakdown.
- P25b is refuted by design 13 alone; the other four designs are above the threshold.

## Manifest

Code commit `39194667f13bb1f799dae7f05d0aca0433b7cb20` (`3919466`), `git_dirty: false`,
2026-10-24T04:35:19 UTC, `benchmarks/stage25/run.py`, runtime 94.86 s (`summary.json` `runtime_seconds`).
`widegap 0.6.0`. Python 3.12.3, NumPy 2.5.4, SciPy 1.18.1, torch 2.14.1+cpu (1 thread), AMD Ryzen 7 5700U,
16 logical CPUs, Linux 7.0.0-38-generic. `seeds` [] (none recorded); `dev: false`.

## Not claimed

- Generality beyond the Niwa 2015 and Hatakeyama 2004 coefficient sets: other coefficient sets are not tested.
- Generality beyond this device family and geometry set (0.8 um and 1.4 um bodies, stage-12 and stage-13 designs).
- That the Hatakeyama-based avalanche BV or blocking voltage is the correct physical value.
- A 2D punch-through voltage, or agreement with measured devices.
