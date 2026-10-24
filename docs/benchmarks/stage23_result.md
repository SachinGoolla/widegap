# Stage 23 result: gate PASSED (T1 to T3); artifact-budget table AB1 to AB5 (2026-10-24)

Raw data: `results/stage23/` (`summary.json`, `results.json`, `manifest.json`). Figure:
`figures/stage23_artifact_budget.png`. Pre-registration: `stage23_artifact_budget.md`. Context:
`stage22_result.md`, `stage19_result.md`, `stage20_result.md`, `stage17b_result.md`, `stage16_result.md`,
`stage13_result.md`. All numbers below are quoted from the JSON files above (`summary.json` `table[...]`
unless another file is named); anything not recorded there is marked "derived". For each entry the artifact
fraction is `1 - real/claimed`, applied to the claimed quantity (a gain for AB1, AB2, AB5; an absolute value
for AB3, AB4), as pre-registered. "Claimed" and "real" are compared on the same device, mesh and design.

## Verdict

`criteria`: **T1 true, T2 true, T3 true**; `gate_passed: true`, `dev: false`. **Stage 23 passes.**

| ID | Criterion | Measured | Result |
|---|---|---|---|
| T1 | AB1 recomputed planar-cut gain reproduces the stage-13 record to <= 1e-6 relative | `T1.rel_errors` bv_initial 0.0, bv_final 0.0, gain 0.0 (512.2726 and 633.1832 V; recorded `results/stage13/summary.json` `values/bv_initial_V`, `values/bv_final_V`) | pass |
| T2 | 1D DD reaches the onset for every design; V_pt,es changes <= 5 % from h to h/2 | `T2.*.lower_bound_only` false for all five designs; `rel_change` 2.997 % (initial), 3.899 % (13), 3.000 % (13c), 3.170 % (13d), 3.192 % (17b); `mesh_ok` true for all | pass |
| T3 | every solve converged; no "On entry to" BLAS lines | `T3.ok` true, `job_errors` [], `on_entry_to_lines` 0, `stage16_solves_converged` true | pass |

Reported, not gated (`results.json` `config.T3_interpretation`): `T3.n_unbracketed_lines` is 0 for every job
except `dense|13`, which has **40** unbracketed dense lines (above twice the first crossing; per the
registered interpretation they cannot change the minimum). The DD's first non-converged target after its
onset (`dd_first_failed_target_V`) is 230 V (initial) and 240 V (the other four), above each onset.

## R1: the artifact-budget table

| ID | Quantity | Claimed | Real | Artifact fraction | Source (file, key) |
|---|---|---|---|---|---|
| AB1 | BV gain initial to final, stage 13, old off-state, h (V) | **120.911** (planar cut 512.273 to 633.183) | **12.414** (dense true min 406.057 to 418.470; 49 lines) | **0.8973** | claimed: recomputed, equals `results/stage13/summary.json` `values/bv_initial_V`, `values/bv_final_V` (T1); real: recomputed dense minimum, theta `results/stage13/results.json` `final/theta`; `table[AB1]` |
| AB2a | BV gain over initial, stage 17 vs 17b, dense, h (V) | **90.737** (607.608 - 516.871) | **55.592** (572.464 - 516.871) | **0.3873** | `results/stage17/results.json` `R1/bv_gain_h_17`, `R1/bv_true_min_h/17`, `R1/bv_true_min_h/initial`; `results/stage17b/results.json` `R4/bv_true_min_h/17b` |
| AB2b | same, stage 17 final under dynamic reclassification (V) | 90.737 | 90.737 | 0.0 | as AB2a plus `results/stage17b/results.json` `R4/reclass_rel_change_h2` = -0.04095 (R_on,ch -4.10 % vs R0 at h/2; frozen 2.4463e-3, reclassified 2.4023e-3 ohm cm2) |
| AB3 [initial] | blocking voltage min(BV, V_pt) vs avalanche BV, 0.8 um body, h (V) | **516.871** | **94.807** | **0.8166** | claimed `results/stage13d/results.json` `dense/initial_L0/bv_min`; real: DD onset J_n = 1 A/cm2, h/2, tied body, worst column x = 1.819 um, theta = 0 |
| AB3 [13] | same, stage-13 final (V) | **181.058** | **110.551** | **0.3894** | claimed `results/stage13c/results.json` `dense/stage13_final_L0/bv_min`; DD onset, column x = 1.691 um, theta `results/stage13/results.json` `final/theta` |
| AB3 [13c] | same, 13c final (V) | **618.847** | **109.854** | **0.8225** | claimed `results/stage13c/results.json` `dense/final_L0/bv_min`; DD onset, x = 1.755 um, theta `results/stage13c/results.json` `final/theta` |
| AB3 [13d] | same, 13d final (V) | **577.940** | **83.156** | **0.8561** | claimed `results/stage13d/results.json` `dense/final_L0/bv_min`; DD onset, x = 1.819 um, theta `results/stage13d/results.json` `final/theta` |
| AB3 [17b] | same, 17b final (V) | **572.464** | **80.858** | **0.8588** | claimed `results/stage17b/results.json` `dense/final_L0/bv_min`; DD onset, x = 1.819 um, theta `results/stage17b/results.json` `final/theta` |
| AB4 [13d initial] | R_on,sp, bulk Hall vs channel mobility, h/2 (ohm cm2) | **2.7929e-4** | **2.4369e-3** | **-7.725** (channel/bulk 8.73x) | claimed: new bulk solve at level 1, equals `results/stage13d/results.json` `dense/initial_L1/R` (rel diff 0.0); real: new channel solve (d_ch = 20 nm) at level 1; `AB4_checks` reproduces stage 16 levels 0 (rel diff 0.0); recorded `results/stage16/results.json` `runs/ron\|13d_initial\|{0,2}/channel/20/R_on_sp_ohm_cm2` |
| AB4 [13d final] | same | **2.7924e-4** | **2.4333e-3** | **-7.714** (8.71x) | as above, `13d_final` |
| AB5 | BV gain over initial, stage 22 D-gradient optimum vs E-gradient optimum, 1.4 um geometry, dense, h (V) | **64.052** (608.055 - 544.003) | **64.094** (608.098 - 544.003) | **-0.00067** | `results/stage22/results.json` `p22b/dense_bv_min_D`, `p22b/dense_bv_min_E`, `initial_dense_bv_min`; BV(D) - BV(E) = -0.0427 V; SLSQP iterations 37 / 10 = 3.70x (E run from stage 21) |

Figure: `figures/stage23_artifact_budget.png` (panels: BV gains AB1/AB2a/AB2b/AB5; AB3 avalanche vs blocking
with V_pt,es(h/2) dashed; AB4 R_on,sp).

### Reading each entry

- **AB1.** Of a planar-cut gain of 120.9 V, 12.4 V remains in the dense true minimum (artifact 89.7 %).
  Note the "real" endpoints (406.06 to 418.47 V) are minima at x = 0.92 and 0.85 um (`note`); the dense
  minimum of the initial design is itself 106 V below its planar-cut value (derived).
- **AB2 has two readings; neither is preferred, and the sign of the artifact is model-dependent.** (a)
  Keeping the surface classification valid (17b), 55.6 V of the 90.7 V gain remains (artifact +38.7 %). (b)
  Under dynamic reclassification of the stage-17 final, the off-state BV is independent of the on-state
  classification, so the gain is unchanged (artifact 0.0 %), and R_on,ch is 4.10 % *lower* at h/2, i.e. the
  design is feasible and the frozen constraint was slack. `real_not_unique: true`.
- **AB4 is a model-error magnitude, not an optimizer artifact.** The channel-mobility model gives an
  R_on,sp 8.7 times the bulk-mobility value for the same two designs; nothing was optimized against the
  bulk value in this entry. The negative "fraction" (-7.7) only reflects that real exceeds claimed.
  Consistency with stage 16: recorded channel/bulk 8.02x (level 0), 8.80x (level 2) for 13d initial;
  8.01x and 8.63x for 13d final; the new level-1 ratios are 8.73x and 8.71x.
- **AB5 is on a different geometry** (d_body = 1.4 um, stage 20/21/22), not the 0.8 um geometry of AB1 to
  AB4. The simplified-adjoint artifact in the achieved optimum is -0.0427 V in dense BV (-0.00067 as a
  fraction of the gain); its cost shows up in iteration count (37 vs 10), not in the design.

### AB3 in detail (0.8 um geometry)

`AB3` entries in `results.json`; `V_pt,es` is the stage-20 electrostatic barrier criterion, reported
alongside; `blocking_V` = min(BV, DD onset) is the DD onset in every case (`limited_by`: punch-through).

| Design | Avalanche BV (V) | DD onset J_n = 1 A/cm2 (V) | V_pt,es h (V) | V_pt,es h/2 (V) | `dd_V_reached` (V) |
|---|---|---|---|---|---|
| initial | 516.871 | 94.807 | 114.074 | 117.493 | 225 |
| 13 final | 181.058 | 110.551 | 95.476 | 99.198 | 235 |
| 13c final | 618.847 | 109.854 | 111.155 | 114.490 | 235 |
| 13d final | 577.940 | 83.156 | 103.296 | 106.571 | 235 |
| 17b final | 572.464 | 80.858 | 101.227 | 104.459 | 235 |

## Key observation (post hoc, not a pre-registered prediction)

Computed after the fact from the pre-registered AB3 definitions; **it is not a pre-registered prediction, and
no prediction of this kind was registered before the data existed.** A prospective test is planned
(experiment A, see below).

Differences are in V and derived from the table above (design minus initial; initial: avalanche BV 516.871,
blocking voltage 94.807):

| Design | Avalanche BV change | Blocking-voltage change |
|---|---|---|
| 13 final | **-335.813** (181.058) | **+15.744** (110.551) |
| 13c final | +101.976 | **+15.047** (109.854) |
| 13d final | +61.069 | **-11.651** (83.156) |
| 17b final | +55.592 | **-13.949** (80.858) |

Per design: the blocking voltages of 13c (109.9 V) and 13 final (110.6 V) are **higher** than the initial
design's 94.8 V; those of 13d (83.2 V) and 17b (80.9 V) are **lower**. For the three designs whose avalanche
BV rose (13c, 13d, 17b), the avalanche-BV gains of +56 to +102 V coincided with blocking-voltage changes of
-14 to +16 V (-13.9, -11.7 and +15.0 V), i.e. the objective's gains did not transfer to the device's blocking
capability. The 13 final design is different in kind: its dense avalanche BV on this geometry and off-state
model is 181.058 V, i.e. 335.8 V *below* the initial design, and its blocking voltage is +15.7 V (the
stage-13 result was a planar-cut gain; see AB1). Across all four designs the blocking-voltage change lies
between -13.9 and +15.7 V, against avalanche-BV changes between -335.8 and +102.0 V. All five designs are
punch-through limited (`limited_by`), with blocking voltages 81 to 111 V against avalanche BVs 181 to 619 V.

## Interpretation (not a measured result)

On the 0.8 um body, the optimizations of stages 13 to 17b moved an avalanche quantity (BV of a zero-current
model) that the device does not reach, because it punches through at 81 to 111 V. The pre-registered
entries show a per-correction artifact that is large for AB1 (89.7 %) and AB3 (39 % to 86 % of the
absolute value), model-dependent in sign for AB2, a model-error magnitude (8.7x) for AB4 and negligible for
AB5. The table is not a ranking of "importance" and the fractions are not additive.

## Caveats

- 1D DD onset: the blocking voltage uses a 1D tied-body drift-diffusion onset (J_n = 1 A/cm2, level h/2) at
  **each design's own worst barrier column** (x = 1.691 to 1.819 um), not a 2D current-carrying breakdown.
- Zero-current BV: the avalanche BV is the zero-current model BV that the earlier stages optimized.
- One device family (4H-SiC DMOSFET study cell), one coefficient set; AB5 on a different geometry
  (1.4 um body); AB1 uses the old (stage-12) off-state; AB3 uses the 12c model on the 0.8 um body.
- The post hoc observation rests on four designs and one mesh; the 13 final design has a different
  avalanche behaviour (the stage-13 design was not constrained as the later ones).
- `dense|13` has 40 unbracketed dense lines (not gated, argued not to change the minimum).
- DD onset reached target voltages 225 V (initial) and 235 V (others), failing at 230 and 240 V
  respectively; every onset (81 to 111 V) is well below.

## Manifest

Code commit `a546d98d75f70f4e0ece0a6e4d36fec54faf9a61` (`a546d98`), `git_dirty: false`,
2026-10-24T00:36:00 UTC, `benchmarks/stage23/run.py`, runtime 58.02 s (`runtime_seconds`), peak RSS 265.3 MiB.
`widegap 0.6.0`. Python 3.12.3, NumPy 2.5.4, SciPy 1.18.1, torch 2.14.1+cpu (1 thread), AMD Ryzen 7 5700U,
16 logical CPUs, Linux 7.0.0-38-generic. `seeds` []; `dev: false`.

## Next: experiment A (planned, not run)

A prospective, pre-registered test: a guard-ablation optimization (with and without the punch-through guard,
body doping as a design variable) with the prediction of the blocking-voltage outcome registered before any
run. Stage 24 in `docs/PROJECT_SPEC.md` §12 (planned).

## Not claimed

- That any one artifact fraction generalizes beyond this device and model family.
- A ranking of the corrections' importance across devices.
- That the blocking-voltage observation is a tested prediction (it is post hoc), or that it holds for other
  devices, objectives or off-state models.
- A 2D punch-through voltage or agreement with measured devices.
