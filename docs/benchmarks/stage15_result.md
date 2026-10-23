# Stage 15 result: gate NOT passed (15d failed L1); 15b passed F1-F5, 15f passed V1-V2 (2026-10-22)

Raw data: `results/stage15b/`, `results/stage15d/`, `results/stage15f/` (each with
`summary.json`, `results.json`, `manifest.json`). Figures: `figures/stage15b_bv_comparison.png`,
`stage15b_field_lines.png`, `stage15d_convergence.png`, `stage15d_hardmax_positions.png`.
Pre-registration: `stage15_caveat_closure.md`, including Amendments 1 and 2 (both dated before
the official runs and disclosed below). Implemented by a subagent from the pre-registered spec,
then reviewed. Label: **new** (post-hoc re-evaluation and caveat-closure study on the stage
12c-13d devices; no optimization was run). All numbers below are quoted from the JSON files above;
anything not recorded there is marked "not recorded".

## Verdict

The stage pass rule is "F1-F5, L1, V1 and V2 hold". **L1 failed, so stage 15 as a whole does
not pass.** The three sub-stages:

| Sub-stage | Gating criteria | Result |
|---|---|---|
| 15b field-line BV | F1-F5 | **pass** (`gate_passed: true`, all five true) |
| 15d local-field adequacy | L1 | **fail** (E_ox and W2 fail; W1 and W3 pass) |
| 15f Ishikawa v_sat | V1, V2 | **pass** (`gate_passed: true`) |

The R items (R1-R4) and L2 are reported whatever the outcome. The L1 failure is not
re-registered: the criterion stands as written (see 15d for a precise description of how it fails).

Headline findings:

- The vertical-cut BV used through stages 12-13d is a good approximation of the field-line BV
  for these designs: field-line minus cut BV ranges from -6.77 V (-1.09%) to +0.29 V (+0.06%)
  over the eight design/mesh cases, and the 13d improvement survives (+60.3 V field-line vs
  +61.1 V cut on h; +61.7 V vs +62.3 V on h/2).
- Local field maxima in this device are not uniformly mesh-resolvable with the present tensor
  meshes: the E_ox and W2 area-mean measures fail L1.
- The saturation-velocity source change (Khan & Cooper 2.2e7 to Ishikawa 1.7e7 cm/s) moves R_on
  and V_T by less than 1% and I_D at V_DS = 15 V by about -12% to -20%, depending on the
  Caughey-Thomas exponent beta.

## 15b: field-line ionization-integral BV (gate PASSED, F1-F5)

Method as pre-registered, with Amendment 2 (stagnation stop rule, below). Settings from the JSON
config: N_s = 400 junction-contour seeds plus the junction crossings of the 49 dense cuts (about
458-459 seeds per evaluation), RK4 with step 0.00125 um, march step 20 V, bracket tolerance
0.02 V, e_frac 1e-3, stage-12c off-state model (phi_n = V_D w, phi_p = 0), rounding r = 0.1 um.
Level 0 is h and level 1 is h/2 (`L0`, `L1` in the JSON). The stage-12c off-state model was
used for all designs, including the 13b final design that was optimized on the old rule.

| ID | Criterion | Measured | Threshold | Result |
|---|---|---|---|---|
| F1 | 1D reduction, field-line BV vs stage-12a 1D BV (N_D = 1e16, Niwa) | h: 1787.8688 vs 1787.8755 V, rel. 3.72e-6; h/2: 1787.8754 vs 1787.8807 V, rel. 2.93e-6; max line lateral drift 5.7e-13 um | <= 1e-3 | pass |
| F2 | halving the step (0.00125 to 0.000625 um), 13d initial design, h | 517.1597 vs 517.1600 V, rel. 4.88e-7 | <= 1e-4 | pass |
| F3 | doubling N_s (400 to 800), h | initial: 517.1597 vs 517.1446 V, rel. 2.92e-5; 13d final: 577.4562 vs 577.4507 V, rel. 9.49e-6 | <= 2e-3 | pass |
| F4 | 13d initial design, h vs h/2 (denominator the h/2 value) | 517.1597 vs 516.0429 V, rel. 2.16e-3 | <= 3% | pass |
| F5 | all root-finds converged; every line ended by a stop rule | 13 root-finds, 0 not converged; 475 798 lines traced; 0 step-limit lines | 0 failures | pass |

F5 stop-rule counts over all 475 798 lines: left SiC 9713, low field 466 015, **stagnation 70**,
step limit 0 (`n_stagnation_stops: 70`). The 70 stagnation stops are spread over six of the 11
traced evaluations: 18 (initial, h), 7 (13b, h), 1 (13c, h), 5 (13d, h), 29 (initial, h, doubled
N_s), 10 (13d, h, doubled N_s); the h/2 evaluations and the halved-step evaluation had none.
The JSON also records 0 step-limit lines at the BV states.

### R1 (reported, not gated): field-line vs dense vertical-cut BV

Dense cuts: 49 cuts, all resolved (`cut_resolved_all: true` in every case). BV in V; the
difference is field-line minus cut-minimum. "Limiting line" is the field line that reaches
I_n = 1 first, given by its seed and its p-side and n-side end points (um).

| Design | Mesh | Field-line BV | Cut BV (limiting cut x, um) | Difference | Relative | Limiting line seed (um) | p end (um) | n end (um) |
|---|---|---|---|---|---|---|---|---|
| initial | h | 517.160 | 516.871 (1.700) | +0.288 | +0.056% | (1.726, -0.488) | (1.746, -0.425) | (1.424, -3.897) |
| initial | h/2 | 516.043 | 518.094 (1.771) | -2.051 | -0.396% | (1.717, -0.903) | (1.774, -0.443) | (1.480, -3.850) |
| 13b | h | 611.383 | 613.477 (1.700) | -2.094 | -0.341% | (1.716, -0.487) | (1.738, -0.424) | (1.399, -3.899) |
| 13b | h/2 | 612.329 | 616.148 (1.771) | -3.819 | -0.620% | (1.757, -0.487) | (1.769, -0.443) | (1.468, -3.852) |
| 13c | h | 614.950 | 618.847 (1.700) | -3.897 | -0.630% | (1.708, -0.487) | (1.732, -0.425) | (1.344, -3.900) |
| 13c | h/2 | 615.766 | 622.536 (1.771) | -6.770 | -1.087% | (1.513, -0.858) | (1.635, -0.421) | (1.188, -3.851) |
| 13d | h | 577.456 | 577.940 (1.700) | -0.483 | -0.084% | (1.666, -0.865) | (1.742, -0.424) | (1.438, -3.899) |
| 13d | h/2 | 577.784 | 580.413 (1.771) | -2.629 | -0.453% | (1.764, -0.487) | (1.775, -0.443) | (1.506, -3.850) |

Reading R1:

- The field-line BV is lower than the cut BV in seven of the eight cases and higher (by 0.29 V,
  +0.056%) in one. No sign was pre-registered; both signs are valid findings.
- The relative differences lie between -1.09% and +0.06%. The 13d and initial designs are within
  0.46% on both meshes; the 13b and 13c designs reach -0.62% to -1.09%, with the largest
  difference (-6.77 V, -1.09%) for the 13c final design on h/2.
- All eight limiting lines run from the p side near y = -0.42 to -0.44 um (the junction beside the
  curved body corner) to the drift region at y about -3.85 to -3.90 um, with lengths 3.42-3.50 um.
  All eight seeds are junction-contour points (kind `contour`) at x = 1.51-1.76 um; the limiting cut is x = 1.700 um on h and 1.771 um on h/2 for every design,
  which is the same location as in the stage-13d record.
- The limiting seed differs between h and h/2 and between designs in y (-0.49 or about -0.86 to
  -0.90 um) while the p-side end is at x = 1.73-1.78 um in seven of the eight cases (1.635 um for 13c on h/2). This was not analysed
  further; the limiting line is not a unique, mesh-converged object.

### R2 (reported, not gated): does the optimization gain survive?

BV gain (final minus initial design), in V, field-line vs dense cut:

| Design | h: field-line | h: cut | h/2: field-line | h/2: cut |
|---|---|---|---|---|
| 13b | +94.22 | +96.61 | +96.29 | +98.05 |
| 13c | +97.79 | +101.98 | +99.72 | +104.44 |
| 13d | +60.30 | +61.07 | +61.74 | +62.32 |

**The 13d gain survives:** field-line +60.30 V (h) and +61.74 V (h/2), i.e. 98.7% and 99.1% of
the cut gain (+61.07 V and +62.32 V). The 13b and 13c gains also survive, at 97.5% and 95.9% (h)
and 98.2% and 95.5% (h/2) of the cut gains; field-line gain ratios (final/initial) are 1.116
(13d, h), 1.120 (13d, h/2), 1.182 and 1.187 (13b), 1.189 and 1.193 (13c). The field-line method
therefore does not overturn the stage 13b/13c/13d conclusions; it removes 0.6 V (13d, h/2) to 4.7 V (13c, h/2: 104.44 to 99.72) from the cut-based gains.

### Amendment 2 (disclosed)

The development run found lines ending in a field-reversal sink (a zero of E in the p body beside
the n+ source junction, near (1.77, -0.44) um) at which |E| stays above the low-field stop
threshold, so the line orbited until the step limit (1 line in 49 586 at step 0.00125 um in the
dev run). A third stop rule (stagnation: net displacement over 200 steps below 10 steps) was added
before the official run; F5 was left unchanged. The step 0.00125 um was chosen on a termination
ground only (the config note records that 0.01, 0.005 and 0.0025 um left a few lines orbiting).
The dev BVs had been seen; F1-F4 passed in dev. As the amendment states, the rule acts only on
paths that have stopped moving and cannot change the BV of a line that reaches the avalanche
region. The official run has 0 step-limit lines and 70 stagnation stops. The JSON records 0
step-limit lines at the BV states; at which march states the 70 stagnation stops occurred, and
whether any is at a BV state, is not recorded in the summary.

## 15d: local-field mesh adequacy (gate FAILED on L1)

Measures, at V_blk = 400 V, 13d initial and final designs, on h, h/2, h/4 (levels 0, 1, 2; 1107,
4293, 16 905 nodes): E_ox, W1, W2, W3. The replacement measure is the area-weighted mean of |E|
over a fixed 0.1 um square centred at the hard-max position of the 13d initial design on h,
clipped to the oxide (E_ox) or SiC (W1-W3). The clipped area fractions are 0.250 (E_ox, centre on
the x = 0 boundary), 1.000 (W1, W2) and 0.253 (W3, centre at y = -0.0006 um on the window edge).
All eight reference values reproduce the 13d R1 smooth-max values exactly
(`consistency_with_13d_R1.max_rel_diff: 0.0`), and 0 failed evaluations.

### L1 (gating): each area-mean measure, 13d initial design

| Measure | h | h/2 | h/4 (MV/cm) | Observed order | Richardson (MV/cm) | Finest vs Richardson | Result |
|---|---|---|---|---|---|---|---|
| E_ox | 5.23856 | 5.24171 | 5.24474 | 0.056 | 5.3217 | 1.45% | **fail** (order < 0.8) |
| W1 | 2.45255 | 2.46536 | 2.47137 | 1.092 | 2.4767 | 0.21% | pass |
| W2 | 2.54957 | 2.55480 | 2.56018 | -0.044 | 2.3813 | 7.51% | **fail** (order < 0.8; 7.5% > 2%) |
| W3 | 1.96975 | 1.97264 | 1.97372 | 1.417 | 1.9744 | 0.03% | pass |

All four sequences increase monotonically with refinement. **L1 fails for E_ox and W2.** Per the
pre-registration, this is recorded as: local field maxima in this device are not mesh-resolvable
with the present tensor meshes, for the measures E_ox and W2 (the area means for W1 and W3 do
converge at order 1.09 and 1.42).

The two failures are not of the same kind, and the difference matters for how they are read:

- **W2 is a genuine non-convergence.** The area mean rises by +0.205% (h to h/2) and +0.211%
  (h/2 to h/4); the increments (5.2e-3, 5.4e-3 MV/cm) do not shrink, the order is negative, and
  the finest value is 7.5% from the Richardson estimate. The measure keeps rising at a constant
  rate per refinement.
- **E_ox's area-averaged value barely changes**: +0.060% from h to h/2, +0.058% from h/2 to h/4,
  +0.118% from h to h/4 (5.23856 to 5.24474 MV/cm). Its finest value is within 1.45% of the
  Richardson estimate, so the 2% distance criterion would pass on its own. It fails on the order
  (0.056), which is computed from two differences of 3.15e-3 and 3.03e-3 MV/cm, of order 0.06%
  of the field: the order is log2 of their ratio, so a 5% change in one difference moves it by
  about 0.07. The pre-registered order criterion is therefore ill-conditioned for this measure,
  and its Richardson value (5.32 MV/cm) is itself a large extrapolation from near-equal
  increments. This does not show that E_ox is converged (the two increments are nearly equal, not
  shrinking); it shows that the order estimate cannot discriminate here. **The failure stands as
  written; no criterion is re-registered.**

### R3 (reported): hard-max positions and pinning

Hard-max edge, initial design. "Pinned" means the argmax is within 2 cells of the same feature on
all three meshes, or on the window boundary. The figure `stage15d_hardmax_positions.png` shows
the positions.

| Measure | Hard-max position h, h/2, h/4 (um) | Nearest feature on h, h/2, h/4 (distance in cells) | Pinned (initial) | Pinned (final) |
|---|---|---|---|---|
| E_ox | (0, 0.0125), (0, 0.00625), (0, 0.003125) | p-body junction (5.5, 11.0, 22.0) | no | no |
| W1 | (1.067, -0.859), (1.067, -0.829), (1.100, -0.844) | window boundary (0.25, 0.50, 0.0); p-body junction also near | yes (p-body junction) | no |
| W2 | (1.691, -0.859), (1.691, -0.888), (1.691, -0.874) | p-body junction (0.30, 0.13, 0.72) | yes (p-body junction) | no |
| W3 | x = 0 on all three meshes, y = -0.0006 (h), -0.0003 (h/2) um, h/4 not listed here (within 0.0006 um of 0) | window boundary (0.0, 0.0, 0.0) | yes (window boundary) | yes (window boundary) |

The initial-design E_ox hard-max sits at the x = 0 symmetry line (the JFET centre) in the oxide, 5.5,
11.0 and 22.0 cells from the nearest listed feature (the p-body junction, 0.73 um away), so it is
not pinned to a listed feature; it moves towards y = 0 as the mesh is refined, with edge
midpoints at 0.0125, 0.00625 and 0.003125 um. W1 and W2 on the initial design are pinned by the
p-body junction (W1 also touches the window boundary); on the final design the JSON marks W1 and W2 as not pinned (for W2 the junction distance is
0.30, 1.14 and 2.79 cells on h, h/2, h/4). W3 is pinned by the window boundary (the SiC surface at the
JFET centre) for both designs.

### Why the 13d smooth-max values crept (reported, from `context_initial`)

13d R1 saw all window fields rise 1-6% from h to h/4. The 15d context data (initial design)
show that this was mainly the smooth-max bias, which grows with the number of edges: the
smooth-max bound ln(N)/kappa (kappa = 20 per MV/cm) is 0.165, 0.199, 0.233 MV/cm for E_ox with
27, 53, 105 edges, and the actual smooth-max minus hard-max gap grows with it.

| Measure | Edges (h, h/2, h/4) | Smooth minus hard (MV/cm) | Smooth h to h/4 | Hard h to h/4 | Hard-max order, distance |
|---|---|---|---|---|---|
| E_ox | 27, 53, 105 | 0.012, 0.035, 0.063 | +1.21% | +0.22% | 3.09, 3.2e-5 (pass) |
| W1 | 81, 342, 1332 | 0.039, 0.072, 0.128 | +6.14% | +2.41% | -0.27, 8.2% (fail) |
| W2 | 127, 495, 1950 | 0.115, 0.182, 0.246 | +5.60% | +0.70% | -0.06, 9.4% (fail) |
| W3 | 66, 240, 925 | 0.109, 0.171, 0.235 | +6.14% | +0.19% | 3.92, 8.2e-6 (pass) |

For E_ox, W2 and W3 the hard maximum moves by 0.2-0.7% while the smooth maximum moves by 1.2-6.1%,
so the creep reported in 13d is mostly the smooth-max bias. W1 is the exception: its hard maximum
also rises 2.4% (order -0.27, 8.2% from Richardson), so W1's hard-max is itself not converged
(the W1 area mean does converge). The hard-max values for W2 (2.5364, 2.5451, 2.5543 MV/cm) are
the same as the stage-12c global-peak values, which failed C1 for the same reason.
`stage15d_convergence.png` shows hard max, 13d smooth max and the area mean for each measure.

### L2 (reported): does any measure increase from the 13d initial to the final design?

Area means at the fixed centres, MV/cm, 400 V:

| Measure | h: initial to final | h/2: initial to final |
|---|---|---|
| E_ox | 5.23856 to 5.23869 (**+2.4e-5**, increases) | 5.24171 to 5.23702 (-0.090%) |
| W1 | 2.4526 to 2.3361 (-4.75%) | 2.4654 to 2.3439 (-4.93%) |
| W2 | 2.5496 to 2.4112 (-5.43%) | 2.5548 to 2.4117 (-5.60%) |
| W3 | 1.9698 to 1.9666 (-0.16%) | 1.9726 to 1.9678 (-0.24%) |

`L2_any_increase: true`, driven by E_ox on h alone, by a relative +2.4e-5 (above the 1e-6 test
tolerance, but a 1e-5-scale change; E_ox is the constrained quantity in 13d and its smooth-max was
held equal to about 1e-12 there; this area mean is a different measure and was not constrained).
On h/2 it decreases. W1-W3 decrease on both meshes, as in 13d D6.

### Manifest

15d manifest: commit `2153246`, `git_dirty: true`, 2026-10-23T00:32:43 UTC (stage 15d code
and the 15f code were committed in `2153246`; the results in `4722eb0`). The dirty flag came from
the stage-15b files, which were untracked at the time (committed afterwards in `0d63445`);
`benchmarks/stage15d/run.py` and `src/widegap/localfield.py` do not import them.
Runtime 55.3 s.

## 15f: saturation velocity from the Ishikawa thesis (gate PASSED, V1, V2)

Source: Ishikawa, Kyoto University doctoral thesis (2025, doi:10.14989/doctor.k26017), Table 3.1 and
Section 3.4.3: v_sat = 1.7e7 cm/s for both //c and perpendicular-to-c directions (isotropic),
recommended after averaging over N_D below 1e16 cm-3. The Khan & Cooper entries (2.2e7 cm/s, beta
1.2) remain in the registry unchanged (`khan_cooper_entries_unchanged` in the JSON).

| ID | Criterion | Result |
|---|---|---|
| V1 | the 10 new registry entries (`vsat_n_par_c`, `vsat_n_perp_c_ishikawa`, and eight `ct_beta_*_ishikawa_*` entries at 2e15, 5e15, 1e16, 3e16 cm-3 for //c and perpendicular) pass `check_entry` | **pass**: `check_entry: "ok"` and `usable_in_registry: true` for all 10; the verifier names are not in the summary; they are in the registry entries themselves (`extracted_by: literature-agent-15f (sonnet)`, `verified_by: verifier-agent-15f (haiku)`, in `sic4h.json`) |
| V2 | all DD solves converge for every R4 evaluation, both Ishikawa variants | **pass**: `V2_errors` null for all six case/variant pairs, `converged_all: true` for all nine runs |

### Amendment 1 (disclosed)

The thesis fit exponents at N_D <= 1e16 are below 1 (0.93 perpendicular, 0.87 //c at 1e16), and the
stage-9a Caughey-Thomas implementation requires beta >= 1 (`ct_factor`; F^(beta-1) is singular at
zero-field edges). The official evaluation therefore uses the Ishikawa v_sat with two model-choice
values, **beta = 1.0** (the model's lower bound) and **beta = 1.2** (Khan & Cooper), recorded as
model choices, not registry values (`beta_note`). **beta < 1 is unsupported by the model** and was
not run. A dev run at beta = 1.0 was seen before the amendment (R_on about +0.7%, I_D at 15 V
about -20%); the official results at beta = 1.0 are consistent with it (below), and the
amendment did not choose beta from that outcome: both bracketing values were run.
Also recorded by the amendment: I_D at V_DS = 15 V uses V_GS = 15 V (the stage-12 on-state gate
bias; the spec did not state this).

### R4 (reported): Khan & Cooper vs Ishikawa

Khan & Cooper: v_sat 2.2e7 cm/s, beta 1.2. Ishikawa: 1.7e7 cm/s at beta 1.0 and 1.2. Ratios are
Ishikawa over Khan & Cooper; V_T shift in volts (Ishikawa minus Khan & Cooper). Mesh level 0.

| Case | Variant | R_on,sp (ohm cm2) | R_on ratio | V_T (V) | V_T shift (V) | I_D at V_DS = 15 V (A/cm) | I_D ratio |
|---|---|---|---|---|---|---|---|
| 12b | Khan & Cooper | 2.9190e-4 | | 7.8551 | | 3.6180 | |
| 12b | Ishikawa, beta 1.0 | 2.9401e-4 | 1.0072 | 7.8571 | +0.0021 | 2.8931 | 0.7996 |
| 12b | Ishikawa, beta 1.2 | 2.9215e-4 | 1.0009 | 7.8555 | +0.0005 | 3.1955 | 0.8832 |
| 13d initial | Khan & Cooper | 2.8555e-4 | | 7.6661 | | 4.4451 | |
| 13d initial | Ishikawa, beta 1.0 | 2.8753e-4 | 1.0069 | 7.6754 | +0.0092 | 3.5280 | 0.7937 |
| 13d initial | Ishikawa, beta 1.2 | 2.8579e-4 | 1.0008 | 7.6684 | +0.0023 | 3.8733 | 0.8714 |
| 13d final | Khan & Cooper | 2.8555e-4 | | 7.6644 | | 4.4672 | |
| 13d final | Ishikawa, beta 1.0 | 2.8742e-4 | 1.0065 | 7.6737 | +0.0093 | 3.5445 | 0.7934 |
| 13d final | Ishikawa, beta 1.2 | 2.8577e-4 | 1.0008 | 7.6667 | +0.0023 | 3.8897 | 0.8707 |

(The 13d Khan & Cooper R_on values equal the stage 13d summary values: 2.855541e-4 for both designs,
`consistency_13d`.)

- **R_on changes by +0.65% to +0.72% (beta 1.0) and +0.08% to +0.09% (beta 1.2); V_T moves by at
  most +0.0093 V (about 0.12%). Both are below 1%.** Lower v_sat raises R_on slightly; at
  V_DS = 0.1 V the channel is not velocity-saturated, so little effect was expected.
- **I_D at V_DS = 15 V falls by 20.0% to 20.7% (beta 1.0) and 11.7% to 12.9% (beta 1.2).** The
  v_sat effect alone (beta 1.2 rows) is -12% to -13%; the extra decrease at beta 1.0 is the beta
  effect.
- The 13d optimization changes R_on_sp by 0.0% (fixed in the optimization) and V_T by -0.0018 V
  (Khan & Cooper), and I_D at 15 V by +0.5%, so the initial-versus-final comparison does not
  depend on the v_sat source at this resolution.
- **The anisotropy caveat is closed** by the thesis finding that v_sat is isotropic (1.7e7 in
  both directions): applying a single perpendicular-to-c value to every edge is not an error in
  direction. **The v_sat value itself remains uncertain**: 2.2e7 (Khan & Cooper) vs 1.7e7
  (Ishikawa); the thesis attributes the higher value to measurement difficulties on heavily doped
  samples and notes the Hall scattering factor may make velocities derived from Hall mobility
  overestimates. The R4 results are a sensitivity bound within these two published values, not a
  determination of the true value.
- **Holes use the electron v_sat and beta** (`sic.py`, `vsat_p = vsat_n`); there is no registered
  hole v_sat. R4 does not separate a hole effect (holes are minority carriers in this n-channel
  device), so the effect of this assumption is not quantified here.
- Note that R_on, V_T and I_D here are from the stage-12/13d on-state models with Hall mobility
  (caveat unchanged, 15e).

### Manifest

15f manifest: commit `2153246`, `git_dirty: true`, 2026-10-23T00:33:41 UTC, same dirty-tree note as
15d (untracked 15b files, committed afterwards; not imported by 15f). `dev_unverified_registry:
false`. Runtime 80.8 s.

## 15a restated (no new run)

The zero-current off-state model (phi_n = V_D w) is not replaced by a full drift-diffusion
off-state. The pre-registration restates the existing stage-12c gates as the bound on the
interpolation error of the electron quasi-Fermi potential, and no new run was made. The recorded
values (`results/stage12c/summary.json`; [12c/13c result](stage12c_13c_result.md)):

- **C2 (insensitivity to the w interpolation): relative difference 0**, bit-identical for peak
  field (2.5364 MV/cm), rounded planar BV (566.26 V) and the dense minimum (516.87 V at
  x = 1.70 um); n < 1e-50 where the two w fields differ.
- **C3 (carrier self-consistency): worst n/N = 1.96e-19** (limit 1e-6; worst case
  `harmonic_march_to_BV_L0`), i.e. at most 2e-19.

These bound the electron-density error that the choice of w can introduce in the zero-current
model. They do not make the model a full DD solution: the hole and current response, and any
difference of a full-DD off-state at hundreds of volts, remain unquantified. A full-DD off-state
at that bias would meet the same hole-row underflow as 15c.

## Items not run in this stage

- **15c (V_DS above 15 V in the on-state model): deferred.** Needed only for output curves. A future
  scheme must be gated on a 1D reverse-biased diode up to at least 1 kV against the existing
  log-space solution.
- **15e (Hall, not channel, mobility): blocked on sources** (paper request sent to the owner).
- **Oxide-field guideline (3 MV/cm): blocked on a source.** It remains unsourced; the 13d E_ox of
  5.25 MV/cm (about 1.75 times the guideline) stays a reported number with no reliability claim.

## Provenance

- 15b: `benchmarks/stage15b/run.py`, manifest commit `0d63445926cd76b354c4faf9b480ca8d0364b1f2`,
  `git_dirty: false`, 2026-10-23T00:53:41 UTC, runtime 384.2 s; results committed in `e09a644`.
- 15d: `benchmarks/stage15d/run.py`, commit `2153246153aff157fd67545094e1af803937538d`,
  `git_dirty: true` (see above), 2026-10-23T00:32:43 UTC, 55.3 s.
- 15f: `benchmarks/stage15f/run.py`, commit `2153246153aff157fd67545094e1af803937538d`,
  `git_dirty: true` (see above), 2026-10-23T00:33:41 UTC, 80.8 s.
- All three: Python 3.12.3, widegap 0.2.0, NumPy 2.5.4, SciPy 1.18.1, torch 2.14.1+cpu (1 torch
  thread), AMD Ryzen 7 5700U, 16 logical CPUs, Linux 7.0.0-38-generic. The deterministic
  pipelines record no seeds (`seeds: []`).

## Not claimed

- A full 2D avalanche simulation: no impact-ionization generation is solved in drift-diffusion;
  the field-line BV is a post-hoc ionization-integral evaluation of the zero-current off-state
  fields, and a different approximation from the vertical cuts, not an exact reference.
- Agreement with measured devices (none claimed for any stage 12-15 result).
- Mesh convergence beyond h/4 (15d) or h/2 (15b); mesh-converged local fields (15d L1 failed).
- That the E_ox area mean has converged: its order estimate is ill-conditioned, not shown good.
- A hole saturation velocity, or any effect of the hole v_sat assumption.
- Any beta below 1 (outside the model's domain).
- A determination of the true saturation velocity: 2.2e7 vs 1.7e7 cm/s remains an uncertainty.
- Channel (interface) mobility, V_DS above 15 V, or an oxide reliability limit.
- A full-DD off-state result (15a is a restated bound, not a new run).
