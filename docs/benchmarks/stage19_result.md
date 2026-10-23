# Stage 19 result: gate NOT passed (P1 and P2 failed, P3 passed); R2 verdict: punch-through CONFIRMED (2026-10-23)

Raw data: `results/stage19/` (`summary.json`, `results.json`, `manifest.json`). Figures:
`figures/stage19_curves.png`, `figures/stage19_doping.png`. Pre-registration:
`stage19_punchthrough_check.md`. All numbers below are quoted from the JSON files above; anything not
recorded there is marked "derived". Mesh levels: L0 = h (stage-13d level 0), L1 = h/2, L2 = h/4
(supplementary, floating body, low bias only, not gated). Columns: x = 2.0 um (mesh column 1.998 um) and
x = 1.6 um. Device: stage-13d initial design at theta = 0, r = 0.1 um. Method: 1D drift-diffusion (stage 9a
edge core, incomplete ionization, no impact ionization, R = 0) through the n+/p/n-/n+ stack, one
`continue_bias` march per run (1 V steps to 150 V, 5 V steps to 400 V).

## Verdict

`criteria`: **P1 false, P2 false, P3 true**; `gate_passed: false`, `dev: false`. **Stage 19 does not pass.**

R2 verdict (`R2.verdict`): **confirmed (DD V_pt below the avalanche BV)**, for both body configurations
(`verdict_by_body`: floating and tied both "confirmed"). Both failed criteria are discussed below: P1a is a
badly designed criterion (finding, not re-registered), P2 fails on numerics. The physical conclusion
(punch-through far below the avalanche BV) rests on the practical-onset currents of R1, not on the failed
criteria and not on the pre-registered V_pt definition alone (see R1).

## Criteria

| ID | Criterion | Measured | Result |
|---|---|---|---|
| P1 | P1a and P1b (floating, gated; tied reported) | P1a fails, P1b passes (`P1a_ok_floating: false`, `P1b_ok_floating: true`; `tied_P1a_ok: false`, `tied_P1b_ok: true`) | **fail** |
| P2 | every continuation step to 150 V converged, no "On entry to" lines | `ok: false`; `on_entry_to_lines: 0`; floating runs stop at 73 to 101 V | **fail** |
| P3 | V_pt mesh change <= 5 % (floating) | 3.04 % (x = 2.0), 0.042 % (x = 1.6); `ok: true` | pass |

### P1a: fails, and the criterion is ill-designed

P1a requires |ln J_h - ln J_h/2| / |ln J_h/2| <= 1e-3 at V_D = 10 V (`P1a_rel_log` tolerance 0.001).

| Config | x (um) | J_n(10 V), h (A/cm2) | J_n(10 V), h/2 (A/cm2) | abs dlnJ | rel_log |
|---|---|---|---|---|---|
| floating | 2.0 | -4.6816e-45 | -3.6917e-45 | 0.2376 | 0.002322 |
| floating | 1.6 | -2.8824e-45 | -2.3070e-45 | 0.2227 | 0.002166 |
| tied | 2.0 | -4.6816e-45 | -3.6559e-45 | 0.2473 | 0.002417 |
| tied | 1.6 | -2.6556e-45 | -2.2826e-45 | 0.1513 | 0.001472 |

All four rel_log values (0.0015 to 0.0024) exceed 1e-3. **Finding (not a re-registration, the gate is not
re-defined):** at V_D = 10 V the "current" is a leakage of about 1e-45 A/cm2, an equilibrium-scale quantity
that has no physical meaning; ln|J| is about -102, so a 0.2 to 0.25 change in ln J (a 15 to 28 % change in J,
`rel_J` 0.249 to 0.281) is divided by about 102 and compared with 1e-3. Comparing log-currents of a
meaningless 1e-45 A/cm2 leakage with a 1e-3 tolerance was an ill-designed criterion: it tests a quantity
that is not the punch-through current. The supplementary h/4 run (floating) gives `abs_dlnJ_h2_h4` 0.109
(x = 2.0) and 0.101 (x = 1.6), i.e. J_n at 10 V is not mesh-converged at that level either. The mesh
dependence of the punch-through onset voltage itself is quantified in R1 (practical-onset voltages) and P3.
No new criterion is registered here.

### P1b: passes

|barrier_1D(10 V) - barrier_12c(10 V)| <= 1 V_T for every (column, mesh); recorded `diff_10V` (floating):
-0.4216 (x = 2.0, L0), -0.2640 (x = 2.0, L1), +0.1873 (x = 1.6, L0), +0.2006 (x = 1.6, L1) V_T. All
`ok: true`; the tied configuration has |diff_10V| <= 0.4216 as well (`tied_P1b_ok: true`). At 40 V the
difference grows to -2.019 V_T (x = 2.0, L0) in the tabulated comparison (informational, not gated).

### P2: fails (numerics)

`ok: false`, `on_entry_to_lines: 0`. Highest converged V_D (`V_last_converged`) and failure data
(`reached_150V`):

| Config | x (um) | Level | V_last_converged (V) | reached_150V | failed at target (V) | residual_test_ok | step_test_ok |
|---|---|---|---|---|---|---|---|
| floating | 2.0 | L0 | 73 | false | 74 | true | false |
| floating | 2.0 | L1 | 76 | false | 77 | true | false |
| floating | 1.6 | L0 | 94 | false | 95 | true | false |
| floating | 1.6 | L1 | 101 | false | 102 | true | false |
| tied | 2.0 | L0 | 108 | false | 109 | false | true |
| tied | 2.0 | L1 | 225 | true | 230 | false | true |
| tied | 1.6 | L0 | 135 | false | 136 | false | true |
| tied | 1.6 | L1 | 200 | true | 205 | false | true |

**Floating body (gated): stalls at 73 to 101 V.** At the failing step the last residual ratio is 5.6e-16 to
1.5e-15 (`residual_test_ok: true`) but the last step `last_dX` is 3.1e-8 to 3.6e-7 against the 1e-8 limit
(`step_test_ok: false`), 200 iterations. The worst row is a phi_n or phi_p row in the depleted body (for
example x = 2.0, L0: `phi_n`, node 16, depth 0.573 um, F = 2.9e-37, S = 5.1e-20; x = 1.6, L1: `phi_n`, node
30, depth 0.527 um; the hole density at the worst nodes is of order 1e-27 to 1e-50 in the recorded rows).
This is a double-precision floor of the same class as the stage-15c / hole-row underflow caveat in the
README (carrier rows with densities far below the double-precision resolution of the Newton step): the
failure is in rows that carry no physical current, not in a divergence. It is the cause recorded by the
diagnostic fields (derived from the `failure` records); no separate run isolates it.

**Tied body (supplementary): reaches 108 to 225 V** (L0: 108 V at x = 2.0 and 135 V at x = 1.6; L1: 225 V at
x = 2.0 and 200 V at x = 1.6). These runs fail differently: `residual_test_ok: false` (ratio 2.2e-14 to
4.4e-14) with `step_test_ok: true`. The two L1 tied runs pass 150 V (`reached_150V: true`); the two L0 tied
runs stop below 150 V. The tied runs get past the floating stall voltages (73 to 101 V); the cause is not isolated by a run
(plausibly the hole contact keeps the body hole row away from the precision floor). They do not make P2
pass, because P2 is gated on the floating configuration.

### P3: passes

V_pt (smallest V_D with |J_n| > 1e3 |J_n(10 V)|, floating):

| x (um) | V_pt, h (V) | V_pt, h/2 (V) | rel_change | ok |
|---|---|---|---|---|
| 2.0 | 57.145 | 58.936 | 0.03039 | true |
| 1.6 | 79.211 | 79.177 | 0.000420 | true |

Tied: V_pt 57.145 to 58.950 V (x = 2.0), 79.191 to 80.712 V (x = 1.6), rel_change 0.0306 and 0.0188.
**Caveat on what this V_pt is (recorded `note_V_pt_definition`):** 1e3 x J_n(10 V) marks a barrier drop of
about ln(1e3) = 7 V_T from a ~1e-45 A/cm2 level (`barrier_drop_at_V_pt_VT` 5.0 to 6.4 V_T; barrier at V_pt
108.8 to 110.3 V_T). It is not an onset of appreciable current; P3 therefore says the mesh dependence of a
tiny-leakage threshold is small, not that the practical onset voltage is mesh-converged. The practical onset
voltages in R1 differ by 2 to 3 % between L0 and L1 at the same column (derived).

## R1: barrier curves, practical onset and comparison with the avalanche BV

### Barrier curves (h/2, L1)

Barrier = psi(surface) - min psi over the barrier window, in V_T (`barrier_VT`); J_n = |J_n| at the barrier in
A/cm2. Floating runs stop before the curve collapses (last converged point shown); tied runs continue.

| V_D (V) | x = 2.0, tied: barrier | J_n | x = 1.6, tied: barrier | J_n |
|---|---|---|---|---|
| 10 | 115.64 | 3.656e-45 | 115.41 | 2.283e-45 |
| 40 | 114.08 | 3.868e-44 | 114.39 | 1.043e-44 |
| 60 | 108.06 | 2.831e-41 | 113.44 | 4.049e-44 |
| 80 | 43.23 | 4.359e-13 | 110.70 | 1.427e-42 |
| 90 | 21.63 | 7.586e-04 | 82.10 | 1.882e-30 |
| 100 | 7.46 | 5.847e+02 | 54.85 | 2.371e-18 |
| 120 | 2.87 | 1.601e+04 | 13.83 | 6.702e-01 |
| 130 | 1.84 | 2.387e+04 | 5.71 | 2.097e+03 |
| 150 | 0.62 | 3.761e+04 | 2.43 | 1.544e+04 |
| 200 | 0.00 | 6.278e+04 | 0.00 | 4.349e+04 |

Floating, h/2, last converged point: x = 2.0, 76 V: barrier 53.83 V_T, J_n 9.674e-18 A/cm2; x = 1.6, 101 V:
barrier 52.30 V_T, J_n 2.842e-17 A/cm2. At x = 1.6 the floating and tied curves coincide at 90 V and 100 V
(barrier 82.10 and 54.85 V_T, J_n 1.882e-30 and 2.371e-18 A/cm2 in both), so over the range where both
exist the body configuration makes no visible difference (derived from the series). The zero-current stage-12c
barrier is 115.9 V_T at 10 V (x = 2.0, L1) and the stage-19 DD barrier is within 0.5 V_T of it (P1b).
Figure: `figures/stage19_curves.png`; doping of the two columns: `figures/stage19_doping.png`.

### Practical onset (tied runs; the only runs that reach it)

Voltages at which |J_n| reaches 1e-3, 1 and 1e3 A/cm2 (`onset_V`), and at which the barrier falls below
10 V_T (`V_barrier_below_10VT`, the stage-17c scoping metric):

| x (um) | Mesh | J_n = 1e-3 (V) | J_n = 1 (V) | J_n = 1e3 (V) | barrier < 10 V_T (V) | V_reached (V) |
|---|---|---|---|---|---|---|
| 2.0 | L0 | 86.75 | 91.52 | 98.30 | 95.56 | 108 |
| 2.0 | L1 | 90.16 | 94.31 | 100.99 | 97.33 | 225 |
| 1.6 | L0 | 112.37 | 117.80 | 125.93 | 121.98 | 135 |
| 1.6 | L1 | 115.63 | 120.31 | 127.85 | 123.46 | 200 |

**Worst column: x = 2.0 um**, in all four of these metrics and in the pre-registered V_pt (58.94 V floating,
58.95 V tied, h/2). Note that the scoping had x = 1.6 um as the worst column in the 2D zero-current barrier
metric; the 1D DD stack finds x = 2.0 um earlier. Over both columns and both meshes the practical onset
(J_n from 1e-3 to 1e3 A/cm2) spans **86.7 to 127.9 V**, and the 10 V_T crossing 95.6 to 123.5 V. Floating
runs have no practical-onset value (`onset_V` all `None`) because they stop first. Comparison with the
scoping estimate: 100 to 120 V for the 10 V_T crossing (stage-17c scoping, 2D zero-current model) and
150 V (1D analytic estimate, `analytic_1D_estimate_V`) against 95.6 to 123.5 V here. The DD result is
consistent with the scoping range and lower than the analytic estimate.

### Against the avalanche BV

Avalanche BV of the same device: `avalanche_dense_min_BV_V` 516.87 V (stage 12c, dense minimum, level 0);
stage-13d initial design `bv_obj` 525.08 V, dense minimum 516.87 V (L0) and 518.09 V (L1), limiting cut
x = 1.6 um. **Key statement: in this model the study DMOSFET body punches through at about 90 to 130 V
(practical onset, tied runs; 86.7 to 127.9 V in the table above), far below the 512 to 620 V avalanche BVs
reported since stage 12. Those BVs are avalanche limits of an electrostatic zero-current model that cannot
carry source-injected electrons; the model could not see punch-through.** The factor between the worst-column
practical onset (J_n = 1e3 A/cm2 at 98.3 V, L0) and the 516.87 V avalanche BV is 5.3 (derived).

## R2: verdict

`verdict: "confirmed (DD V_pt below the avalanche BV)"`; floating and tied both confirmed. By the registered
definition (finest-mesh floating V_pt of the worst column, 58.94 V, against 516.87 V) and, independently of
that definition, by the practical-onset currents of the tied runs (86.7 to 127.9 V). The floating-body
practical onset is not observed (the runs stall first); the floating barrier is already at 52 to 66 V_T at
the stall points and the floating and tied curves coincide where both exist (x = 1.6), so no difference
between the configurations is seen, but a floating-body onset voltage is not measured.

## Why it punches through (derived from the doping figure, not a separate run)

The column doping (`figures/stage19_doping.png`) has a p-body of N_A about 2.7e17 to 3e17 cm-3 under the n+
source that ends at about 0.8 to 0.9 um depth with a drift doping of about 5e16 below (the spec's 1D estimate
used about 0.3 um of body under the source at 2.7e17 or more and drift doping 5e16). The body under the
source is too thin: the drain-side depletion reaches the source-side barrier region before any avalanche
field is reached.

## Caveats

- 1D column model with two columns; no 2D punch-through voltage is computed (stated in the pre-registration).
- The gated floating runs stall at 73 to 101 V and do not reach the 150 V minimum of the spec; the practical
  onset results come from the tied supplementary runs, whose body contact is a hole-only ohmic contact at the
  node nearest mid-body depth (0.55 um, `body_tie_node_depth_um`).
- Tied runs end on a residual-test floor (`residual_test_ok: false`, ratio 2e-14 to 4e-14) after the
  reported range; the curves up to V_last_converged passed the unchanged stage-8.1 convergence test.
- Non-finite values appear in rejected Newton trial steps (`nonfinite_solves` 40 to 145 per run); 0 "On entry
  to" lines (`on_entry_to_lines: 0`).
- No impact ionization, no recombination, isothermal 300 K, vertical electron mobility with
  Caughey-Thomas saturation; no 2D effects (JFET-region pinch-off, lateral spreading).
- The J_n values below about 1e-30 A/cm2 are not physically meaningful; only the voltages at 1e-3 A/cm2 and
  above are used for the statements above.

## Manifest

Code commit `d568cc0cd0f7f520f4fbb6ff23a216aef6500404` (`d568cc0`), `git_dirty: false`,
2026-10-23T19:04:46 UTC, `benchmarks/stage19/run.py`, runtime 124.1 s (`runtime_seconds`
124.1043361639604; `parallel_seconds` 124.07), 10 workers, peak RSS 263.1 MiB. Python 3.12.3, NumPy 2.5.4,
SciPy 1.18.1, torch 2.14.1+cpu (1 thread), AMD Ryzen 7 5700U, 16 logical CPUs, Linux 7.0.0-38-generic.
`widegap 0.5.0` (installed metadata 0.5.0). `seeds: []`; `dev: false`.

## Not claimed

- A 2D punch-through voltage (1D is a column model; the true 2D V_pt may differ).
- Agreement with measured devices.
- A floating-body punch-through onset voltage (floating runs stall at 73 to 101 V).
- That P1a or P2 passed, or that the P1a criterion is re-registered; P1a is reported as ill-designed only.
- That the avalanche BV results of stages 12 to 17b are wrong as avalanche objectives; they are not the
  device blocking voltage in this model.
- A design fix: a body redesign against punch-through is the next stage (planned as stage 20) and is not
  specified here.
