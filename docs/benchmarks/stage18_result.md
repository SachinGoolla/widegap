# Stage 18 result: gate NOT passed (A1-A3 hold; A4 mesh and A5 d_ch fail) (2026-10-23)

Raw data: `results/stage18/` (`summary.json`, `results.json`, `manifest.json`). Figures:
`figures/stage18_bracket.png`, `figures/stage18_r2.png`. Pre-registration:
`stage18_accumulation_bracket.md`. Comparison data: `results/stage17/` ([stage 17 result](stage17_result.md)),
`results/stage17b/` ([stage 17b result](stage17b_result.md)) and `results/stage13d/`. All numbers below are
quoted from the JSON files above; anything not recorded there is marked "not recorded" or "derived"
(computed from recorded values). Level 1 = h/2, level 2 = h/4 (`config.level_names`). The L end is
mu_acc = 17.5 cm2/Vs (Zhao 2007), the U end is the stage-16 default (bulk Hall mobility).

## Verdict

The stage pass rule is "A1-A5 hold". A1, A2 and A3 hold; **A4 (mesh) and A5 (d_ch) fail**
(`gate_passed: false`, `dev: false`), so **stage 18 does not pass**. R1 and R2 are reported as pre-registered.

The bracket numbers below are therefore *not* mesh- or d_ch-converged for the L end and must be quoted with
that sensitivity as uncertainty. The failure is in the L end only: A2 shows the U end reproduces stage 16
exactly.

## Criteria

| ID | Criterion | Measured | Threshold | Result |
|---|---|---|---|---|
| A1 | registry entry for Zhao's mu_s passes `check_entry` incl. independent re-verification | `check_entry: "ok"`, `verified_by: "verifier-agent-18 (haiku)"`, `verified_value_match: true`; value 17.5 cm^2/(V s), doi 10.1016/j.sse.2007.06.006, locator p. 1141 Sec. 3 | ok | pass |
| A2 | (U) setting reproduces stage-16 channel R_on (12b, h/2) | `R_on_U_this_run` 0.003215047237922234 = `R_on_stage16` 0.003215047237922234; `rel_diff` 0.0 | <= 1e-12 | pass |
| A3 | every solve converges; no "On entry to" BLAS lines | `failed: []`, `on_entry_to_lines` 0, `stderr_bytes` 0 | none | pass |
| A4 | R_on(L) h/2 to h/4 (12b) | 0.003338162994539712 to 0.0034287652781895077 ohm cm2, `rel_diff` 0.027141360022861605 (2.71 %) | <= 2 % | **fail** |
| A5 | R_on(L) d_ch 20 vs 40 nm (12b, h/2) | 0.003338162994539712 to 0.0033924496701982615 ohm cm2, `rel_diff` 0.01626244007477995 (1.63 %) | <= 1 % | **fail** |

## R1: bracket [R_on(U), R_on(L)] (ohm cm2)

`complete_results_file` for the 17b design: **true** (`designs.s17b_final`: `present: true`,
`complete_results_file: true`, `dev: false`, source `results/stage17b/results.json`), so the 17b rows come
from the official, complete stage-17b run. The stage-17 final is from `results/stage17/results.json`
(`present: true`). "Bulk" is the stage-12/13 bulk-mobility R_on (no channel model), shown for reference.

| Design | Level | R_on(bulk) | R_on(U) | R_on(L) | L/U - 1 | acc. edges / columns |
|---|---|---|---|---|---|---|
| 12b | 1 (h/2) | 2.8526e-4 | 3.215047e-3 | 3.338163e-3 | +3.83 % | 110 / 12 |
| 12b | 2 (h/4) | 2.8419e-4 | 3.255829e-3 | 3.428765e-3 | +5.31 % | 460 / 24 |
| 13d initial | 1 | 2.7929e-4 | 2.436882e-3 | 2.539322e-3 | +4.20 % | 100 / 11 |
| 13d initial | 2 | 2.7803e-4 | 2.446694e-3 | 2.598883e-3 | +6.22 % | 420 / 22 |
| 13d final | 1 | 2.7924e-4 | 2.433282e-3 | 2.523156e-3 | +3.69 % | 100 / 11 |
| 13d final | 2 | 2.7792e-4 | 2.398423e-3 | 2.608255e-3 | +8.75 % | 440 / 23 |
| stage-17 final | 1 | 2.8730e-4 | 2.337085e-3 | 2.550439e-3 | +9.13 % | 110 / 12 |
| stage-17 final | 2 | 2.8594e-4 | 2.402331e-3 | 2.597068e-3 | +8.11 % | 440 / 23 |
| 17b final | 1 | 2.7971e-4 | 2.436882e-3 | 2.537635e-3 | +4.13 % | 100 / 11 |
| 17b final | 2 | 2.7837e-4 | 2.446599e-3 | 2.597122e-3 | +6.15 % | 420 / 22 |

Figure: `figures/stage18_bracket.png`. Mesh sizes: 4293 nodes (level 1), 16905 nodes (level 2).

Observations (derived from the table):

- **Bracket width, L/U - 1: +3.7 % to +9.1 % over all rows** (level 1: +3.69 % to +9.13 %; level 2:
  +5.31 % to +8.75 %). The exploratory scoping bracket (20-50 cm2/Vs) was +2.0 % to +6.2 %.
- **Mesh sensitivity of the width** is large and not monotonic across designs: 12b 3.83 to 5.31 %, 13d
  initial 4.20 to 6.22 %, 13d final 3.69 to 8.75 %, 17b 4.13 to 6.15 %, stage-17 final 9.13 to 8.11 %. The
  13d-final level-2 value (8.75 %) is out of line with the 13d-initial (6.22 %) and 17b (6.15 %) rows; the
  cause is not analysed here.
- **d_ch sensitivity (12b, h/2):** R_on(L) at d_ch = 40 nm is 3.392450e-3; against R_on(U) = 3.215047e-3 that
  is L/U - 1 = +5.52 % (derived), against +3.83 % at 20 nm.
- The U values on the 17b and 13d-initial designs equal the stage-17b/16 R_on,ch values (e.g. 2.436882e-3 at
  h/2, 2.446694e-3 at h/4) as expected from A2.
- The U column is 7.1 to 10.5 times the bulk column (`U_over_bulk_minus_1`: 7.13 to 10.46), the stage-16
  channel-model effect and unrelated to the accumulation mobility.
- Ordering of designs under L is not resolved: the spread between designs at one level (e.g. level 2: L
  2.5971e-3 to 2.6083e-3 across the four designs other than 12b) is comparable to or smaller than the mesh
  and d_ch sensitivity above.
- Not recorded: the edge counts for the d_ch = 40 nm run.

## R2: JFET-level sensitivity dR_on/dtheta_J (FD, 13d initial design, h/2)

FD at eps = 1e-3 and 5e-4 (the two agree; `spread` <= 7.9e-11 in all three), 5 solves each,
all converged.

| Model | R_on at theta = 0 | dR_on/dtheta_J (ohm cm2) | |dR/dtheta_J| / R_on (derived) |
|---|---|---|---|
| bulk (no channel model) | 2.792927e-4 | -4.336342e-5 | 0.155 |
| U (channel model, stage 16) | 2.436882e-3 | -9.426121e-5 | 0.0387 |
| L (channel model + mu_acc 17.5) | 2.539322e-3 | -2.819609e-4 | 0.111 |

`R2.ratios` (absolute values): **U/bulk = 2.1737, L/bulk = 6.5023, L/U = 2.9913**. Figure:
`figures/stage18_r2.png`.

Conclusion (R2 question: does the accumulation mobility explain the stage-17 factor 2.2-2.4?):

- **No.** The factor already appears with the U setting (no accumulation mobility change): U/bulk = 2.17,
  essentially the stage-17 value (stage 17 R2 reports 2.17 to 2.53 times the bulk value). It comes from the
  stage-16 channel model (the Noguchi mu(N_S) mobility), not from the accumulation mobility. The U/bulk
  value 2.17 is just below the "2.2-2.4" range quoted in the spec; the stage-17 record is 2.17-2.53.
- Applying the pessimistic accumulation mobility (L) raises the sensitivity a **further factor of 2.99**
  (about 3x; verified, L/U = 2.9913), so L/bulk is 6.50.
- Caveat: the L number inherits the A4/A5 failures (it is a mesh- and d_ch-sensitive model-definition quantity,
  see Diagnosis).

## Diagnosis (interpretation, not a gate)

This section is interpretation. It is not tested by a pre-registered criterion.

The A4 and A5 failures plausibly stem from the model definition. The constant Zhao mu_s is applied to **all**
electrons on lateral edges within d_ch of the interface in the gated n-type columns, including the
equilibrium (bulk, N_D-donated) electrons. Zhao's mu_s, however, is a mobility of the accumulated
(gate-induced) charge in a lumped I_ACC + I_CH model. The conductance replaced by the L setting therefore
grows with d_ch (more layer volume) and with the number of nodes falling inside the layer. The recorded
counts are consistent with this: accumulation edges go 110 to 460 (12b), 100 to 420 (13d initial and 17b), 100
to 440 (13d final), 110 to 440 (stage-17 final) from h/2 to h/4, i.e. about 4x, while the column count only
doubles (12 to 24 for 12b). The L-vs-U difference grows from h/2 to h/4 for 12b (+3.83 % to +5.31 %) and
from 20 to 40 nm (+3.83 % to +5.52 %), the direction this explanation predicts. This is a hypothesis; no
run isolates it.

A follow-up (**18b, not pre-registered yet**) could apply mu_acc only to the excess charge n - n_eq. It would
need its own pre-registration and sourcing check.

Despite the failures, the bracket width (L/U - 1, +3.7 % to +9.1 %) is the best available estimate of the size
of the effect, with its mesh and d_ch sensitivity (12b: +3.83 to +5.31 % across one refinement, +3.83 to
+5.52 % across d_ch 20 to 40 nm) as uncertainty. It is not a best-estimate mobility effect and no
accumulation mobility for NO-nitrided oxide is claimed.

## Caveats

- Zhao's oxide is not nitrided and its gate overdrive (about 3-7 V above flat band) is below our 15 V; L is a
  pessimistic, sourced lower bound, not a best estimate.
- The L model is mesh- and d_ch-dependent (A4, A5 failed); L values are not converged.
- All stage-16 caveats carry over (Hall factor 1, NO-annealed parameters, no N_it dependence, 300 K; channel
  R_on not mesh-converged, H5 failed).
- R2 is evaluated at one design (13d initial, h/2) and at theta = 0 only.

## Manifest

Code commit `68bf85672466b33d529cd425d725839f10534de1`, `git_dirty: false`, 2026-10-23T07:32:13 UTC,
`benchmarks/stage18/run.py --workers 12`, runtime 627.3 s (`runtime_seconds` 627.2863504680572). Python
3.12.3, NumPy 2.5.4, SciPy 1.18.1, torch 2.14.1+cpu (1 thread), AMD Ryzen 7 5700U, 16 logical CPUs, Linux
7.0.0-38-generic. `widegap 0.4.0` (installed metadata 0.4.0). `seeds: []`; `dev: false`.

## Not claimed

- Any best-estimate accumulation mobility for NO-nitrided oxide; any field, N_S or temperature dependence.
- That the bracket numbers are mesh- or d_ch-converged (A4, A5 failed).
- That the physical accumulation effect equals the L end; L applies mu_s to all edge electrons in the layer.
- Any ranking of the 13d, stage-17 and 17b designs under L (spread comparable to the mesh sensitivity).
