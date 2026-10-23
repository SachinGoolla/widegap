# Stage 16 result: gate NOT passed (H5 failed; H1-H4 and H6 pass) (2026-10-23)

Raw data: `results/stage16/` (`summary.json`, `results.json`, `manifest.json`). Figures:
`figures/stage16_ron_comparison.png`, `figures/stage16_mu_ch_profile.png`. Pre-registration:
`stage16_channel_mobility.md`, including its dated pre-official Note (written after the dev run and
before the official run; no criterion changed). Label: **new** (forward re-evaluation with a
sourced inversion-layer mobility; no optimization was run). All numbers below are quoted from the
JSON files above; anything not recorded there is marked "not recorded".

## Verdict

The stage pass rule is "H1-H6 hold". **H5 failed, so stage 16 does not pass.**

| ID | Criterion | Measured | Threshold | Result |
|---|---|---|---|---|
| H1 | the 11 registry entries pass `check_entry` (plus independent re-verification against rendered page images) | `check_entry: ok` and `usable_in_registry: true` for all 11; `all_usable: true` | all pass | pass |
| H2 | code path with mu_ch = bulk mobility reproduces stage-12 R_on,sp (12b, h) | 2.9189823529917334e-4 vs 2.918982352991735e-4 ohm cm2, rel. 5.6e-16 | <= 1e-12 | pass |
| H3 | long-channel lateral n-MOSFET, mu_eff vs mu_ch at mid-channel | rel. 1.11%, 0.77%, 0.78% for N_A = 1e16, 1e17, 3e17 cm-3 (details below) | <= 2% | pass |
| H4 | Picard convergence, every inner Newton converged | 67 evaluations, `failed: []` | none failed | pass |
| H5 | 12b R_on,sp with the channel model, h vs h/2 | 3.116922e-3 vs 3.215047e-3 ohm cm2, rel. **3.148%** (denominator: the h value) | <= 3% | **fail** |
| H6 | d_ch insensitivity, 20 vs 40 nm (12b, h/2) | rel. 1.54e-8 | <= 1% | pass |

H5 had already failed in the dev run (3.15% vs 3%, disclosed in the spec Note) and the Note
states that it "stands as written". The criterion is not re-registered and no value was tuned.
The only check that this is a mesh effect and not a bug in the Picard loop is indirect: all
Picard iterations converged to < 1.2e-11 relative change (4 iterations each), and H2 shows the
channel code path reduces exactly to the bulk model.

**Consequence.** The channel-model R_on,sp is not mesh-converged to the 3% the stage asked for at
h. It should be quoted at h/2 or h/4 with the mesh change as an uncertainty (below), not at h.

### H5 detail: R_on,sp(mesh), 12b

| Level | Nodes | R_on,sp bulk (ohm cm2) | R_on,sp channel (ohm cm2) | Change from previous level, bulk | Change, channel |
|---|---|---|---|---|---|
| h | 1107 | 2.918982e-4 | 3.116922e-3 | | |
| h/2 | 4293 | 2.852553e-4 | 3.215047e-3 | -2.28% | **+3.15%** |
| h/4 | 16 905 | 2.841923e-4 | 3.255829e-3 | -0.37% | +1.27% |

- With the channel model R_on,sp **increases** with refinement (+3.15% then +1.27%; the change
  shrinks by a factor of about 2.5). The bulk model alone **decreases** (-2.28% then -0.37%). The
  cumulative h to h/4 changes are +4.46% (channel) and -2.64% (bulk).
- Two increments do not establish an order or a limit. A Richardson estimate was not
  pre-registered for this quantity and is not computed here.
- Practical rule: quote the channel R_on,sp at h/2 or h/4; use the last refinement change (+1.27%
  from h/2 to h/4) as the mesh uncertainty. The h/4 value is itself not shown converged.

## R1 (reported): bulk vs channel R_on, V_T, channel share

V_T is the stage-12 transfer-curve threshold, recorded at level 0 only. "Share (drop)" is
definition (i) of the Note: the surface phi_n drop across the gated p-surface span divided by
V_DS. "Share (1-Rb/Rc)" is definition (ii). R_on at V_GS = 15 V, V_DS = 0.1 V (config `V_ON`,
`V_LIN`).

| Case | Level | R_on,sp bulk (ohm cm2) | R_on,sp channel (ohm cm2) | Ratio | V_T bulk (V) | V_T channel (V) | Share (drop) | Share (1-Rb/Rc) |
|---|---|---|---|---|---|---|---|---|
| 12b | h | 2.918982e-4 | 3.116922e-3 | 10.68 | 7.8551 | 9.1527 | 0.9414 | 0.9064 |
| 12b | h/2 | 2.852553e-4 | 3.215047e-3 | 11.27 | not recorded | not recorded | 0.9464 | 0.9113 |
| 12b | h/4 | 2.841923e-4 | 3.255829e-3 | 11.46 | not recorded | not recorded | 0.9478 | 0.9127 |
| 13d initial | h | 2.855541e-4 | 2.289339e-3 | 8.02 | 7.6661 | 8.6544 | 0.9109 | 0.8753 |
| 13d initial | h/4 | 2.780339e-4 | 2.446694e-3 | 8.80 | not recorded | not recorded | 0.9221 | 0.8864 |
| 13d final | h | 2.855541e-4 | 2.286567e-3 | 8.01 | 7.6644 | 8.6492 | 0.9107 | 0.8751 |
| 13d final | h/4 | 2.779204e-4 | 2.398423e-3 | 8.63 | not recorded | not recorded | 0.9196 | 0.8841 |

(13d was not run at h/2.) `share_bulk_drop`, the same drop measure under the bulk model, is
0.374 (12b, h), 0.285 (13d initial, h) and 0.285 (13d final, h).

- **R_on,sp is 8.0 to 11.5 times the bulk-Hall value.** The channel share is about 0.88 to 0.95
  depending on the definition and the case: the channel dominates (about 90%) the on-resistance.
- V_T rises by 1.30 V (12b), 0.99 V (13d initial) and 0.98 V (13d final) with the channel model.
  The registered D_it and Q_f set V_T and N_S in both models; the Noguchi mobility itself does not
  depend on them (caveat below), so the shift is a consequence of the lower channel mobility in
  the V_T extraction; the mechanism was not isolated in this stage.

### Consequence for stages 10-13d

Stages 10-13d optimized the drift and JFET doping using the bulk mobility. Under the channel model
that optimization acted on a minor part of R_on (the non-channel part is about 5 to 13% of the total by the two share definitions). It was not re-optimized with the channel model (an exact adjoint
through a nonlocal mobility was not built).

**13d R_on parity, channel model.** 13d fixed R_on,sp at its initial value in the optimization
(bulk model). Under the channel model:

| Level | Initial (ohm cm2) | Final (ohm cm2) | Final vs initial |
|---|---|---|---|
| h | 2.289339e-3 | 2.286567e-3 | -0.12% |
| h/4 | 2.446694e-3 | 2.398423e-3 | **-1.97%** |

At h the parity holds to 0.12%; at h/4 the final design has a 2.0% lower channel R_on,sp, while
the bulk values differ by -0.04% (2.780339e-4 vs 2.779204e-4). The difference between levels
(0.12% vs 1.97%) is of the order of the mesh change itself, so a parity statement at the 1% level
is not supported by this data; the final design is not worse than the initial design in R_on under
the channel model on either mesh. No h/2 run exists to check the trend.

## R2 (reported): clamps, extrapolation, d_ch, profile

At V_GS = 15 V, level h, d_ch = 20 nm (`counts_15V`):

| Case | Channel columns | Gated columns | Channel edges | N_S clamped (all / gated) | E_eff clamped | N_S > 1e13 | gamma < 0 | N_A out of fit range |
|---|---|---|---|---|---|---|---|---|
| 12b | 15 | 9 | 65 | 6 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| 13d initial | 10 | 6 | 40 | 4 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| 13d final | 10 | 6 | 40 | 4 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |

- The N_S clamp (N_S < 1e11 cm-2) fires only in columns that are not gated (6 of 15 for 12b, 4 of
  10 for 13d), none in the gated span. No column has N_S above 1e13 cm-2; the maximum N_S is
  4.36e12 (12b), 4.97e12 (13d initial) and 4.97e12 cm-2 (13d final). E_eff was never clamped.
- Range over the profile columns (gated spans, x = 0.80 to 1.57 um for 12b and 0.80 to 1.40 um for
  13d), V_GS = 15 V, h:

| Case | mu_ch (cm2/Vs) | E_eff (MV/cm) | N_S (cm-2) |
|---|---|---|---|
| 12b | 17.83 to 24.97 | 0.657 to 0.824 | 2.82e12 to 4.36e12 |
| 13d initial | 17.95 to 24.96 | 0.656 to 0.862 | 2.87e12 to 4.97e12 |
| 13d final | 17.95 to 24.98 | 0.658 to 0.862 | 2.87e12 to 4.97e12 |

  The record's `E_eff_range_MV_cm` field gives 0.574 to 0.824 (12b) and 0.574 to 0.862 (13d); its
  lower end (0.574) is below the profile minimum (0.656-0.658) and is not reconciled in the JSON
  (not recorded which column); the table shows the profile values. The highest mu_ch is at the first column (x = 0.8 um,
  where N_A = 1.5e17 cm-3 and the p-surface doping is lowest); the lowest is mid-channel, where
  N_A = 3e17 cm-3. The N_S values are below the 1e13 cm-2 extrapolation
  threshold of the spec (count 0).
- **d_ch = 10 nm** (12b, h/2): R_on,sp = 3.215028e-3 ohm cm2 vs 3.215047e-3 at 20 nm, relative
  difference 6.0e-6 (and 1.5e-8 for 40 nm vs 20 nm, H6). The R_on is insensitive to the
  thickness within which edges are replaced and N_S is integrated, at this mesh (5, 7, 9 nodes
  within 10, 20, 40 nm at h/2; 3, 4, 5 at h).

## R3 (reported): gamma(N_A) and mu_ph reproduction

Pre-registration check with the registry epsilon_r = 10.2, n_i = 1.0032e-8 cm-3 (recomputed value
in the JSON: 1.003227748e-8), kT/q = 0.025852 V, V_B = 0.

| N_A (cm-3) | N_dep (cm-2) | gamma, log10 (code) | gamma, spec log10 | gamma, ln reading (code) | gamma, spec ln | N20 measured |
|---|---|---|---|---|---|---|
| 1e16 | 5.675e11 | 0.2287 | 0.229 | 0.0916 | 0.092 | 0.20 |
| 3e16 | 9.927e11 | 0.3326 | 0.333 | 0.3309 | 0.331 | 0.34 |
| 1e17 | 1.832e12 | 0.4465 | 0.447 | 0.5931 | 0.593 | 0.46 |
| 4e17 | 3.707e12 | 0.5776 | 0.578 | 0.8948 | 0.895 | 0.54 |

- The code reproduces the spec's log10 values to a maximum absolute difference of 4.95e-4 (rounding
  of the spec values).
- Against the N20 measured gamma, max |difference| is 0.0376 for log10 and 0.355 for ln. The log10
  reading is the one used; the ln reading is rejected by the data, not by assumption.
- N20 Fig. 4b anchor at N_dep = 3.5e12: code gamma (log10) 0.5669 (ln: 0.8702), mu_C0 = 18.20; the
  figure reads about 0.55 and about 19 (read by eye, per the JSON).
- mu_ph (N17 Fig. 6): code gives 163.24 cm2/Vs at E_eff = 0.1 MV/cm, matching the plotted line
  (about 163, per the spec), and 66.5 at 1 MV/cm. Had E_eff been in V/cm the value would be 0.746
  cm2/Vs. Full table (MV/cm to cm2/Vs): 0.02 to 305.8; 0.05 to 213.9; 0.1 to 163.2; 0.2 to 124.6;
  0.5 to 87.1; 1.0 to 66.5.

## H3 detail: long-channel lateral n-MOSFET

L = 10 um, V_DS = 0.05 V, V_GS set by root-finding so that N_S(mid) = 3e12 cm-2 (all three within
1e-4 relative of 3e12), x_mid = 6.0 um.

| N_A (cm-3) | V_GS (V) | E_eff,mid (MV/cm) | mu_ch(mid) | mu_eff = I_D L/(q W N_S V_DS) | Rel. diff. | Local-gradient mu | Rel. diff. (local) |
|---|---|---|---|---|---|---|---|
| 1e16 | 9.244 | 0.2781 | 53.66 | 54.26 | 1.11% | 53.82 | 2.9e-3 |
| 1e17 | 12.198 | 0.5024 | 27.22 | 27.43 | 0.77% | 27.21 | 4.8e-5 |
| 3e17 | 15.382 | 0.7456 | 18.12 | 18.26 | 0.78% | 18.12 | 9.3e-5 |

Reading: mu_eff is defined with the external V_DS and the mid-channel (E_eff, N_S) point, while
the real channel has a non-uniform N_S, E_eff and mu_ch along its length and the n+ regions drop
part of V_DS. The recorded channel voltage drop is 0.04970 V (1e16), 0.04984 V (1e17) and
0.04989 V (3e17) of 0.05 V, so 0.59%, 0.32% and 0.21% of V_DS falls outside the gated channel
(n+ series resistance). That series drop acts to *lower* mu_eff relative to mu_ch, whereas mu_eff
is observed *higher* (by 1.1% to 0.8%); the series resistance therefore explains the order of
magnitude but not the sign of the 0.8-1.1% difference, and the cause of the sign was not
isolated here (candidates: variation of mu_ch along the channel, and the mid-channel point not
being the channel average). The *local-gradient* mobility (current over the local channel field
and N_S) agrees with mu_ch(mid) to within 2.9e-3 for all three, so the implementation is
consistent with the formula to better than 0.3% when the same local quantities are used.
Picard iterations: 1 for each H3 case (the uniform channel starts converged).

## Solver robustness fix (commit d85438d)

The first official-run attempt crashed with a glibc "double free" in a worker (dev-run crash of the
same kind is disclosed in the spec Note; the Note said the root cause was not found at that time).
It has since been diagnosed (by a Sonnet subagent; reviewed):

- **Root cause.** exp overflow in the MOS and quasi-Fermi densities on diverging Newton *trial*
  iterates (potentials up to about 1200 V_T) gave inf/NaN Jacobian entries. The non-finite matrix
  was passed to SuperLU, which made illegal BLAS calls ("On entry to DTRSV ... illegal value"; this
  warning is present in the stage 10, 12 and 13 logs as well) and, once, corrupted memory.
- **Fix.** A `safe_spsolve` wrapper (`src/widegap/linalg.py`) never passes a non-finite matrix to
  SuperLU, and five Newton loops that could have accepted a NaN step now stop and report
  non-convergence (files changed: `adjoint.py`, `adjoint2d.py`, `adjoint_sic.py`,
  `fv/avalanche.py`, `fv/drift_diffusion.py`, `fv/qf_drift_diffusion.py`, `fv2d/equilibrium.py`,
  `fv2d/mos.py`, `fv2d/poisson.py`, `fv2d/qf_dd.py`, `offstate.py`, `sic.py`).
- **Evidence the earlier results are unaffected.** For finite inputs the solve is bit-identical.
  A stage-13d re-run after the fix reproduced all 494 non-timing summary values exactly. **Limit
  of this evidence:** only stage 13d was re-run. Stages 10, 12, 12c/13c, 13, 13b and 15 were not
  re-run; their results are unaffected only in the sense that finite solves cannot change and that
  13d, which shares the code, reproduced exactly. A NaN step accepted in a loop that is now
  guarded could not have been detected from the earlier results alone.

## Interpretation

- The channel resistance is about 90% of R_on,sp under this model (share 0.88 to 0.95 by the two
  definitions). The stage 10-13d optimizations, which optimized the drift and JFET doping with the
  bulk mobility, therefore optimized a minor part of R_on (about 10% by the 1-Rb/Rc measure).
- The 13d R_on parity (initial vs final) holds at h to 0.12% under the channel model; at h/4 the
  final design is 1.97% lower (see above), so the parity is "not worse" rather than "equal" at
  the finest level. This is a statement about the same model-mesh pair, not a measured device.
- The 8.0 to 11.5 times increase in R_on,sp is a model result for Hall factor 1, no N_it
  dependence, and NO-annealed Si-face parameters. **No agreement with measured devices is
  claimed**, in either direction.

## Caveats (from the pre-registration)

- Hall factor taken as 1, as Noguchi does.
- Parameters are for NO-annealed devices only (Si-face, NO-nitrided oxide).
- The Coulomb term does not depend on N_it or Q_f; the registered D_it and Q_f set V_T and N_S
  only.
- 300 K only.
- No mobility model on the n-type accumulation surface (JFET and n+ regions use the bulk model).
- The 13d optimum was found with the bulk model and is not re-optimized here.
- Clamps (N_S 1e11, E_eff 1e-2 MV/cm) are model choices, not paper values; they fired only outside
  the gated span at this operating point.
- Discretization choices recorded in the spec Note: harmonic mean of the two column values per
  edge, N_S integrated over nodes at depth <= d_ch, N_dep with total N_A.

## Manifest

Commit `d85438d0921264e34dc20a9ec6725d902250af3f`, `git_dirty: false`, 2026-10-23T03:35:56 UTC,
`benchmarks/stage16/run.py`, runtime 258.26 s (summary), 4 workers. Python 3.12.3, widegap 0.2.0,
NumPy 2.5.4, SciPy 1.18.1, torch 2.14.1+cpu (1 thread), AMD Ryzen 7 5700U, 16 logical CPUs, Linux
7.0.0-38-generic. Deterministic pipeline, `seeds: []`. `dev_unverified_registry: false`.

## Not claimed

- That stage 16 passed: H5 failed.
- A mesh-converged channel R_on,sp: it moves +3.15% (h to h/2) and +1.27% (h/2 to h/4); the h/4
  value is not shown converged.
- Agreement with measured R_on of any device.
- Hall-effect-corrected, annealed-process, temperature-dependent, N_it-dependent or accumulation-
  layer mobility.
- Any re-optimization with the channel model.
- That stages 10 to 13 and 15 (other than 13d) were re-run after the solver fix.
