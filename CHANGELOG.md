# Changelog

All notable changes to this project are documented in this file. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added

- **Stage 19**, punch-through check of the study DMOSFET by 1D drift-diffusion through the n+/p/n-/n+ stack at x = 2.0 and 1.6 um, floating and tied body, meshes h and h/2 ([spec](docs/benchmarks/stage19_punchthrough_check.md), [result](docs/benchmarks/stage19_result.md)). **The stage gate is NOT passed: P1 and P2 failed, P3 passed** (commit d568cc0, `git_dirty: false`, widegap 0.5.0). P1a failed (rel_log 0.0015 to 0.0024 against 1e-3) at J_n of about 1e-45 A/cm2, an ill-designed criterion (reported as a finding, not re-registered); P1b passed (|diff| <= 0.42 V_T at 10 V). P2: floating-body runs stall at 73 to 101 V (precision floor in depleted-body hole/electron rows), tied-body runs reach 108 to 225 V. P3: V_pt mesh change 3.04 % and 0.042 %. **R2 verdict: confirmed** for both body configurations. R1: practical onset (J_n = 1e-3 to 1e3 A/cm2, tied runs) 86.7 to 127.9 V, barrier below 10 V_T at 95.6 to 123.5 V, worst column x = 2.0 um, against the 516.87 V avalanche BV. The 512 to 620 V BVs of stages 12 to 17b are avalanche BVs of a zero-current model and not the device's blocking voltage; dated notes were added to the affected result documents. Not claimed: a 2D punch-through voltage, agreement with measured devices.
- **Stage 20**, p-body depth redesign against punch-through, d_body = 0.8 to 1.6 um ([spec](docs/benchmarks/stage20_body_redesign.md), [result](docs/benchmarks/stage20_result.md)). **The stage gate is PASSED** (G1-G3; commit 27fecce, `git_dirty: false`, widegap 0.5.0). G2: BV mesh change <= 0.24 %, V_pt,es mesh change <= 3.00 %; G3: V_pt,es(0.8 um, h/2) 117.49 V inside 78 to 141 V. R2: smallest safe design **d_body = 1.4 um** (V_pt,es 855.4 V (h) / 866.2 V (h/2) against BV 544.0 / 543.1 V, ratios 1.57 / 1.59): R_on,ch +1.39 %, BV +5.2 % (h) / +4.8 % (h/2), V_T unchanged, E_ox at 400 V -26 % (3.91 MV/cm on h/2). 1.2 um narrowly misses (ratios 1.06 / 1.07). DD runs for d_body >= 1.2 um reached only 400 V (numerics), so the claim rests on the electrostatic 10 V_T barrier criterion calibrated in stage 19. Not claimed: a 2D current-carrying breakdown, fabrication feasibility, measured devices; stages 12 to 17b results still refer to the 0.8 um body.

## [0.5.0] - 2026-10-23

### Added

- **Stage 17**, channel-aware re-optimization of the 13d DMOSFET (stage-16 channel-model R_on on h/2 with an exact adjoint through mu(N_S); off-state on h) ([spec](docs/benchmarks/stage17_channel_aware_opt.md), [result](docs/benchmarks/stage17_result.md)). **The stage gate is NOT passed: J4 failed** (113 of 173 evaluations changed the sign of N_D - N_A in a frozen-classification surface column). J1-J3 and J5-J7 passed: adjoint vs Richardson FD 5.19e-9; feasibility R/R0 - 1 = 6.5e-13, E_ox 4.5e-12; dense true-min BV 516.87 to 607.61 V (h) and 518.09 to 609.71 V (h/2), 13d final 577.94/580.41 V; h/4 R_on re-check -1.27e-4; SLSQP exit 0; 119 non-finite calls all in rejected Newton trial steps, 0 "On entry to" lines. R2: lagged-gradient error 5.07% (initial), 5.25% (final); the JFET-level sensitivity is 2.17 to 2.53 times the bulk value, with the excess in the node (Poisson/ionization) part near the JFET surface, not the mobility edges. The optimum numbers are not trustworthy as design results (the binding R_on constraint was partly evaluated outside the frozen-classification validity). Manifest commit `a2425f0`, `git_dirty: false`.
- **Stage 17b**, stage 17 with the surface-classification guard enforced as 35 SLSQP constraints (margin 0.02; J4b, J8, R4) ([spec](docs/benchmarks/stage17b_classification_constraint.md), [result](docs/benchmarks/stage17b_result.md)). **The stage gate is PASSED** (J1-J3, J4b, J5-J8; commit 69eb1c3, `git_dirty: false`). Dense true-min BV 516.87 to 572.46 V (h) and 518.09 to 574.91 V (h/2) at R_on,ch active and E_ox inactive; one active classification constraint (h/4 corner column, x = 0.733 um). The constraint costs 35.1 V (h) against the stage-17 final, whose reclassified R_on,ch is -4.1% (h/2) and -1.8% (h/4). The 17b optimum is 5.5 V below the 13d final (577.94/580.41 V), so no improvement over 13d is claimed. m = 0.02 is not physical; the compensated corner region is unmodelled; accumulation-layer mobility is bracketed in stage 18 (gate failed).
- **Stage 18**, accumulation-layer mobility bracket on the JFET surface (Zhao 2007, 17.5 cm2/Vs, against the bulk model) ([spec](docs/benchmarks/stage18_accumulation_bracket.md), [result](docs/benchmarks/stage18_result.md)). **The stage gate is NOT passed: A4 and A5 failed** (A1-A3 passed; commit 68bf856, `git_dirty: false`). A4: 12b R_on(L) changes 2.71% from h/2 to h/4 (limit 2%); A5: 1.63% from d_ch 20 to 40 nm (limit 1%). R1 bracket L/U - 1 = +3.7% to +9.1% over 12b, 13d initial/final, stage-17 final and 17b final (h/2, h/4), mesh- and d_ch-sensitive. R2: dR_on/dtheta_J is 2.17 times bulk with the U setting (so the stage-17 factor comes from the stage-16 channel model) and a further 2.99 times larger with L. Interpretation (not a gate): the failures stem from applying mu_s to all edge electrons in the layer, including equilibrium ones; an 18b applying it to the excess charge is not pre-registered.
- **Stage 18b**, accumulation mobility applied to the gate-induced (excess) charge only (pre-registered 2026-10-23 after stage 18 failed A4 and A5; [spec](docs/benchmarks/stage18b_accumulation_excess.md), [result](docs/benchmarks/stage18b_result.md)). **The stage gate is NOT passed: B3 failed narrowly** (B1, B2, B4 passed; commit f09529b, `git_dirty: false`). B3: 12b R_on(L') changes 2.098% from h/2 to h/4 (limit 2%; stage 18 A4: 2.71%). B4: 0.048% from d_ch 20 to 40 nm (limit 1%; stage 18 A5: 1.63%), consistent with the stage-18 diagnosis that the d_ch dependence came from the model definition. R_on(U) itself moves +1.27% h/2 to h/4 (derived), the rest of the L' mesh change is unexplained. R1 bracket L'/U - 1 = +2.6% to +7.2% over five designs on h/2 and h/4; R2 excess sheet charge 0.90-0.92 of C_ox(V_GS - V_FB,local)/q; R3 dR_on/dtheta_J L'/U = 2.39.
- Registry entry for the Zhao et al. (2007) accumulation-layer mobility (17.5 cm2/Vs, stage-18 pessimistic bound); the registry now has 109 entries.

### Changed

- README: stage 17 and 17b rows and updated mobility caveat. `docs/PROJECT_SPEC.md`: v2.3 stage-table revision (16, 17 and 18 not passed, 17b passed).

### Fixed

- Run manifests now record the source `__version__` of `widegap` (commit d2f1d9c); the stage 17 manifest still records 0.2.0, the stage 17b manifest records 0.4.0.

### Known limitations

- Stage 17b passes but its optimum is not better than the 13d final in dense true-min BV (5.5 V lower); the 17b margin m = 0.02 is not physical.
- Stage 17 does not pass (J4); its optimum depends on a frozen surface classification the optimum crosses (theta_J about 0.074, optimum 0.147); the reclassified R_on of the final design is not recorded in the official run (4.1% lower in the dev run, per the 17b spec). No accumulation-layer mobility; h/2 on-state mesh uncertainty of about 1.3%.

## [0.4.0] - 2026-10-23

### Added

- **Stage 16**, channel (inversion-layer) mobility, Noguchi et al. 2017/2020 (`widegap` registry entries, sheet model with Picard iteration), forward re-evaluation of 12b and 13d ([spec](docs/benchmarks/stage16_channel_mobility.md), [result](docs/benchmarks/stage16_result.md)). **The stage gate is NOT passed: H5 failed.** H1-H4 and H6 passed (11 registry entries verified; H2 5.6e-16; H3 0.77-1.11% for three N_A; 67 evaluations converged; d_ch 20 vs 40 nm 1.5e-8). H5: 12b channel R_on,sp 3.117e-3 (h) vs 3.215e-3 ohm cm2 (h/2), +3.148% against a 3% limit; h/2 to h/4 +1.27% (3.256e-3); the bulk model alone moves -2.28% and -0.37%. R1: channel R_on,sp is 8.0-11.5 times the bulk-Hall value, the channel carries about 90% of it, and V_T rises by 0.98-1.30 V. 13d initial vs final channel R_on,sp: -0.12% (h), -1.97% (h/4). R2: 0 columns with N_S above 1e13, N_S clamp only outside the gated span, d_ch = 10 nm changes R_on by 6.0e-6. Manifest commit `d85438d`, `git_dirty: false`. No agreement with measured devices is claimed; the 13d design was not re-optimized.
- **Stage 15**, caveat closure for the SiC device results, pre-registered with two disclosed amendments ([spec](docs/benchmarks/stage15_caveat_closure.md), [result](docs/benchmarks/stage15_result.md)). **The stage gate is NOT passed: L1 failed.** 15b and 15f passed their gating criteria; 15d did not.
- **15b**, field-line ionization-integral BV (`widegap.fieldline`), a post-hoc re-evaluation of the 12c off-state fields for the initial, 13b, 13c and 13d designs on h and h/2. F1-F5 passed (1D reduction 3.7e-6 relative, step 4.9e-7, seeds 2.9e-5, mesh 0.22%, 0 step-limit lines). Field-line minus dense vertical-cut BV is between -6.77 V (-1.09%) and +0.29 V (+0.06%); the 13d gain survives (+60.3 V field-line vs +61.1 V cut on h, +61.7 V vs +62.3 V on h/2), as do the 13b and 13c gains. Manifest commit `0d63445`, `git_dirty: false`.
- **15d**, local-field mesh adequacy (`widegap.localfield`). **L1 failed** for E_ox (order 0.056; the area mean changes about 0.06% per refinement, so the pre-registered order is ill-conditioned; the failure stands) and W2 (order -0.044, 7.5% from Richardson, monotone increasing: genuine non-convergence); W1 (order 1.09) and W3 (order 1.42) pass. L2: E_ox area mean rises by a relative 2.4e-5 from the initial to the final design on h; W1-W3 decrease. R3 pinning table and the finding that the 13d smooth-max creep is mainly the ln(N)/kappa smooth-max bias growing with edge count are in the result.
- **15f**, Ishikawa v_sat option (1.7e7 cm/s, isotropic) with 10 new registry entries that pass `check_entry` (V1), and an option in `widegap.sic`. All solves converged (V2). R4: R_on changes by +0.65% to +0.72% (beta 1.0) and +0.08% to +0.09% (beta 1.2), V_T by at most 0.0093 V, I_D at V_DS = 15 V by -20.0% to -20.7% (beta 1.0) and -11.7% to -12.9% (beta 1.2), for the 12b device and the 13d initial and final designs. The Khan and Cooper entries are unchanged.
- 15a restated as a bound from the 12c gates C2 (relative difference 0) and C3 (n/N at most 2e-19) with no new run; 15c deferred; 15e blocked on sources. The 3 MV/cm oxide guideline is replaced by the 4 MV/cm gate-bias TDDB convention of Shi et al. (2024), p. 6, as a reference only: 13d E_ox is 1.3x that, and it is an off-state, drain-induced field with a different stress polarity.

### Fixed

- **Solver robustness** (commit `d85438d`): exp overflow on diverging Newton trial iterates produced inf/NaN Jacobians that were passed to SuperLU (illegal BLAS calls, seen in the stage 10, 12, 13 and 16 logs, and once a glibc double free that crashed the first stage-16 official attempt). A `safe_spsolve` wrapper now never passes non-finite matrices to SuperLU, and five Newton loops reject NaN steps and report non-convergence. Finite solves are bit-identical; a stage-13d re-run reproduced all 494 non-timing summary values exactly. Only 13d was re-run, so earlier stages' results are unaffected on that evidence only.

### Changed

- README caveat list updated with the stage 15 status of each caveat: BV method difference quantified (15b); v_sat anisotropy closed but the value uncertain (15f); local-field convergence diagnosed and partly open (15d); mobility, V_DS limit and the oxide-field reference comparison open. Updated after stage 16: the mobility caveat is quantified (channel R_on,sp 8.0-11.5 times the bulk value) but the stage-16 gate was not passed (H5).

### Process

- Amendment 1 (15f, before the official run): the thesis exponents below 1 are outside the model's domain (beta >= 1 is required), so the Ishikawa set was evaluated at beta = 1.0 and 1.2; a dev run at beta = 1.0 had been seen and the amendment did not choose beta from it.
- Amendment 2 (15b, before the official run): a stagnation stop rule was added after the dev run found lines orbiting a field-reversal sink; F5 is unchanged. The official run has 70 stagnation stops and 0 step-limit lines.
- 15d L1 failed and was not re-registered into a pass.
- The 15d and 15f manifests record `git_dirty: true` only because the (then untracked) 15b files existed at run time; they were committed afterwards in `0d63445` and 15d and 15f do not import them. Commit `2153246`.

### Known limitations

- Stage 16 does not show a mesh-converged channel-model R_on,sp (12b: +3.15% then +1.27% per refinement; H5 failed); Hall factor 1, NO-annealed parameters, no N_it dependence, 300 K, no accumulation-layer model, and 13d not re-optimized under the channel model.
- Stage 15 does not show mesh-converged local fields (E_ox and W2 area means fail L1); no agreement with measured devices is claimed.
- The field-line BV is a post-hoc ionization-integral evaluation of zero-current fields, not a full 2D avalanche simulation.
- The v_sat value remains uncertain between the two published values (2.2e7 and 1.7e7 cm/s); holes use the electron v_sat; beta below 1 is unsupported by the model.

## [0.3.0] - 2026-10-22

### Added

- **Stage 13**, full-physics adjoint optimization of the 4H-SiC DMOSFET (gates T1-T7 passed, but the objective was gamed). Adjoint gradients agree with Richardson FD to 6.7e-9 (R_on) and 2.0e-8 (BV). The optimizer raised the planar-cut BV from 512 V to 633 V while the peak field at the JFET/body corner rose from 4.54 to 6.32 MV/cm and the corner-column BV fell from 417 V to 338 V. This is a specification flaw, not a better device ([result](docs/benchmarks/stage13_result.md)).
- **Stage 13b**, corner-aware soft-min BV objective with a rounded junction (r = 0.1 um). Gate failed on U1 (the corner peak field did not converge: 4.02, 4.48, 4.77 MV/cm on h, h/2, h/4, 9.3% from Richardson against a 2% gate); U2-U7 passed. The failure traced to the stage-12 off-state model, not the rounding ([result](docs/benchmarks/stage13b_result.md)).
- **Stage 12c**, corrected zero-current off-state model (phi_n = V_D w, with w the harmonic extension across the p-body). C2, C3, C5 and C6 passed; C1 and C4 failed (see Fixed and Known limitations) ([result](docs/benchmarks/stage12c_13c_result.md)).
- **Stage 13c**, the 13b problem re-run on the corrected model. True-minimum BV rises 516.9 to 618.8 V at unchanged R_on (518.1 to 622.5 V on h/2). Gate failed on U6 (JFET-centre surface field 2.94 to 3.18 MV/cm, under the gate oxide). The stage-13 optimum evaluated on the corrected rounded device has a true-minimum BV of only 181 V, confirming that stage 13 was gamed ([result](docs/benchmarks/stage12c_13c_result.md)).
- **Stage 13d**, oxide-field-constrained, corner-aware optimization at V_blk = 400 V. All of D1-D7 passed: true-minimum BV 516.9 to 577.9 V on h (518.1 to 580.4 V on h/2) at fixed R_on and fixed E_ox; adjoint vs Richardson FD 1.7e-8 (BV), 1.1e-9 (R_on), 1.3e-12 (E_ox). The constraint costs about 41 V against 13c (618.8 to 577.9 V at R/R_0 = 1), 1.4 V at R/R_0 = 1.2 ([result](docs/benchmarks/stage13d_result.md)).
- **Stage 11**, PINN-vs-FV comparison study, labelled new (comparison study). Gate M1-M5 passed. Q-A: PINN relative L2 error over 5 seeds ranges from 7.71e-4 (best) to 1.71e-2 (worst); FV reaches the best PINN error with N=128 in 1.46e-3 s versus a median PINN training time of 448 s. Q-B: worst-point surrogate dJ/da errors 1.75e-1 to 5.33e-1 across seeds; no break-even number of queries exists, because the PINN inference cost per query (about 2.7 ms) is not below FV at matched accuracy (about 1.1 ms), so training cost is never amortized. Conclusion: no advantage found for this problem ([result](docs/benchmarks/stage11_result.md)).

### Fixed

- **Stage 12 planar-cut BV was biased low.** The artifact in the stage-12 off-state model (phi_n set by local doping type, giving an artificial electron reservoir on the p side, about 2e18 cm-3 at 400 V) also affected the planar junction. Stage 12c gives 556.46 V against the old 512.27 V (+8.6%), within 3.4% of an independent 1D flat-quasi-Fermi bound of 576.16 V (the old value was 11.1% below it). The caveat in `stage12_result.md` was corrected; "mesh-stable" had not meant "correct". The corner-column values (375-417 V in stage 12, 417 to 338 V in stage 13) are artifacts and must not be used quantitatively.
- A corrected off-state model (stage 12c) is added and used by stages 13c and 13d; the original stage-12 model is retained and reproduces 512.27 V. With the corrected model the corner field no longer runs away with refinement (global peak 2.536, 2.545, 2.554 MV/cm; variation 0.7% against 19% before), and a local corner-window measure converges at order 1.07.

### Process

- Stage 13's gates passed while its objective was flawed; this was reported as a specification flaw by the lead, with the lesson that breakdown objectives must cover every location that can break down first.
- Stage 13b failed U1 and stage 12c/13c failed C1, C4 and U6. The C1 criterion (global peak-field convergence order) was ill-posed as written and the 13d criteria replaced it with defined local measures; no failed gate was re-registered into a pass. 13c recorded `git_dirty: true` only because of the untracked 12c results folder.
- Stage 13b disclosed process deviations: the rounding length r = 0.1 um was chosen after r = 0.2 um gave a path-dependent on-state (R_on varied by 1.7e-4 between cold starts), and more than one development run was made.
- Stage 13 disclosed deviation: the T3 tolerance was applied to R/R_0 because the literal tolerance is unit-dependent and vacuous for R_on in ohm cm2.
- Stage 11 was amended before its official run, with disclosure: Q-B as written was degenerate (symmetric doping gives J = 0 and dJ/da = 0 identically), so only the n side is scaled and the evaluation points avoid a = 1. The PINN configuration was unchanged. Official run: `--jobs 1`, total runtime 6968 s, commit f19ccce, `git_dirty: false`.
- The stage 13d manifest records `git_dirty: true` because untracked stage-11 files (not imported by 13d) existed at run time; they were committed afterwards.

### Known limitations

- Stage 13d: absolute E_ox is 5.25 MV/cm at 400 V, about 1.75 times an unsourced 3 MV/cm guideline, both before and after optimization, and rises to 6.20 MV/cm at the new rated BV (5.90 MV/cm for the initial design at its own BV). The constraint only prevents worsening at 400 V; no oxide reliability claim is made.
- Stage 13d local window fields and E_ox are not mesh-converged (they rise about 1-6% from h to h/4; all eight convergence flags fail), so absolute values are mesh-dependent and only same-mesh comparisons are used.
- Stage 12c C1 and C4 failed as defined, and the 13c U6 failure motivated 13d; the stage 13 and 13b optima are not better devices.
- BV uses vertical-cut ionization-integral approximations of the 2D device (dense true minimum over at least 40 cuts), not a full 2D avalanche solution; SLSQP gives local optima from one starting design in a 5-parameter family. The V_blk = 400 V reference is a choice.
- Stage 11 covers one 1D problem, one fixed untuned PINN configuration, CPU only (16 torch threads, one job at a time). It does not show that PINNs are generally uncompetitive. The parent-process manifest lists torch's default of 8 threads; each PINN subprocess set 16 and recorded it.
- No agreement with measured devices is claimed for any stage 12-13d result.

## [0.2.0] - 2026-10-21

### Added

- **Stage 12**, first 4H-SiC device studies with sourced parameters (pre-registered; gate passed after one disclosed re-registration). Avalanche breakdown of 4H-SiC p+n diodes (non-punch-through, parallel to c, 300 K) at N_D = 1e15 to 1e17 cm⁻³, with Niwa 2015 and Hatakeyama 2004 coefficients giving 11.2 kV to 322 V and 14.2 kV to 310 V; the two coefficient sets differ by up to 21% in breakdown voltage. Vertical 4H-SiC DMOSFET half-cell: V_T = 7.9 V, R_on,sp ≈ 2.9e-4 Ω·cm², planar-cut avalanche BV = 512 V (Niwa; 522 V Hatakeyama), mesh-stable to 0.03%. Oxide-interface sensitivities: V_T = 7.6 V to 10.6 V over D_it ×0 to ×10, and 4.3 V to 11.3 V for Q_f = ±1.5e12 cm⁻² ([result](docs/benchmarks/stage12_result.md)).
- Material registry grown to 87 entries, including impact ionization (Hatakeyama 2004, Niwa 2015), high-field transport (Khan & Cooper 2000), interface-trap densities (Saks 2000, Dhar 2004) and SiC/SiO₂ band offsets (Afanas'ev 1996) ([stage 7 amendments](docs/benchmarks/stage7_sic_materials.md)).
- Owner-approved-convention parameter class (`approved_by`, `approved_on`, `rationale` required). Thermal SiO₂ static relative permittivity ε_r = 3.9 is the single documented use of this class; no other parameter uses it.

### Fixed
- A commit briefly pushed a failing test because a piped `pytest | tail` masked pytest's exit code; fixed in d954698, and commits are now gated on pytest's own exit status.

- Material registry raised a KeyError for owner-approved entries that have no DOI or locator; the registry now tolerates them.

### Process

- Stage 12 run 1 failed S2′ (mesh-convergence order 2.90 at N_D = 1e15 Niwa, superconvergence). The criterion was re-registered as S2″ with no order upper bound, and run 2 passed S1–S6. The amendment is disclosed in the stage 12 spec.
- Stage 13 (full-physics adjoint optimization of a 4H-SiC DMOSFET) pre-registered, not yet run.

### Known limitations

- Stage 12 uses Hall, not drift, mobility, with no interface-scattering (channel) mobility model, so R_on is optimistic.
- The Khan & Cooper v_sat (perpendicular to c, 23 °C) is applied to all edges.
- D_it is uniform, taken from an order-of-magnitude text value.
- ε_ox is an owner-approved convention, not a read primary source.
- Avalanche BV uses a 1D cut, not the full 2D ionization integral; the corner cut (375–417 V) is mesh-dependent and reported only.
- V_DS is limited to 15 V by double-precision underflow in hole rows.
- At N_D = 1e17 cm⁻³ the Niwa peak field (3.4 MV/cm) exceeds its stated 1–3 MV/cm validity.
- No agreement with measured devices is claimed; outputs are model results from traceable parameters.

## [0.1.0] - 2026-10-21

### Stages

- **Stage 1**, linear 1D Poisson with heterogeneous permittivity: finite-volume reference solver, max-norm order 2.00 on uniform meshes, interface case exact to round-off ([result](docs/benchmarks/poisson_1d.md)).
- **Stage 2**, nonlinear equilibrium Poisson for a 1D abrupt pn junction: Newton solver, MMS order 2.00, bulk neutrality and peak field checked against the depletion approximation; gate passed after re-registration ([result](docs/benchmarks/stage2_result.md)).
- **Stage 2b**, re-registered quadratic-Newton criterion S2′ used for the passing stage 2 rerun ([result](docs/benchmarks/stage2b_newton_criterion.md)).
- **Stage 3**, 1D steady-state drift-diffusion with Scharfetter-Gummel fluxes and optional SRH; method of manufactured solutions and Jacobian checks; current agrees with an independent DEVSIM cross-check to 3.7e-8 at V = 1 V_T ([result](docs/benchmarks/stage3_result.md)).
- **Stage 3b**, re-registered Newton stopping tolerance for the stage 3 gate; passing run uses the unchanged physics ([result](docs/benchmarks/stage3_result.md)).
- **Stage 4**, field-based breakdown index in 1D as a numerical quantity; BV extraction |B−1| ≤ 2.4e-14 and junction field within 6e-5 of the exact first integral ([result](docs/benchmarks/stage4_result.md)).
- **Stage 4b**, re-registered criteria against the exact first-integral reference; passing run ([result](docs/benchmarks/stage4_result.md)).
- **Stage 5**, discrete adjoint doping sensitivities for the stage 2–4 objectives; worst adjoint-vs-finite-difference errors 4.7e-7 (diode current) and 6.8e-7 (field-based BV); one adjoint gradient costs 0.98 Newton iterations ([result](docs/benchmarks/stage5_result.md)).
- **Stage 6**, constrained 1D doping optimization maximizing field-based BV at fixed drift-resistance proxy; exact field-based BV rises 80.66 → 89.80 V_T (+11.3%) on h = 0.025 ([result](docs/benchmarks/stage6_result.md)).
- **Stage 6b**, re-registered gradient-check criterion for the stationary optimum; passing run ([result](docs/benchmarks/stage6_result.md)).
- **Stage 7**, 4H-SiC material models with 51 independently verified parameters, each traced to a read primary full text; incomplete-ionization equilibrium at realistic n_i/N0 ([result](docs/benchmarks/stage7_result.md)).
- **Stage 7b**, re-registered geometry for the stage 7 M5 check after a failed run; passing run, M5(b) at 2.1e-14 V_T ([result](docs/benchmarks/stage7_result.md)).
- **Stage 7.1**, drift-diffusion at realistic 4H-SiC n_i/N0 using quasi-Fermi variables; minority currents exact to 3.9e-14 relative, Jacobian vs finite differences ≤ 4.7e-7 ([result](docs/benchmarks/stage7_1_result.md)).
- **Stage 8**, 2D finite-volume box method on tensor meshes (linear, equilibrium, quasi-Fermi DD); MMS orders ≈ 2, exact 1D reduction at n_i/N0 = 1e-24 ([result](docs/benchmarks/stage8_result.md)).
- **Stage 8.1**, 2D MOS structures and a MOSFET half-cell with study-parameter oxide; MOS-capacitor reference agreement ≤ 8.4e-5 relative, drain-current mesh order 2.12 ([result](docs/benchmarks/stage8_1_result.md)).
- **Stage 9a**, drift-diffusion physics: non-equilibrium incomplete ionization, field-dependent mobility, and interface charge; Jacobians vs finite differences ≤ 3.6e-8 ([result](docs/benchmarks/stage9a_result.md)).
- **Stage 9b**, impact ionization and avalanche breakdown in 1D; discrete multiplication matches 1/(1−I_n) to 5.5e-4 and BV mesh order 2.02, with a disclosed B1 caveat ([result](docs/benchmarks/stage9b_result.md)).
- **Stage 10**, 2D adjoint sensitivities and doping optimization of a vertical MOSFET half-cell: 2D adjoint + vertical DMOSFET doping optimization: field-based BV +36% at fixed on-current (study parameters); adjoint vs FD 6e-8 / 3.8e-7; Pareto set and 2D sensitivity maps ([result](docs/benchmarks/stage10_result.md)).

### Known limitations

- All device constants in stages 8.1, 9a, 9b and 10 are study parameters (for example oxide permittivity ratio 0.38 and flat-band voltage in stage 8.1); no SiC-specific avalanche or MOSFET prediction is claimed ([stage 8.1](docs/benchmarks/stage8_1_result.md), [stage 9b](docs/benchmarks/stage9b_result.md)).
- Stage 7 material model: N_v uses parabolic band-edge masses, so N_v and n_i are underestimates; the g_D = 2 degeneracy ignores the 1s(E) valley-orbit state; mobility fits are Hall mobilities, not drift mobilities; the hole μ_min⊥c conflict (40 vs 20) is unresolved; E_g(T) at 500 K is extrapolated ([result](docs/benchmarks/stage7_result.md)).
- Stage 7 model outputs (E_g, n_i, built-in potential, ionized fraction) are not validated against device measurements.
- Stage 8.1 has no interface traps, no fixed charge, constant mobility, and complete ionization ([result](docs/benchmarks/stage8_1_result.md)).
- Stage 4 field-based BV is a numerical proxy with no impact ionization, and the max-element field definition converges only at first order ([result](docs/benchmarks/stage4_result.md)).
- Stage 9b avalanche B1 check passed only after the study parameters were scanned; with diode-like coefficients the mismatch was 1–9%. Avalanche BV is 0.77–0.81 times the stage 4 proxy with study coefficients. Cold-start Newton fails above about 60–150 V_T on p-i-n structures ([result](docs/benchmarks/stage9b_result.md)).
- Stage 7.1 resolves majority-carrier current only for |j| ≳ 1e-6 in scaled units (about 1e-2 A/cm²) ([result](docs/benchmarks/stage7_1_result.md)).
- Stage 10 uses the field-based BV as a proxy: the breakdown index is smoothed, excludes a box around the source/body junction because of study scaling, and no quantitative comparison with the published devices is claimed ([stage 10 pre-registration](docs/benchmarks/stage10_mosfet_optimization.md)).
- Stages 8 and 8.1 use tensor meshes only; unstructured meshes are not implemented ([result](docs/benchmarks/stage8_result.md)).
- Benchmarks were run on CPU only (stated in the stage 8 result).

### Process

- Every stage was pre-registered before its code was evaluated; each gate was run from a documented command and a recorded manifest ([spec §3](docs/PROJECT_SPEC.md)).
- Several criteria were re-registered with disclosure after failed gates: stage 2 (S2′, [stage2b](docs/benchmarks/stage2b_newton_criterion.md)), stage 3 (tolerance, [stage3b](docs/benchmarks/stage3b_tolerance.md)), stage 4 (first-integral criteria, [stage4b](docs/benchmarks/stage4b_first_integral.md)), stage 6 (gradient check, [stage6b](docs/benchmarks/stage6b_gradient_check.md)), and stage 7 (geometry, [stage7b](docs/benchmarks/stage7b_geometry.md)).
- Some criteria were amended before their official runs, after development runs, and disclosed in the stage documents: stage 8 (P9), stage 8.1 (G2, G6), and stage 10 (O3′ and two deviations).
