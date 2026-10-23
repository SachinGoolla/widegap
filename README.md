# widegap

Verified finite-volume solvers (1D and 2D) for semiconductor device equations, with discrete adjoint sensitivities and constrained doping optimization, built in pre-registered stages toward optimization of wide-bandgap (4H-SiC) power devices. A physics-informed neural network (PINN) is included as a comparison study. **Verified so far (stages 1-13d, see the table below): linear and nonlinear Poisson, drift-diffusion with Scharfetter-Gummel fluxes (1D and 2D, including realistic 4H-SiC n_i/N0), impact ionization and avalanche breakdown, incomplete ionization, field-dependent mobility, a traced 4H-SiC material registry, sourced 4H-SiC device studies (p+n diodes, vertical DMOSFET half-cell), and adjoint-based doping optimization of that DMOSFET with a corner-aware breakdown objective and a gate-oxide-field constraint. All device results are model results from traceable parameters; no agreement with measured devices is claimed. A PINN-vs-FV comparison (stage 11) found no PINN advantage on the problem studied.**

Current release: v0.4.0 (see [CHANGELOG.md](CHANGELOG.md)). This README describes the repository as of 2026-10-23, after v0.4.0 (stage 15 caveat closure, stage 16 channel mobility), with the stage 17 result and the stage 17b and 18 pre-registrations added.

## Status

Stage numbering follows [docs/PROJECT_SPEC.md](docs/PROJECT_SPEC.md) §4. Each stage must pass its own pre-registered gate before the next starts.

| Stage | Content | Status |
|---|---|---|
| 1 | Linear Poisson with heterogeneous permittivity (MMS, FV reference, PINN) | verified (gate passed) |
| 2 | Nonlinear equilibrium Poisson, Boltzmann carriers, 1D abrupt pn junction | verified after one failed run and a disclosed re-registered criterion ([result](docs/benchmarks/stage2_result.md)) |
| 3 | Drift-diffusion with Scharfetter-Gummel fluxes (1D, steady state) | verified after one failed run and a disclosed re-registered tolerance ([result](docs/benchmarks/stage3_result.md)); DEVSIM cross-check reported separately |
| 4 | Field-based breakdown-voltage index in 1D (numerical proxy, no impact ionization) | verified after one failed run and a disclosed re-registration against an exact first-integral reference ([result](docs/benchmarks/stage4_result.md)) |
| 5 | Discrete adjoint doping sensitivities (potential, diode current, breakdown index, field-based BV) | verified on the first official run ([result](docs/benchmarks/stage5_result.md)); adjoint gradient costs 0.98 Newton iterations |
| 6 | Constrained doping optimization (1D: maximize field-based BV at fixed drift-resistance proxy) | verified after one failed run and a disclosed re-registered gradient check ([result](docs/benchmarks/stage6_result.md)); exact field-based BV +11.3% via graded drift doping |
| 7 | 4H-SiC material models (each parameter traced to a primary source that was read) | verified after a provenance block (resolved with user-supplied full texts) and one disclosed geometry re-registration ([result](docs/benchmarks/stage7_result.md)); 51 independently verified parameters |
| 7.1 | Drift-diffusion at realistic 4H-SiC n_i/N0 (quasi-Fermi variables) | verified on the first official run ([result](docs/benchmarks/stage7_1_result.md)); minority currents exact to 4e-14 at n_i/N0 = 1e-24 |
| 8 | 2D finite volume (tensor meshes; linear, equilibrium, quasi-Fermi DD) | verified ([result](docs/benchmarks/stage8_result.md)); exact 1D reductions at n_i/N0 = 1e-24 |
| 8.1 | MOS structures in 2D: oxide, gate, MOSFET half-cell (study-parameter oxide) | verified ([result](docs/benchmarks/stage8_1_result.md)); exact MOS-capacitor reference to 8e-5, I_D mesh order 2.1 |
| 9a | DD physics: non-equilibrium incomplete ionization, field-dependent mobility, interface charge | verified ([result](docs/benchmarks/stage9a_result.md)) |
| 9b | Impact ionization and avalanche breakdown (1D) | verified with a disclosed B1 caveat ([result](docs/benchmarks/stage9b_result.md)) |
| 10 | 2D adjoint sensitivities and doping optimization of a vertical MOSFET half-cell | verified after a disclosed O3 re-registration ([result](docs/benchmarks/stage10_result.md)); field-based BV +36% at fixed I_D (study parameters) |
| 11 | PINN / AD sensitivities vs the FV adjoint (research question) | **new (comparison study)**, gate passed (M1-M5) ([result](docs/benchmarks/stage11_result.md)); **no PINN advantage found**: FV reaches the best PINN error (7.71e-4) in 1.46 ms vs a median 448 s of PINN training; no break-even M for the parametric surrogate; one 1D problem, one untuned config, CPU only |
| 12 | 4H-SiC device studies with sourced parameters (avalanche BV vs doping; DMOSFET on/off) | verified after one disclosed re-registration ([result](docs/benchmarks/stage12_result.md)); BV 1.79 kV at 1e16 (Niwa); the planar-cut DMOSFET BV of 512 V was later found biased low by the stage-12 off-state model (556 V with the stage-12c model, see [12c/13c result](docs/benchmarks/stage12c_13c_result.md)) |
| 13 | Full-physics adjoint optimization of a 4H-SiC DMOSFET | machinery verified; **objective flaw found**: planar-cut BV was gamed while the corner BV fell ([result](docs/benchmarks/stage13_result.md)); 13b fixed the blind spot but exposed an off-state artifact; 12c/13c ([result](docs/benchmarks/stage12c_13c_result.md)) fixed the artifact (true-min BV 517→619 V at fixed R_on) but **failed** on ill-posed peak criteria and a new gate-oxide-field issue; 13d ([result](docs/benchmarks/stage13d_result.md)) adds an oxide-field constraint at 400 V and **passes** D1-D7: true-min BV 517→578 V at fixed R_on and E_ox, costing ~41 V vs 13c; E_ox stays ~5.25 MV/cm (1.3x the 4 MV/cm gate-bias TDDB convention of Shi et al. 2024; off-state drain-induced field, a reference and not a design limit) |
| 14 | GaN | out of scope for now |
| 15 | Caveat closure for the SiC device results ([result](docs/benchmarks/stage15_result.md)); **gate NOT passed (L1 failed)** | **15b** field-line ionization-integral BV passed F1-F5: field-line BV is within -1.09% to +0.06% of the dense vertical-cut BV for the initial, 13b, 13c and 13d designs on h and h/2, and the 13d gain survives (+60.3 V field-line vs +61.1 V cut on h). **15d** local-field mesh adequacy **failed L1** for E_ox (order 0.06, ill-conditioned: the area mean changes ~0.06% per refinement) and W2 (order -0.04, 7.5% from Richardson: genuine non-convergence); W1 and W3 pass. **15f** Ishikawa v_sat (1.7e7 cm/s, isotropic) passed V1-V2: R_on and V_T change by <1%, I_D at V_DS = 15 V by -12% to -20%. 15a restated from 12c C2/C3 (no new run); 15c deferred; 15e blocked on a source; the oxide-field comparison is now against the 4 MV/cm TDDB convention (Shi et al. 2024), a reference only |
| 16 | Channel (inversion-layer) mobility, Noguchi et al. 2017/2020 ([result](docs/benchmarks/stage16_result.md), [spec](docs/benchmarks/stage16_channel_mobility.md)); **gate NOT passed (H5 failed)** | H1-H4 and H6 pass (11 registry entries verified; code path reproduces stage-12 R_on to 5.6e-16; long-channel mu_eff vs mu_ch within 0.77-1.11%; all 67 evaluations converged; d_ch 20 vs 40 nm 1.5e-8). **H5 failed:** the 12b channel-model R_on,sp changes +3.15% from h to h/2 (limit 3%) and +1.27% from h/2 to h/4 (bulk model alone: -2.28%, -0.37%). With the channel model R_on,sp is 8.0 to 11.5 times the bulk-Hall value and the channel carries about 90% of it (0.88 to 0.95 by two definitions); quote at h/2 or h/4 with the mesh change as uncertainty. Stages 10-13d optimized a minor part of R_on; 13d initial vs final channel R_on,sp: -0.12% (h), -1.97% (h/4). Also fixes a solver robustness bug (non-finite matrices passed to SuperLU; commit d85438d). Hall factor 1, NO-annealed parameters, no N_it dependence, 300 K, no accumulation-layer model, 13d not re-optimized; no agreement with measured devices is claimed |
| 17 | Channel-aware re-optimization of the 13d DMOSFET: R_on from the stage-16 channel model, exact adjoint through mu(N_S) ([result](docs/benchmarks/stage17_result.md), [spec](docs/benchmarks/stage17_channel_aware_opt.md)); **gate NOT passed (J4 failed)** | J1-J3 and J5-J7 pass (adjoint vs Richardson FD 5.2e-9; feasible to 1e-12; dense true-min BV 516.9 to 607.6 V on h and 518.1 to 609.7 V on h/2 at fixed R_on,ch and E_ox, against 577.9/580.4 V for the 13d final; h/4 R_on re-check -1.3e-4; 119 non-finite calls, all in rejected Newton trial steps). **J4 failed** on 113 of 173 evaluations: a gated surface column at the JFET/body corner changes classification (p-type to n-type) at theta_J about 0.074 and the optimum (theta_J = 0.147) is past it, so the binding R_on constraint was partly evaluated outside the frozen-classification validity; **the optimum numbers are not trustworthy as design results**. Lagged-gradient error 5.1-5.3% (7-8% for the JFET level); the JFET-level sensitivity is 2.2-2.5 times the bulk value, the excess in the Poisson/ionization (node) part near the JFET surface, not the mobility edges. No accumulation-layer mobility. **17b** (classification guard as an SLSQP constraint, [spec](docs/benchmarks/stage17b_classification_constraint.md)) and **18** (accumulation-layer mobility bracket, [spec](docs/benchmarks/stage18_accumulation_bracket.md)) are pre-registered and not run |

## Installation

```bash
python -m venv .venv
source .venv/bin/activate
pip install torch --index-url https://download.pytorch.org/whl/cpu   # optional, CPU-only
pip install -e ".[test]"
```

Python >= 3.10. The package is named `widegap`. torch is optional: the finite-volume reference solver depends only on NumPy and SciPy, and the PINN requires torch (CI runs the test suite and example both with and without torch).

## Minimal example

```bash
python examples/minimal_poisson.py
```

It solves benchmark case P1 with the finite-volume solver and, if torch is installed, a shortened PINN run. It prints the max error of the FV solution and the relative L2 error of the short PINN run (with its training time). Without torch it prints the FV line and a message that the PINN was skipped. These numbers are not the benchmark results below.

## Benchmark results (stage 1)

The tables below are the stage 1 results only (the generated block is rewritten by `benchmarks/poisson1d/report.py`). Results for later stages are in `docs/benchmarks/stage*_result.md` and the table above.

Problem definitions, error measures and acceptance criteria were fixed before any solver result existed: [docs/benchmarks/poisson_1d.md](docs/benchmarks/poisson_1d.md). One dated amendment (R5 flux balance, restated as R5b/R5c) was also made before any solver output. All quantities are dimensionless.

<!-- results:poisson1d:begin -->
### Finite-volume reference solver (uniform meshes)

Max-norm nodal error against the exact solution (`results/poisson1d/fv.json`):

| Case | N=32 | N=64 | N=128 | N=256 | N=512 | Observed orders |
|---|---|---|---|---|---|---|
| P1 | 8.0e-04 | 2.0e-04 | 5.0e-05 | 1.3e-05 | 3.1e-06 | 2.00, 2.00, 2.00, 2.00 |
| P2 | 6.5e-04 | 1.6e-04 | 4.0e-05 | 1.0e-05 | 2.5e-06 | 2.00, 2.00, 2.00, 2.00 |
| P4 | 1.7e-01 | 3.9e-02 | 9.7e-03 | 2.4e-03 | 6.1e-04 | 2.13, 2.01, 1.99, 2.00 |

- P3 (dielectric interface): max error 5.0e-15 (N=64, interface on a node), 5.3e-14 (N=63, interface inside an element).
- Linear exact solution: max error 5.7e-15.
- R5b discrete flux balance (P1, N=512): relative imbalance 2.1e-13.
- R5c end-element boundary-flux error, N=32..512: 2.5e-03, 6.3e-04, 1.6e-04, 3.9e-05, 9.9e-06; observed orders 2.00, 2.00, 2.00, 2.00.
- FV solve time at N=512: 4.3 ms.

### PINN (3x32 tanh, float64, exact Dirichlet ansatz, Adam then L-BFGS)

Every seed that was run is listed. Relative L2 errors are evaluated on a fixed 2001-point grid.

| Case | Seed | rel. L2 vs exact | rel. L2 vs FV (N=512) | final residual MSE | train time (s) | peak RSS (MiB) |
|---|---|---|---|---|---|---|
| P1 | 0 | 3.2e-07 | 1.6e-06 | 1.6e-08 | 28.2 | 401 |
| P1 | 1 | 2.0e-06 | 3.6e-06 | 9.3e-07 | 26.6 | 402 |
| P1 | 2 | 5.5e-06 | 7.3e-06 | 7.7e-07 | 29.4 | 402 |
| P1 | 3 | 2.5e-06 | 4.1e-06 | 7.5e-07 | 27.7 | 402 |
| P1 | 4 | 4.0e-07 | 2.1e-06 | 4.8e-08 | 29.4 | 402 |
| P2 | 0 | 1.2e-06 | 2.2e-06 | 1.3e-06 | 25.6 | 402 |
| P2 | 1 | 3.1e-07 | 1.7e-06 | 4.2e-08 | 25.8 | 401 |
| P2 | 2 | 1.1e-06 | 2.4e-06 | 5.4e-07 | 26.5 | 403 |
| P2 | 3 | 1.3e-06 | 2.5e-06 | 1.1e-06 | 27.7 | 402 |
| P2 | 4 | 1.0e-06 | 2.1e-06 | 6.5e-07 | 26.2 | 403 |
| P4 | 0 | 3.0e-04 | 3.2e-04 | 1.1e-03 | 26.0 | 401 |
| P4 | 1 | 9.7e-06 | 1.2e-04 | 5.7e-05 | 26.9 | 401 |
| P4 | 2 | 4.2e-04 | 4.3e-04 | 1.6e-03 | 26.8 | 401 |
| P4 | 3 | 2.8e-05 | 1.2e-04 | 1.1e-04 | 27.4 | 401 |
| P4 | 4 | 1.1e-04 | 1.6e-04 | 1.5e-03 | 26.8 | 401 |

**Gate:** passed (R1_P1 pass, R1_P2 pass, R2 pass, R3 pass, R4 pass, R5_linear pass, R5b pass, N1 pass). Pre-registered configuration: True.

Provenance: commit `ce1373d` (clean tree: True), AMD Ryzen 7 5700U with Radeon Graphics, 8 torch threads, CPU only, torch 2.14.1+cpu, NumPy 2.5.4, SciPy 1.18.1, run at 2026-10-21T04:38:06 UTC.
<!-- results:poisson1d:end -->

Reading the tables:

- The "vs FV" column is dominated by the FV solution's own error (about 3e-6 for P1 at N=512), so it is larger than the "vs exact" column. The two comparisons agree to within that level, as pre-registered (N2).
- P4 (steep layer) is a stress test and does not gate. Its PINN error varies by more than an order of magnitude across seeds.
- Cost-to-accuracy: on P1 the PINN needs about 25 s of training to reach the error that the FV solver reaches in milliseconds. The gate shows the PINN code is correct, not that it is competitive.
- Repeated runs of the same configuration on this machine and build gave bit-identical PINN errors. This is not a cross-platform determinism claim.

Figures are generated from the same JSON (`python benchmarks/poisson1d/make_figures.py`):

![FV convergence](figures/poisson1d_fv_convergence.png)
![PINN errors per seed](figures/poisson1d_pinn_errors.png)

Reproduce everything:

```bash
python benchmarks/poisson1d/run.py --out results/poisson1d
python benchmarks/poisson1d/make_figures.py
python benchmarks/poisson1d/report.py   # rewrites the tables above from the JSON
```

## Supported features

Stage 1 core:

- 1D vertex-centred (box-method) finite-volume solver for $-(\varepsilon u')' = f$ with Dirichlet conditions, length-weighted harmonic-mean element permittivity, banded direct solve.
- 1D PINN for the same equation: strong-form residual loss via automatic differentiation, Dirichlet conditions imposed exactly by construction.
- Manufactured-solution problems P1 to P4 and error measures (max error, relative L2, observed order).
- Reproducible benchmark runner that writes per-case and per-seed JSON and a manifest.

Later stages (details in the result documents linked in the status table):

- 1D Newton solvers for equilibrium pn junctions and steady drift-diffusion with Scharfetter-Gummel fluxes, plus quasi-Fermi-variable formulations valid at realistic 4H-SiC n_i/N0 (`widegap.fv`).
- 2D tensor-mesh finite-volume core (linear, equilibrium, quasi-Fermi drift-diffusion, MOS structures) and a vertical DMOSFET half-cell (`widegap.fv2d`, `widegap.devices`).
- Physics models: non-equilibrium incomplete ionization, field-dependent mobility, interface charge, impact ionization and avalanche breakdown.
- A 4H-SiC material registry in which each parameter is traced to a primary source (`widegap.materials`).
- Discrete adjoint sensitivities and constrained SLSQP doping optimization in 1D and 2D, including the full-physics SiC DMOSFET with a corner-aware soft-minimum BV objective and a gate-oxide-field constraint (`widegap.adjoint*`, `widegap.optimize`).
- Corrected zero-current off-state model (stage 12c, `widegap.offstate`).
- PINN-vs-FV comparison study (`widegap.pinn_study`).

## Limitations

Scope and models:

- The solvers are 1D and 2D (tensor meshes only), isothermal (300 K in the SiC device studies), and use a zero-current off-state model for blocking; there is no self-heating, no switching or transient simulation, and no unstructured meshes.
- Device results use model parameters from traced sources (registry) plus study parameters where stated; optimized profiles are model results, not fabrication-ready designs. Agreement with measured devices is not claimed for any stage.
- GaN (stage 14) is out of scope.

Known caveats of the SiC device results (stages 12-13d; status after stages 15, 16 and 17, see [stage 17 result](docs/benchmarks/stage17_result.md), [stage 16 result](docs/benchmarks/stage16_result.md), [stage 15 result](docs/benchmarks/stage15_result.md), [CHANGELOG.md](CHANGELOG.md) and the result documents):

- Mobility is Hall mobility in the bulk model; no interface-scattering (channel) mobility in stages 12-13d, so their R_on is optimistic. **Quantified, gate not passed (stage 16, H5 failed):** with the Noguchi et al. channel mobility (Hall factor 1, NO-annealed parameters, no N_it dependence, 300 K) R_on,sp is 8.0 to 11.5 times the bulk value (12b 10.7x at h, 11.5x at h/4; 13d 8.0x at h, 8.6-8.8x at h/4) and the channel carries about 90% of it; the 12b channel R_on,sp moves +3.15% (h to h/2, above the 3% limit) and +1.27% (h/2 to h/4), so quote it at h/2 or h/4 with the mesh change as uncertainty. A channel-aware re-optimization was run in stage 17 (gate not passed, J4: the frozen surface classification is crossed at the optimum, so its optimum numbers are not trustworthy; 17b is the pre-registered follow-up). The accumulation-layer (n-type surface under the gate) mobility is still the bulk model; a sourced bracket (stage 18, Zhao 2007 lower end 17.5 cm2/Vs) is pre-registered and not run, and a scoping bracket gave +2.0% to +6.2% in R_on. No agreement with measured devices is claimed ([stage 16 result](docs/benchmarks/stage16_result.md), [stage 17 result](docs/benchmarks/stage17_result.md)).
- The off-state uses a zero-current model (phi_n = V_D w with w a harmonic extension across the p-body in the corrected stage 12c model); it is not a full drift-diffusion off-state solution. **Open, bound restated (15a, no new run):** the 12c gates C2 (insensitivity to the w interpolation, relative difference 0, bit-identical) and C3 (carrier self-consistency, n/N at most 2e-19) bound the electron-density error of the interpolation; the hole and current response is not quantified. The original stage-12 model had an artificial electron reservoir on the p side; its corner-column BV values (375-417 V, and 417 to 338 V in stage 13) are artifacts and must not be used quantitatively.
- BV is evaluated with vertical-cut ionization integrals (dense true minimum over at least 40 cuts). **Quantified (15b, gate passed):** a field-line ionization-integral BV on the same off-state fields differs from the dense cut minimum by -1.09% to +0.06% over the initial, 13b, 13c and 13d designs on h and h/2 (13d: -0.08% on h, -0.45% on h/2), and the 13b/13c/13d gains survive (13d: +60.3 V field-line vs +61.1 V cut on h). It is still a post-hoc ionization-integral evaluation of zero-current fields, not a full 2D avalanche simulation (no impact-ionization generation in drift-diffusion).
- The Khan and Cooper v_sat (perpendicular to c, 23 C) was applied to all edges. **Anisotropy closed, value uncertain (15f, gate passed):** the Ishikawa thesis finds v_sat isotropic (1.7e7 cm/s, //c and perpendicular) versus Khan and Cooper 2.2e7 cm/s, a difference the thesis attributes to measurement difficulties and the Hall scattering factor. Switching to 1.7e7 (beta 1.0 and 1.2 evaluated; beta below 1 is outside the model) changes R_on and V_T by less than 1% and I_D at V_DS = 15 V by -12% to -20%. Holes still use the electron v_sat and beta (no registered hole value).
- V_DS is limited to 15 V in the on-state model by double-precision underflow in hole rows. **Open** (15c deferred; needed only for output curves).
- Stage 13d local window fields and gate-oxide field E_ox are not mesh-converged (they rise about 1-6% from h to h/4), so only same-mesh comparisons are used. **Diagnosed, partly open (15d, L1 failed):** the 13d creep is mainly the smooth-maximum bias ln(N)/kappa growing with edge count (e.g. W3: smooth max +6.1% h to h/4, hard max +0.19%). Replacement area-mean measures converge for W1 and W3 (orders 1.09 and 1.42) but not for E_ox (order 0.06, ill-conditioned, change about 0.06% per refinement) or W2 (order -0.04, 7.5% from Richardson). The pre-registered L1 failure stands as written.
- Stage 13d absolute E_ox is 5.25 MV/cm at 400 V, about 1.3 times the 4 MV/cm gate-bias TDDB convention stated by Shi et al. (2024), Mater. Sci. Semicond. Process. 174, 108194, doi:10.1016/j.mssp.2024.108194, p. 6 (no citation given there), before and after optimization. The off-state, drain-induced field has a different stress polarity from positive-gate-bias TDDB, so this is a reference, not a design limit; no oxide reliability claim is made. **Open** (reference comparison only). The 15d area mean of E_ox on h rises by a relative 2.4e-5 from the initial to the final design (it falls 0.09% on h/2); W1-W3 fall on both meshes.
- SLSQP gives local optima from one starting design in a 5-parameter family; the 400 V blocking reference is a choice. Unchanged.
- D_it is uniform from an order-of-magnitude text value; epsilon_ox = 3.9 is an owner-approved convention. Unchanged.
- Stage 11 covers one 1D problem, one untuned PINN configuration, CPU only. Unchanged.

Earlier-stage notes:

- The stage 1-6 solvers are 1D, in scaled (dimensionless) units, and their breakdown and optimization results are proxy results.
- PINN timings are CPU-only and were measured on one machine.
- The PINN error on the steep layer P4 varies by about 40x across seeds.
- Stages 2, 3, 4 and 6 each failed their first run on a pre-registered criterion that was itself flawed (stages 2-3: below the round-off floor; stage 4: depletion-approximation references inadequate at ~100 V_T; stage 6: relative gradient error is ill-posed at a stationary point); each passed after a re-registered criterion disclosed as written after the failure.
- The stage 4-6 field-based breakdown voltage is a proxy: no impact ionization is modelled there, and E_crit is a study constant, not a material value. Avalanche breakdown is modelled from stage 9b on.
- Realistic 4H-SiC n_i/N0 (~1e-24) is verified for equilibrium electrostatics with incomplete ionization (stage 7) and for steady drift-diffusion in quasi-Fermi variables (stage 7.1; complete ionization, constant mobility). The original stage-3 (psi, n, p) solver is valid only for n_i/N0 >= 1e-3.

## Scientific attribution

The method context comes from the following works. **No published device result is reproduced.**

- C. Zhu, P. Andrei, AIP Advances 9, 095301 (2019), [doi:10.1063/1.5113764](https://doi.org/10.1063/1.5113764): adjoint doping-sensitivity method and field-based breakdown index (silicon IGBT).
- P. Brungi, P. Andrei, AIP Advances 16, 035014 (2026), [doi:10.1063/5.0320204](https://doi.org/10.1063/5.0320204): 4H-SiC MOSFET context and comparison targets.
- M. Raissi, P. Perdikaris, G. E. Karniadakis, J. Comput. Phys. 378, 686-707 (2019), [doi:10.1016/j.jcp.2018.10.045](https://doi.org/10.1016/j.jcp.2018.10.045): PINNs.
- P. J. Roache, J. Fluids Eng. 124, 4-10 (2002), [doi:10.1115/1.1436090](https://doi.org/10.1115/1.1436090): method of manufactured solutions.

Classification of this project's work (`docs/literature/requirements_report.md`, section 4):

- *Reproduction*: none claimed. Both published device-level optimization studies cannot be reproduced from the published information (missing supplementary material, model parameters, meshes and impact-ionization coefficients).
- *Independent implementation inspired by published work*: the discrete adjoint doping-sensitivity method and the field-based breakdown index (stages 4 and 5, 1D; stage 10, 2D) and their application to constrained doping optimization of a vertical MOSFET (stage 10) and a 4H-SiC DMOSFET with sourced parameters (stages 12-13d). These are implemented here from the method descriptions and verified against finite differences and exact references; they are not copies of the authors' code and do not reproduce their numbers.
- *New contribution*: the PINN-vs-FV adjoint/sensitivity comparison (stage 11, labelled new, comparison study), whose result is that no PINN advantage was found on the one 1D problem studied; and the corner-aware, oxide-field-constrained objective sequence (stages 13, 13b-13d, 12c), including the finding that the stage 13 planar-cut objective was gamed. These are this project's own studies, not verified against any external result.
- *Unverified*: any statement that PINN or hybrid solvers are competitive for these equations (not supported by stage 11), and any agreement of SiC device results with measurements.

This project is not affiliated with and has not been endorsed by the cited authors.

## Tests

```bash
pytest
```

The test suite covers all stages (`tests/`). The stage 1 checks are listed here as an example of the gate style:

- R1: P1 and P2 observed max-norm orders in [1.9, 2.1].
- R2, R3: max-error thresholds for P1 and P2 at N=512.
- R4: P3 reproduced to round-off with the interface on a node and inside an element.
- R5: a linear solution is reproduced to round-off.
- R5b: discrete flux balance holds to round-off (P1, N=512).
- N1: PINN on P1 reaches rel. L2 error of at most 1e-3 against the exact solution for all five seeds.

The R1 to R5b checks run on the reference solver; the full N1 check is part of the benchmark runner.

## Documentation

- [Project specification v2 (stages, gates, rules)](docs/PROJECT_SPEC.md)
- [Benchmark P1-P4 pre-registration and amendment](docs/benchmarks/poisson_1d.md)
- [Literature and requirements report](docs/literature/requirements_report.md)
- [Stage result documents](docs/benchmarks/) (`stage*_result.md`, with the pre-registrations beside them)
- [Changelog](CHANGELOG.md)
- [Raw results and manifest (stage 1)](results/poisson1d/)

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) and the [Code of Conduct](CODE_OF_CONDUCT.md).

## Citation

See [CITATION.cff](CITATION.cff). Please also cite the original works for the methods used.

## License

Apache License 2.0, see [LICENSE](LICENSE). Cited papers are the property of their authors and publishers and are not redistributed here.
