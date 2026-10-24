# widegap project specification

Version 2, 2026-10-20. This revises the original "master autonomous research
and software engineering prompt" (version 1, drafted with a general-purpose
assistant). Version 1's mission, scientific rules and reporting format are
kept. This version corrects facts that turned out to be wrong, turns goals
into gates that can be checked, and removes structure that produced
paperwork without evidence. §9 lists every change with its reason.

Owner: Sachin Chandra Goolla ([@SachinGoolla](https://github.com/SachinGoolla)).
Licence: Apache-2.0.

## 1. Mission

Build an open, installable Python framework that solves a **declared, verified
subset** of the semiconductor device equations for wide-bandgap power devices,
starting with 4H-SiC. Use it to compute **adjoint doping sensitivities** whose
correctness is demonstrated rather than assumed. PINNs are one backend among
others. They have to earn their place against the conventional solver on
cost-to-accuracy. They are not the project's identity.

The project does not replace commercial TCAD. At every release the README
lists exactly which equations, geometries, models and operating regimes are
supported.

## 2. Scientific anchor (corrected)

| Work | Role | What is actually available |
|---|---|---|
| C. Zhu, P. Andrei, *AIP Adv.* **9**, 095301 (2019), [doi:10.1063/1.5113764](https://doi.org/10.1063/1.5113764) | Origin of the discrete adjoint doping-sensitivity method for power devices, and of the field-based BV index (eq. 1) | Full text (CC BY 4.0). The supplementary material holding models, parameters and the geometry script is **not available to this project** |
| P. Brungi, P. Andrei, *AIP Adv.* **16**, 035014 (2026), [doi:10.1063/5.0320204](https://doi.org/10.1063/5.0320204) | 4H-SiC MOSFET application and comparison targets | Full text (CC BY 4.0). **No supplementary material and no code.** All numerical model parameters are "taken from the Atlas device simulator" and not printed |

Consequences:

- A *reproduction* of either paper's device results is **not feasible** from
  public information. The project must never claim one.
- The achievable and valuable target is an **independent, open
  implementation of the published method**, verified against finite
  differences and analytic limits, applied to devices whose every parameter
  has a cited primary source.
- Secondary sources (LinkedIn, ResearchGate summaries, abstracts) are never
  evidence. Only full texts that have been read count, and the claim-evidence
  matrix records which text each claim came from.

## 3. Gate model

Each stage below follows the same cycle:

1. **Pre-register.** Commit `docs/benchmarks/<stage>.md` with the problem,
   exact or reference solution, error measures, and the **gating** and
   **reported** criteria. This happens before the code that is evaluated. The
   commit timestamp is the evidence.
2. **Implement** the conventional reference first. Then add the PINN or
   other backend, which shares only the problem specification with the
   reference (enforced by a test).
3. **Run** `benchmarks/<stage>/run.py` from a clean tree. It writes raw JSON,
   a manifest (commit, versions, CPU or GPU, threads, seeds, command,
   timestamp) and a `pip freeze`. Every seed is reported.
4. **Amend only forward.** A criterion that turns out to be wrong is changed
   by a dated amendment that explains why, *before* the results it would
   affect. Original text is never rewritten.
5. **Stop on failure.** A failed gate is reported with its cause. Work does
   not move on to the next stage.

## 4. Stages

| # | Stage | Gate, in summary (full criteria pre-registered per stage) | Status |
|---|---|---|---|
| 1 | Linear 1D Poisson, heterogeneous ε (MMS) | FV order 2 ± 0.1; interface exact to round-off; PINN rel. L² ≤ 1e-3, all seeds | **passed** |
| 2 | Nonlinear equilibrium Poisson, Boltzmann carriers, 1D abrupt pn junction (SI units, scaled internally) | Newton converges quadratically; built-in potential equals the analytic $V_{bi}$; agrees with the depletion approximation in its validity regime; MMS order 2 | next |
| 3 | 1D drift-diffusion, Scharfetter–Gummel, SRH | coupled MMS; current constant to round-off; ideal-diode limit; carrier positivity; cross-check against DEVSIM (Apache-2.0) on an identical structure | planned |
| 4 | Field-based breakdown index (ZA19 eq. 1), 1D punch-through diode | matches the analytic field profile; reported as a *proxy*, never as validated avalanche BV | planned |
| 5 | Discrete adjoint sensitivities for stage 2–4 objectives | adjoint equals central finite differences to a pre-registered tolerance at random nodes; measured cost about one linear solve | planned |
| 6 | Constrained doping optimization (bounds, masks, smoothness) | monotone objective; constraints satisfied at every step; final design re-verified by an independent forward solve | planned |
| 7 | 4H-SiC material models | every parameter has a primary source that was read, with units and validity range; anisotropy and incomplete ionization of Al acceptors addressed explicitly | planned |
| 8 | 2D (FV on unstructured meshes), then a SiC MOSFET half-cell | per-stage gates as above | planned |
| 9 | PINN or AD sensitivities compared with the FV adjoint | **research question**, not an assumption: cost-to-accuracy compared against the FV adjoint | hypothesis |
| 10 | GaN (polarization charge, heterointerfaces) | only after stage 8 | out of scope for now |

## 5. Scientific rules

Version 1's rules are kept in full: no fabricated results, parameters or
citations; PINN loss is not evidence of accuracy; no optimization on an
unverified forward solver; no BV claim without a verified breakdown model; no
undocumented clipping or normalization; runtime comparisons record hardware,
tolerances, mesh, budget and stopping criteria; every material parameter has
a traceable source; no merge without passing tests; features are labelled as
implemented, experimental, planned or unsupported. Additions:

- **Cost-to-accuracy is the PINN metric.** Each PINN result is reported next
  to the conventional solver's error and wall time on the same problem. On
  stage 1 the PINN needs about 25 s to reach the error the FV solver reaches
  in about 4 ms.
- **Units.** SI at the API boundary. Internal scaling (thermal voltage,
  intrinsic density, Debye length) is documented in `docs/equations.md` when
  stage 2 introduces it, and never hidden inside a solver.
- **Copyleft hygiene.** GPL and AGPL simulators (ChargeTransport.jl,
  SIMsalabim) may be read for method understanding, but their code is not
  copied.

## 6. Team and model tiering

Version 1's seven agents are reduced to roles that actually produce artifacts:

| Role | Model tier | Owns |
|---|---|---|
| Lead (planning, physics decisions, gate review, integration) | flagship | stage specs, acceptance criteria, final review of every diff |
| Literature | mid tier | bibliography with Crossref-verified DOIs; claim-evidence matrix |
| Implementation | flagship for numerics design, mid tier for routine code | `src/widegap/` |
| Independent verifier | mid tier | tests written from the spec, not from the implementation |
| Mechanical work (fresh-env installs, templates, figure scripts, docs from raw JSON) | small tier | everything that can be checked mechanically |

Rules: one writer per file at a time. Agents report files changed, commands
run, and actual outputs. The lead verifies anything surprising against
primary evidence before accepting it.

## 7. Compute policy

Version 1's staged policy (static checks, small, standard, research) is kept.
Concrete budgets:

- unit tests: under 30 s total on a laptop CPU (CI default);
- stage benchmarks: under 15 min CPU on the development laptop (Ryzen 7
  5700U, 12 GiB, no GPU), with each PINN run in its own subprocess so that
  peak RSS can be attributed to it;
- research sweeps: only after the stage gate passes, with a written budget
  stated before the sweep runs;
- float64 by default; mixed precision only as a tested experiment.

## 8. Repository and release

- The repository grows on demand. A module, document or example is added
  only together with its tests or its evidence (no empty scaffolding).
- Releases are tagged only when every gate up to that stage has passed and CI
  has been observed green. Each release carries its manifest and raw results.
- `docs/repository_metadata.md` holds the GitHub description and topics.

## 9. Changes from version 1, and why

| Version 1 | Version 2 | Reason |
|---|---|---|
| Brungi–Andrei 2026 as the sole foundation, with "supplementary material" to be reviewed | ZA19 is the method origin; BA26 is the application; neither has public parameters or code | Established from the full texts (see `docs/literature/requirements_report.md`) |
| Implicit goal of reproducing published results | Independent implementation only; reproduction declared infeasible | Model parameters are unpublished (Atlas defaults) |
| Package `wbgpinn` | Package `widegap` | The PINN is one backend; naming the project after it prejudges stage 9 |
| "Initial engineering targets" for thresholds | Pre-registered gating and reported criteria, forward-only amendments, stop rule | Turns intentions into evidence that can be audited |
| Seven specialist agents | Five roles with model tiers | Fewer integration conflicts; lower cost |
| LinkedIn and ResearchGate cited as background | Primary sources only | Secondary summaries are not evidence |
| Nine physics and process docs created up front | Docs created by the stage that needs them | Avoids specification drift ahead of code |
| No explicit PINN comparison metric | Cost-to-accuracy against the FV solver | Stage 1 already shows a gap of about 6,000× in wall time |

## 10. Stage-table revision (v2.1, 2026-10-21)

Work after version 2 inserted stages, so the table in §4 is superseded
for stages 7 onwards:

| # | Stage | Status at v0.1.0 |
|---|---|---|
| 7 | 4H-SiC material models (provenance registry) | passed |
| 7.1 | Drift–diffusion at realistic $n_i/N_0$ (quasi-Fermi variables) | passed |
| 8 | 2D tensor-mesh FV (shared edge core) | passed |
| 8.1 | MOS structures in 2D (study-parameter oxide) | passed |
| 9a | DD physics: incomplete ionization, field-dependent mobility, interface charge | passed |
| 9b | Impact ionization and avalanche breakdown (1D) | passed (B1 caveat) |
| 10 | 2D adjoint + vertical MOSFET doping optimization | see result |
| 11 | PINN / AD sensitivities vs the FV adjoint (research question; was 9) | planned |
| 12 | Sourced SiC MOS parameterization (oxide, interface traps, v_sat) and SiC device studies | planned (needs papers) |
| 13 | GaN (was 10) | out of scope for now |

## 11. Stage-table revision (v2.2, 2026-10-22)

Work after v0.1.0 changed several statuses and added stages. This table
supersedes §10 (which is kept unchanged as history). Statuses are taken from
the result documents in `docs/benchmarks/` and from `CHANGELOG.md`.

| # | Stage | Status | Result |
|---|---|---|---|
| 7 | 4H-SiC material models (provenance registry) | passed (after a disclosed geometry re-registration, 7b) | `stage7_result.md` |
| 7.1 | Drift-diffusion at realistic $n_i/N_0$ (quasi-Fermi variables) | passed | `stage7_1_result.md` |
| 8 | 2D tensor-mesh FV (shared edge core) | passed | `stage8_result.md` |
| 8.1 | MOS structures in 2D (study-parameter oxide) | passed | `stage8_1_result.md` |
| 9a | DD physics: incomplete ionization, field-dependent mobility, interface charge | passed | `stage9a_result.md` |
| 9b | Impact ionization and avalanche breakdown (1D) | passed with caveat (B1 passed only after study parameters were scanned) | `stage9b_result.md` |
| 10 | 2D adjoint + vertical MOSFET doping optimization | passed (after a disclosed O3 re-registration) | `stage10_result.md` |
| 11 | PINN / AD sensitivities vs the FV adjoint (research question) | passed (gate M1-M5); no PINN advantage found (one 1D problem, one untuned configuration, CPU only) | `stage11_result.md` |
| 12 | Sourced SiC parameterization and SiC device studies | passed with caveat (after one disclosed re-registration; planar-cut DMOSFET BV later found biased low by the off-state model, see 12c) | `stage12_result.md` |
| 12c | Corrected zero-current off-state model | C2, C3, C5, C6 passed; C1 and C4 failed as defined (C1 ill-posed; C4 showed the stage-12 value was biased low). Model adopted by 13c and 13d | `stage12c_13c_result.md` |
| 13 | Full-physics adjoint optimization of a 4H-SiC DMOSFET | gates passed but the objective was flawed (planar-cut BV gamed while the corner BV fell); superseded | `stage13_result.md` |
| 13b | Corner-aware soft-min BV objective | failed (U1; traced to the stage-12 off-state model), then superseded by 12c/13c | `stage13b_result.md` |
| 13c | 13b re-run on the corrected off-state model | failed (U6: JFET-centre surface field under the gate oxide rose); superseded by 13d. Not re-registered into a pass | `stage12c_13c_result.md` |
| 13d | Oxide-field-constrained, corner-aware optimization | passed (D1-D7), with caveats: absolute $E_{ox}$ about 1.75 times an unsourced 3 MV/cm guideline; local fields not mesh-converged | `stage13d_result.md` |
| 14 | GaN (was 13 in §10, 10 in §4) | out of scope for now | n/a |
| 15 | Caveat closure (planned; sub-items below are independent) | planned, none started | n/a |

Stage 15 sub-items (each needs its own pre-registered spec and gate before
any run; none is pre-registered yet):

| # | Item | Addresses | Dependency |
|---|---|---|---|
| 15a | Off-state with full drift-diffusion (replace the zero-current model) | zero-current off-state limitation | none known |
| 15b | Field-line ionization integral (replace vertical-cut integrals) | vertical-cut BV approximation | none known |
| 15c | High-$V_{DS}$ numerics (remove the 15 V cap) | double-precision underflow in hole rows | none known |
| 15d | Local-field mesh adequacy (converged window fields and $E_{ox}$) | stage 13d non-converged local fields | none known |
| 15e | Channel / interface mobility model | Hall mobility makes $R_{on}$ optimistic | sourced papers required |
| 15f | Anisotropic $v_{sat}$ | $v_{sat}$ (perpendicular to c) applied to all edges | source required |

Numbering note: §10 listed GaN as stage 13; 13 is now the DMOSFET
optimization series and GaN is 14.

## 12. Stage-table revision (v2.3, 2026-10-23)

Work after v2.2 added stages 16 to 18b. These rows extend §11 (kept unchanged as history). Statuses
are taken from the result documents in `docs/benchmarks/` and from `CHANGELOG.md`. Stage 15 (caveat
closure) ran as 15b, 15d and 15f; its gate was not passed (L1 failed), see `stage15_result.md`.

| # | Stage | Status | Result |
|---|---|---|---|
| 16 | Channel (inversion-layer) mobility, Noguchi et al., forward re-evaluation of 12b and 13d | not passed: H5 failed (R_on,sp +3.15% from h to h/2 against a 3% limit); H1-H4 and H6 passed | `stage16_result.md` |
| 17 | Channel-aware re-optimization of the 13d DMOSFET (exact adjoint through the nonlocal mobility) | not passed: J4 failed (113 of 173 evaluations crossed the frozen surface classification); J1-J3 and J5-J7 passed; optimum numbers not trustworthy | `stage17_result.md` |
| 17b | Stage 17 with the classification guard as an SLSQP constraint | passed: J1-J3, J4b, J5-J8 hold; dense true-min BV 572.5 V (h), 35.1 V below the stage-17 final and 5.5 V below the 13d final | `stage17b_result.md` |
| 18 | Accumulation-layer mobility bracket on the JFET surface (Zhao 2007 pessimistic bound) | failed (A4, A5): A1-A3 passed; R_on(L) mesh change 2.71% (limit 2%) and d_ch change 1.63% (limit 1%); bracket L/U - 1 = +3.7% to +9.1% reported with that sensitivity | `stage18_result.md` |
| 18b | Accumulation mobility on the gate-induced (excess) charge only (Zhao 2007 pessimistic bound) | failed narrowly (B3): B1, B2, B4 passed; R_on(L') mesh change 2.098% (limit 2%); d_ch change 0.048% (limit 1%); bracket L'/U - 1 = +2.6% to +7.2% reported | `stage18b_result.md` |
| 19 | Punch-through check of the study DMOSFET (independent 1D drift-diffusion through the n+/p/n-/n+ stack) | not passed: P1 and P2 failed, P3 passed; verdict confirmed (both body configurations): practical onset 86.7 to 127.9 V against the 516.87 V avalanche BV; P1a ill-designed (log-currents at ~1e-45 A/cm2), floating runs stall at 73 to 101 V | `stage19_result.md` |
| 20 | Body-depth redesign against punch-through (d_body 0.8 to 1.6 um) | passed: G1-G3 hold; safe design d_body = 1.4 um (V_pt,es 1.57 / 1.59 x BV on h / h/2; R_on,ch +1.39 %, BV +5.2 % / +4.8 %, V_T unchanged, E_ox(400 V) -26 %); 1.2 um narrowly misses; DD for >= 1.2 um reached only 400 V | `stage20_result.md` |
| 21 | Guarded re-optimization on the punch-through-safe geometry (d_body = 1.4 um) | passed: S1-S6 hold; dense true-min BV 544.00 to 608.10 V (h), 543.15 to 607.76 V (h/2); V_pt,es / BV >= 1.442 at all 81 iterates; R2 comparison with 17b on a different geometry | `stage21_result.md` |
| 22 | Cost of common adjoint simplifications (variants A-D) | passed: Q1-Q3 hold; P22a confirmed (R_on,ch gradient error 27.7 %), P22b refuted (same optimum, 37 vs 10 iterations) | `stage22_result.md` |
| 17c (unblocked) | Channel-aware re-optimization with the p-body doping as a design variable | unblocked by stage 20; to be pre-registered on the 1.4 um geometry | none yet |
