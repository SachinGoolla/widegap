# Claims ledger (working; no prose until the experiments below are recorded)

Each claim gives: its label (reproduction / independent implementation / new / unverified), its evidence (result docs), and its status. Novelty statements rest on a **targeted, not systematic** literature check (2026-10-23, Crossref-verified DOIs). They must be worded accordingly.

## Context facts (with sources)
- **Brungi & Andrei 2026** (doi:10.1063/5.0320204) is the closest prior work: adjoint doping optimization of a 1.2 kV 4H-SiC MOSFET.
  - It uses a Lombardi-type surface mobility with μ_max = 20 cm²/Vs, and BV from the critical-field index max|E|/E_crit = 1.
  - It has no oxide-field constraint, no V_T constraint (V_T fell from 2.7 to 0.7 V in the full-device run) and no punch-through check.
  - It reports no adjoint-vs-FD check.
  - Read via a fetched full-text summary; **locators must be re-checked from the PDF before citation.**
- **Zhu & Andrei 2019** (doi:10.1063/1.5113764): Si IGBT with the same method, read directly. It uses a field-index BV, imposes no constraints on the doping values (Sec. II.D), and reports no FD check.
- **Practice survey:** too thin to claim these simplifications are "common". Say "at least one recent study".

## Candidate claims
| # | Claim | Label | Evidence | Status |
|---|---|---|---|---|
| C1 | Adjoint optimizers in SiC MOSFET design exploit unmodelled failure modes: a planar-cut objective was gamed (13), a frozen model classification was exploited for 35 V of a 91 V gain (17 to 17b), and the blocking-voltage claims of a zero-current off-state model were invalidated by punch-through (19) | new (systematic characterization; the phenomenon is known in general) | stage13, 13b, 12c, 17, 17b, 19 and 23 results | artifact budget done (stage 23, gate passed): AB1 120.9 to 12.4 V (89.7 %), AB2 90.7 to 55.6 V or unchanged (sign model-dependent), AB3 blocking voltage 80.9 to 110.6 V vs avalanche BV 181 to 619 V on the 0.8 um body, AB4 8.7x model error (not an optimizer artifact), AB5 -0.04 V on a different geometry. Observation that avalanche-BV gains (+56 to +102 V for 13c, 13d, 17b) coincided with blocking-voltage changes of -14 to +16 V is **post hoc** (computed from pre-registered definitions, not a pre-registered prediction); the prospective test, experiment A (guard ablation with body doping as a variable), is planned, not run |
| C2 | Neglecting mobility and ionization doping-dependence in the adjoint (variants A-D) gives a 27.7 % R_on,ch gradient error at theta = 0 (BV, E_ox about 1e-6) but, in this problem, does not change the optimum: the D-gradient SLSQP run reaches the same design (dense BV 608.055 vs 608.098 V, both feasible, KKT ratio 0.9 %) at 37 vs 10 iterations. The cost is convergence speed, not the optimum. | new (one problem; constraint-gradient corruption only) | stage 22 (P22a confirmed, P22b refuted) | done; do not claim the simplifications are harmless in general or that any published result is wrong |
| C3 | Exact adjoint through a nonlocal sheet-charge inversion mobility (Noguchi), verified against FD to ~5e-9 | new (likely; no prior found) | stages 17 and 17b J1 | done |
| C4 | A fair PINN-vs-FV comparison on a drift-diffusion problem found no PINN advantage | new for semiconductor DD (partially) | stage 11 | done; one 1D problem only |
| C5 | Open, verified, provenance-enforced SiC TCAD framework | new (plausibly) | repo; 109 verified parameters | done; check DEVSIM and ChargeTransport.jl before claiming |
| C6 | Channel-dominated R_on (~90 %) means drift/JFET doping optimization addresses a minor part of R_on | model result; textbook-consistent | stages 16 and 17 R2 | done; do not present as a discovery |
| C7 | Punch-through-safe body depth 1.4 µm at +1.4 % R_on for the study cell | model result (design study) | stage 20 | done; study geometry only |
| C8 | Guarded optimum on the punch-through-safe geometry (d_body 1.4 um): dense true-min BV 544.00 to 608.10 V (h), V_pt,es / BV >= 1.442 at every iterate (guard 1.2), Pareto at R/R_0 1.05 to 1.2 | model result (design study) | stage 21 | done; electrostatic punch-through criterion only (DD stopped at 400 V for >= 1.2 um); 17b comparison is on a different geometry |

## Gaps a reviewer will raise
- No measured-device validation. Mitigation: compare against sourced published devices; papers to be requested from the owner.
- Generality: one device and one coefficient set. Mitigation: a Hatakeyama arm of an exploit demonstration.
- Narrow gate failures (16 H5, 18b B3, 19 P1/P2) are reported as they are.
