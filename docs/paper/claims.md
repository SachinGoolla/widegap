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
| C1 | Adjoint optimizers in SiC MOSFET design exploit unmodelled failure modes: a planar-cut objective was gamed (13), a frozen model classification was exploited for 35 V of a 91 V gain (17 to 17b), and the blocking-voltage claims of a zero-current off-state model were invalidated by punch-through (19) | new (systematic characterization; the phenomenon is known in general) | stage13, 13b, 12c, 17, 17b and 19 results | needs the pre-registered artifact-budget metric (planned) and a pre-registered prediction experiment (experiment A, planned) |
| C2 | Neglecting mobility and ionization doping-dependence in the adjoint gives an X % gradient error and changes the optimum | new (if confirmed) | stage 22 (pre-registered) | planned |
| C3 | Exact adjoint through a nonlocal sheet-charge inversion mobility (Noguchi), verified against FD to ~5e-9 | new (likely; no prior found) | stages 17 and 17b J1 | done |
| C4 | A fair PINN-vs-FV comparison on a drift-diffusion problem found no PINN advantage | new for semiconductor DD (partially) | stage 11 | done; one 1D problem only |
| C5 | Open, verified, provenance-enforced SiC TCAD framework | new (plausibly) | repo; 109 verified parameters | done; check DEVSIM and ChargeTransport.jl before claiming |
| C6 | Channel-dominated R_on (~90 %) means drift/JFET doping optimization addresses a minor part of R_on | model result; textbook-consistent | stages 16 and 17 R2 | done; do not present as a discovery |
| C7 | Punch-through-safe body depth 1.4 µm at +1.4 % R_on for the study cell | model result (design study) | stage 20 | done; study geometry only |

## Gaps a reviewer will raise
- No measured-device validation. Mitigation: compare against sourced published devices; papers to be requested from the owner.
- Generality: one device and one coefficient set. Mitigation: a Hatakeyama arm of an exploit demonstration.
- Narrow gate failures (16 H5, 18b B3, 19 P1/P2) are reported as they are.
