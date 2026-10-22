# Stage 9a result: gate PASSED (2026-10-21)

Raw data: `results/stage9a/` (commit `924f47c`). The manifest's
`git_dirty: true` reflects only untracked, unimported work-in-progress files
from concurrently running stage-9b/10 agents; all code under test was
committed. Implemented by a Sonnet subagent from the pre-registered spec,
then reviewed.

| ID | Result | Value |
|---|---|---|
| A1 incomplete-ionization DD at equilibrium = stage-7 solver | pass | 7.1e-15. This is a consistency check: both use the same box scheme |
| A2 MMS with incomplete ionization, 1D and 2D | pass | orders 1.994–2.009 |
| A3 v_sat → ∞ equals constant mobility | pass | 0.0 (bitwise) |
| A4 uniform bar vs closed form $j=qn\mu(F)F$, including saturation | pass | 1.95e-8, the expected δ²/N_D² offset at N_D = 10; saturation reproduced |
| A5 MOS with $Q_f$ and $Q_{it}$ vs exact reference | pass | ≤ 8.4e-5 relative over 24 points; $Q_{it}$ swings from +5 to −3 |
| A6 Jacobian, all models, vs FD | pass | ≤ 3.6e-8 (27 random states, δ down to 1e-24) |
| A7 2D MOSFET with all models: current balance | pass | 4.3e-10 |

**Reported (R1, study parameters).** The threshold shifts −6.9 $V_T$ with
$Q_f$ and −4.5 with $Q_f+Q_{it}$. Field-dependent mobility makes the
output current saturate: $I_D$ = 16.2 at $V_{DS}$ = 16, against 83.5 with
constant mobility.

**Implementation notes, accepted as limitations:**
- A4 needs $N_D\gtrsim10$ for the 1e-6 tolerance, because the hole
  density offset scales as δ²/N_D².
- At $F=0$ the β = 1 mobility derivative uses sign(0) = 0. The Jacobian
  stays continuous there.
- `PhysicsMOSDD` must be driven through `solve_bc`.
- Overflow warnings appear in rejected Newton trial steps during R1. Step
  halving recovers, and every point converged.

**Not claimed:** SiC-specific values. All model constants here are study
parameters; registry-sourced values come in a later step.
