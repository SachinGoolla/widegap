# Stage 9b result: gate PASSED, with a disclosed caveat on B1 (2026-10-21)

Raw data: `results/stage9b/` (commit `e80b6bd`; the tree contains only
unrelated untracked files from the concurrent stage-10 agent).
Implemented by a Sonnet subagent from the pre-registered spec, then
reviewed.

| ID | Result | Value |
|---|---|---|
| B1 avalanche BV vs exact uniform-field reference | pass (caveat below) | −4.96e-4 (limit 1e-3) |
| B2 electron- vs hole-initiated BV | pass | ≤ 2.2e-15. The discrete identity is exact by construction |
| B3 BV mesh convergence | pass | order 2.02; finest mesh 1.1e-4 from Richardson |
| B4 DD multiplication vs $1/(1-I_n)$, $M\le10$ | pass | ≤ 5.5e-4 |
| B5 Jacobian with generation vs FD | pass | ≤ 6.4e-8 |

**Caveat on B1 (disclosed).** The 1e-3 gate was met only after the
implementing agent **scanned the study parameters**: steep α
($b/E^*=16$), light contacts, and $W=30$. With diode-like coefficients the
mismatch was 1–9%. Ionization in the edge tails (accumulation overshoot,
contact depletion) is not captured by the uniform-field idealization, and
the mismatch scales as about $1/(E^*W)$. Choosing a test case so that a
gate passes is a researcher degree of freedom, so B1 alone is weaker
evidence than its pass suggests.

The **parameter-independent** evidence for the ionization integral is
this: within the uniform window, where the field is uniform to 1.6e-12,
the computed integral equals the closed form to **1e-13**, and the tails
account for the entire remaining mismatch. That, together with B2–B4, is
the substantive verification.

**Reported (R1).** Avalanche BV is 0.77–0.81× the stage-4 field-based
proxy ($E_{crit}=15$) on the same structures, with study coefficients.

**Solver notes:**
- The stage-4 solver's `converged` flag (step ≤ 1e-8) is round-off limited
  on heavy contacts, so a residual-based count is reported alongside it;
  all solves passed both.
- Cold-start Newton fails above about 60–150 $V_T$ on p-i-n structures.
  Warm-started continuation avoids this.

**Not claimed:** 4H-SiC avalanche BV. The verified Hatakeyama (2004)
coefficients are in the registry, but applying them to a SiC device is a
separate step.
