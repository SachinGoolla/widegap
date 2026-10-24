# Stage 25: generality of the artifact budget under a second impact-ionization coefficient set

Status: **pre-registered** (2026-10-23). No stage-25 code exists. The stage-23
artifact budget used Niwa 2015 coefficients throughout. This stage asks
whether the gain-inflation findings (claims ledger C1a) depend on that choice.
It re-evaluates the **same designs** with **Hatakeyama 2004** coefficients,
which are already in the registry and were used in stages 9b and 12.

## 1. Scope (no re-optimization)

- **AB1′:** the stage-13 initial and final designs, on the sharp device under the stage-12 off-state, mesh h.
  - Claimed: the planar-cut BV gain (Hatakeyama).
  - Real: the dense true-minimum BV change (Hatakeyama, the same 49 + 1 lines as stage 23).
- **AB3′:** the initial, 13, 13c, 13d and 17b final designs on the 0.8 µm geometry.
  - Avalanche BV: the dense true minimum (Hatakeyama, 12c model, h).
  - Blocking voltage: $\min(BV, V_{pt})$, with $V_{pt}$ taken from stage 23 unchanged, since punch-through does not depend on the ionization coefficients.
- **AB6 (new):** the stage-21 initial and final designs (1.4 µm geometry): the dense true-minimum BV gain under Hatakeyama versus Niwa (from stage 21).
- **The Hatakeyama coefficients and orientation convention** are those already used in stage 12.

## 2. Pre-registered predictions

1. **P25a:** the AB1′ artifact fraction is ≥ 0.5. Most of the claimed planar-cut gain is again absent from the dense minimum.
2. **P25b:** every AB3′ design is punch-through-limited ($V_{pt}$ < BV), with an artifact fraction ≥ 0.7.
3. **P25c:** the AB6 dense-BV gain under Hatakeyama is positive and within ±30 % of the Niwa gain (+64.1 V). The stage-21 improvement is not an artifact of the coefficient choice.

Each prediction is reported as confirmed or refuted.

## 3. Criteria

| ID | Criterion | Gating |
|----|-----------|--------|
| H1 | code path: with the Niwa coefficients selected, the recomputed stage-23 AB1 and AB3 avalanche values and the stage-21 dense BVs are reproduced to ≤ 1e-6 relative | yes |
| H2 | every dense line is resolved or explicitly flagged as above 2× its first crossing (stage-23 convention); every solve converged; no "On entry to" BLAS lines | yes |
| H3 | the Hatakeyama dense minimum BV of the initial 0.8 µm design differs by ≤ 3 % between h and h/2 | yes |
| R1 | tables AB1′, AB3′ and AB6 beside their Niwa counterparts, and P25a–c outcomes | reported |

Stage 25 passes if H1–H3 hold.

**Not claimed:** generality beyond these two coefficient sets, or beyond this device family.
