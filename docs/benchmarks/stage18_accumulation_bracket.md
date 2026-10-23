# Stage 18: bracketing the accumulation-layer mobility on the JFET surface

Status: **pre-registered** (2026-10-23). No stage-18 code exists. Stage 16
applies the Noguchi model only to p-type (inversion) surface columns. The
n-type surface under the gate (the JFET accumulation layer) keeps the bulk Hall
mobility. A scoping bracket (stage-17 spec §0, exploratory) found +2.0 % to
+6.2 % in R_on for fixed 20–50 cm²/Vs.

## 1. Sources (full text read; survey on rendered pages)

- **Zhao et al., Solid-State Electron. 51, 1139 (2007)**, doi:10.1016/j.sse.2007.06.006.
  - Device: Si-face 4H-SiC depletion-mode MOSFET, n-channel N_D = 2.4e17 cm⁻³ (SIMS), thermal oxide with **no nitridation**, $C_{ox}$ = 6.9e-4 F/m² (about 50 nm, our inference), room temperature, DC I–V extraction.
  - Accumulation-layer mobility $\mu_s$ = **17.5 ± 0.1 cm²/Vs**, averaged over $V_{GS}$ > 7 V (p. 1141, text). Fig. 4 shows about 20 → 16 cm²/Vs over $V_{GS}$ = 7–11 V.
- Chatty et al., IEEE EDL 22, 212 (2001), and Mater. Sci. Forum 389–393, 1041 (2002) were read and **not used**. Their Hall accumulation values (200–350 cm²/Vs) mix bulk-implant conduction, on 200–250 nm non-nitrided oxides, at N ≲ 1e12 cm⁻².

**Mismatch, recorded now.**
- Zhao's oxide is not nitrided, while the Noguchi inversion model (and the intended device) uses NO-nitrided oxide. Nitridation is known to raise interface-limited mobility. The size of that effect on accumulation mobility is **not sourced here**.
- Zhao's gate overdrive (about 3–7 V above flat band) is below our 15 V.
- Zhao's value is therefore treated as a **pessimistic, sourced lower bound**, not as a best estimate.

## 2. Model (bracket, not a fit)

- **Edges.** On lateral electron edges within $d_{ch}$ = 20 nm of the interface, in n-type surface columns under the gate (the JFET surface), the mobility is set to a constant $\mu_{acc}$.
- **The two bracket ends:**
  - **(L)** $\mu_{acc}$ = 17.5 cm²/Vs (Zhao), the pessimistic end;
  - **(U)** the bulk Hall mobility (stage-16 default), the optimistic end.
- **Velocity saturation** (stage 9a) multiplies the mobility as before.
- **The bracket reported is [R_on(U), R_on(L)].** No best estimate is claimed.

## 3. Acceptance criteria

| ID | Criterion | Gating |
|----|-----------|--------|
| A1 | the registry entry for Zhao's μ_s passes `check_entry`, including independent re-verification on the rendered page | yes |
| A2 | code-path limit: with the (U) setting, stage-16 channel-model R_on is reproduced to ≤ 1e-12 relative (12b, h/2) | yes |
| A3 | every solve (Picard and Newton) converges; no "On entry to" BLAS lines | yes |
| A4 | mesh: R_on(L) changes by ≤ 2 % from h/2 to h/4 (12b) | yes |
| A5 | $d_{ch}$ insensitivity: R_on(L) changes by ≤ 1 % between $d_{ch}$ = 20 and 40 nm (12b, h/2) | yes |
| R1 | bracket [R_on(U), R_on(L)] on h/2 and h/4 for 12b, the 13d initial and final designs, and the stage-17 final design (if stage 17 has run) | reported |
| R2 | the JFET-level sensitivity $dR_{on}/d\theta_5$ (FD) under (L) vs (U) at the 13d initial design: does the accumulation mobility explain the stage-17 §0 factor 2.2–2.4? | reported |

Stage 18 passes if A1–A5 hold.

**Not claimed:**
- any best-estimate accumulation mobility for NO-nitrided oxide;
- any field or $N_S$ dependence;
- temperature dependence.

The caveat "no accumulation-layer mobility model" becomes **"bracketed; NO-nitrided value unsourced"**.
