# Stage 18b: accumulation mobility applied to the gate-induced (excess) charge only

Status: **pre-registered** (2026-10-23). No stage-18b code exists. This is a
disclosed re-registration after stage 18 failed A4 (mesh) and A5 ($d_{ch}$).

## 1. Diagnosis (from stage 18)

- **Stage 18 applied the constant Zhao $\mu_s$ to every electron on the JFET-surface edges within $d_{ch}$.** That includes the flat-band electrons that the donors supply.
- **Zhao's $\mu_s$ belongs to a lumped model, $I_{DS}=I_{ACC}+I_{CH}$.** There $I_{ACC}\propto\mu_s C_{ox}(V_{GS}-V_{FB})$ carries the **gate-induced** charge, and the bulk channel carries the rest at the bulk mobility (Zhao Eq. 1).
- **Consequence for stage 18.** The amount of conductance replaced grew with $d_{ch}$ and with the number of nodes in the layer. Edge counts rose about 4× per refinement. This is consistent with the A4 and A5 failures.

## 2. Model (bracket end L′)

On the same edges as stage 18 (lateral, both nodes within $d_{ch}$, the contiguous n-type surface columns under the gate):

$$\mu_{edge}=\mu_{bulk}\,\frac{n_{fb}}{\bar n}+\mu_{acc}\,\frac{\max(\bar n-n_{fb},0)}{\bar n}$$

- $\bar n$ is the edge-mean electron density.
- $n_{fb}$ is the flat-band (neutral-bulk) electron density at the edge's local doping. It is the ionized donor density from the stage-9a incomplete-ionization model, $N_D^+-N_A^-$, evaluated at the neutral bulk Fermi level.
- If $\bar n\le n_{fb}$ (depletion), $\mu_{edge}=\mu_{bulk}$.
- **Coupling.** $\mu_{edge}$ depends on $n$, so it is lagged in the stage-16 Picard loop with the same tolerance (1e-10) and iteration cap.
- **Velocity saturation** multiplies as before.
- **The other bracket end, U, is unchanged:** bulk mobility everywhere.
- **The reported bracket is [R_on(U), R_on(L′)].**

## 3. Acceptance criteria

| ID | Criterion | Gating |
|----|-----------|--------|
| B1 | code-path limit: with $\mu_{acc}:=\mu_{bulk}$ edge by edge, the stage-16 channel-model R_on of 12b (h/2) is reproduced to ≤ 1e-12 relative | yes |
| B2 | every Picard and Newton solve converges; no "On entry to" BLAS lines | yes |
| B3 | mesh: R_on(L′) of 12b changes by ≤ 2 % from h/2 to h/4 | yes |
| B4 | $d_{ch}$ insensitivity: R_on(L′) of 12b changes by ≤ 1 % between $d_{ch}$ = 20 and 40 nm (h/2) | yes |
| R1 | bracket [R_on(U), R_on(L′)] on h/2 and h/4 for 12b, the 13d initial and final designs, the stage-17 final and the 17b final; comparison with stage-18 L | reported |
| R2 | gate-induced sheet charge on the JFET surface: $\int(\bar n-n_{fb})^+dy$ per column vs $C_{ox}(V_{GS}-V_{FB})/q$, as a consistency check of the "excess" definition | reported |
| R3 | $dR_{on}/d\theta_J$ (FD, 13d initial, h/2) under L′ vs U | reported |

Stage 18b passes if B1–B4 hold.

**Caveats (unchanged):**
- Zhao's oxide is not nitrided.
- Zhao's gate overdrive was lower than ours.
- No field dependence is modelled.
- L′ is a pessimistic sourced bound, not a best estimate.
