# Stage 20: p-body depth redesign against punch-through

Status: **pre-registered** (2026-10-23). No stage-20 code exists. Stage 19
confirmed that the study DMOSFET's p-body under the source punches through at
about 87–128 V (1D DD, body tied). That is far below its avalanche BV of
about 517 V. This stage looks for the smallest change to the body geometry
that removes punch-through as the limiting mechanism, and measures its cost.

## 1. Design family (fixed now)

- The stage-12/13d study device at θ = 0, with the rounded body ($r$ = 0.1 µm).
- **Only the body depth changes:** $d_{body}\in\{0.8, 1.0, 1.2, 1.4, 1.6\}$ µm. 0.8 µm is the current device.
- Unchanged: the surface doping, the channel, $N_A$ = 3e17, the source depth (0.3 µm), the drift thickness and doping, the JFET width and doping, and the oxide.
- **Expected effects.** A deeper body makes the JFET region deeper, which raises R_on. It also moves the body/drift junction deeper. V_T should not change.
- **1D estimate (motivation only, not a criterion).** $V_{pt}\approx(qN_At^2/2\varepsilon)(1+N_A/N_D)$ suggests that an effective full-doping body thickness of about 0.6 µm under the source, $d_{body}\approx$ 1.1–1.2 µm, is needed for $V_{pt}$ of about 700 V.

## 2. Quantities per design

1. **Avalanche BV.** The dense true-minimum BV (≥ 40 vertical cuts, Niwa, stage-12c off-state) on h and h/2.
2. **Electrostatic punch-through voltage $V_{pt,es}$.** The smallest $V_D$ at which the minimum electron barrier $\Delta\psi=\psi_{surface}-\min_y\psi$ falls below **10 $V_T$**. The minimum is over all mesh columns crossing the p-body under the n⁺ source and up to the body-contact edge. The barrier comes from the stage-12c off-state.
   - The 10 $V_T$ level is calibrated by stage 19. There, the 10 $V_T$ crossing (95.6–123.5 V) and the DD current onset (86.7–127.9 V) agreed to within about 10 V on the same columns.
3. **DD confirmation (reported).** The stage-19 1D tied-body DD at the worst column, reporting the $J_n$ onsets (1e-3, 1 and 1e3 A/cm²) up to where the numerics reach.
4. **On-state.** R_on,ch on h/2 (stage-16 channel model), the electrostatic V_T (sheet-charge definition, $N_S^*$ = 1e12 cm⁻²), and E_ox at 400 V.

## 3. Criteria

| ID | Criterion | Gating |
|----|-----------|--------|
| G1 | every solve converged for every design; no "On entry to" BLAS lines | yes |
| G2 | mesh adequacy for every design: dense BV changes by ≤ 3 % and $V_{pt,es}$ by ≤ 5 % between h and h/2 | yes |
| G3 | consistency with stage 19: for $d_{body}$ = 0.8 µm, $V_{pt,es}$ (h/2) lies within the stage-19 practical-onset range widened by 10 % (78–141 V) | yes |
| R1 | per design: BV, $V_{pt,es}$, the ratio $V_{pt,es}/BV$, the DD onsets, R_on,ch, V_T, E_ox, and the worst column | reported |
| R2 | **selection (rule fixed now):** the smallest $d_{body}$ with $V_{pt,es}\ge1.2\,BV$ on both meshes is the "punch-through-safe" design. Report its R_on,ch penalty and BV change against 0.8 µm. If none qualifies, report that | reported |

Stage 20 passes if G1–G3 hold. Whether a safe design exists in this family is a finding.

**Not claimed:**
- a 2D current-carrying breakdown simulation (avalanche and punch-through are evaluated separately);
- fabrication feasibility of a deeper body implant;
- agreement with measured devices.

The Pareto question (re-optimizing doping on the safe geometry) is **stage 17c**, which stays on hold until this stage is recorded.
