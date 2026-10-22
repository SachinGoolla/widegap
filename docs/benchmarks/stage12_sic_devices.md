# Stage 12: 4H-SiC device studies with sourced parameters

Status: **pre-registered** (2026-10-21). No stage-12 code exists. This is the
first stage in which the device physics uses **registry values** for 4H-SiC.
The only exception is the owner-approved convention $\varepsilon_{SiO_2}=3.9$.
The adjoint optimization with full physics is **not** in scope; it is a
later stage.

## 1. Parameter mapping (physical → scaled), recorded in the results

| Quantity | Source |
|---|---|
| $T$ | 300 K |
| $N_0$ | $10^{16}$ cm⁻³ |
| $V_T$ | $k_BT/q$ |
| $L_D$ | from $\varepsilon_{\parallel c}$ = 10.20 (Li 2024) for vertical (∥c) devices |
| $n_i(300\,\text{K})$ | registry (stage 7), giving δ |
| Incomplete ionization | stage-7 levels (Götz N h/k, Koizumi Al) |
| Bulk low-field mobility | Ishikawa Caughey–Thomas fits, ∥c for vertical, ⊥c for the lateral channel (Hall mobility; flagged) |
| High-field mobility (⊥c channel) | Khan & Cooper: $v_{sat}$, β, 23 °C |
| Impact ionization (∥c) | Niwa 2015 Okuto–Crowell (primary) and Hatakeyama 2004 Chynoweth (cross-check) |
| Interface | $D_{it}$ uniform near $E_c$ from Dhar 2004 (NO anneal, order of magnitude; flagged); $Q_f$ = 0 with a reported sensitivity at the Afanas'ev upper bound |
| Flat band (n⁺ poly-Si gate) | $V_{FB}=(\Phi_{C,Si}-\Phi_{C,SiC})/q-(E_c-E_F)_{bulk}/q$, using Afanas'ev offsets and the stage-7 statistics. The derivation goes in `docs/equations.md` |
| Gate oxide | $\varepsilon_{ox}$ = 3.9 (owner-approved convention) |

## 2. Studies

- **12a. 1D 4H-SiC p⁺n diodes** (Al p⁺ 10¹⁹, N drift $N_D\in\{10^{15},3\cdot10^{15},10^{16},3\cdot10^{16},10^{17}\}$ cm⁻³, non-punch-through thickness). Compute the avalanche BV by the ionization integral (stage 9b) on the zero-current reverse-bias field (stage 4/7 electrostatics with incomplete ionization), using both coefficient sets.
- **12b. 2D vertical DMOSFET** (the stage-10 geometry rescaled to physical units, with drift doping and thickness set so that BV is in a range the solver handles). Compute the on-state transfer and output characteristics with all sourced physics (stage 9a), and the off-state BV by the ionization integral along the vertical line through the off-state peak-field point.

## 3. Acceptance criteria

| ID | Criterion | Gating |
|----|-----------|--------|
| S1 | Every numerical input traces to a usable registry entry. A machine-generated table (name, value, source DOI or owner approval) is written to the results; nothing is hard-coded | yes |
| S2 | 12a: BV converges under mesh refinement (3 meshes, order in [0.8, 2.5], finest within 1% of Richardson) for each $N_D$ and coefficient set | yes |
| S3 | 12a: electron- and hole-initiated ionization integrals agree to ≤ 1e-8 relative (stage-9b identity) | yes |
| S4 | 12b: every solve converges; source + drain + body current balance ≤ 1e-8 relative where $|I_D|\ge10^{-6}$ (scaled) | yes |
| S5 | 12b: on-current at one on-state bias on meshes $h$ and $h/2$: change ≤ 3% | yes |
| S6 | 12b: analytic Jacobian with all models vs FD, per row, ≤ 1e-6 | yes |
| R1 | 12a: BV vs $N_D$ table and plot for both coefficient sets, and the Niwa/Hatakeyama ratio. Compared qualitatively, with **no fitting**, against the general published scaling; no quantitative claim against any measured device | reported |
| R2 | 12b: $V_T$, $R_{on}$-type metric, $I_D$–$V_{DS}$ with saturation, BV; sensitivity to $D_{it}$ (0, Dhar value, 10×) and to $Q_f$ (0 or the Afanas'ev upper bound) | reported |
| R3 | Every flagged assumption is listed alongside the results: Hall vs drift mobility, $D_{it}$ uniformity, the 23 °C $v_{sat}$, ⊥c/∥c choices, the owner-approved $\varepsilon_{ox}$, and the 1D-cut BV approximation | reported |

Stage 12 passes if S1–S6 hold. A pass means sourced-parameter 4H-SiC
simulations are **numerically verified and fully traceable**. It does not
mean they are **validated against measured devices**, which this stage does
not attempt.
