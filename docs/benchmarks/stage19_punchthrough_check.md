# Stage 19: does the study DMOSFET punch through before avalanche?

Status: **pre-registered** (2026-10-23). No stage-19 code exists. This is
prompted by an exploratory finding from the stage-17c scoping, and checks it
independently.

## 0. The exploratory finding (to be checked, not yet claimed)

**Measure.** The stage-17c scoping measured the electron barrier in the zero-current (stage-12c) off-state of the study DMOSFET (stage-12 geometry and doping, θ = 0, h and h/2). The barrier is $\psi_{surface}-\min_y\psi$ along vertical columns through the p-body under the n⁺ source ($x\in[1.6,2.6]$ µm).

**Result.** The barrier falls below 10 $V_T$ at **≈ 100–120 V** and is about 1–2 $V_T$ from 150 V up. A 1D estimate gives about 150 V: $(qN_At^2/2\varepsilon)(1+N_A/N_D)$, with about 0.3 µm of body at 2.7e17 or more under the source and drift doping 5e16.

**What this would mean.** If confirmed, the source–body–drift n-p-n reaches punch-through far below the avalanche BVs reported for this device in stages 12–17b (about 512–620 V). Those were computed with a zero-current off-state model that cannot carry source-injected electrons, so the model could not see punch-through. **This is a blind spot of the same kind as the stage-13 planar cut.**

## 1. Independent check (different method)

**Method: 1D drift–diffusion through the n⁺/p/n⁻/n⁺ stack.**
- The stack is the doping along the vertical line at $x$ = 2.0 µm and at the worst column found ($x$ = 1.6 µm), taken from the stage-12 device's doping functions.
- The solver is the stage-7.1 quasi-Fermi DD at realistic $n_i/N_0$ with incomplete ionization (stage 9a), no impact ionization, ohmic contacts at both ends, and the p-body **floating** (no body contact in 1D).
  - Floating the body is worst-case for the source–body bias. A second run ties the body to the source by a hole-only contact at the body node layer if feasible. If it is not feasible, that is recorded.
- Sweep $V_D$ from 0 to 400 V (0 to 150 V at minimum, if numerics limit) with continuation.
- Record the electron current density $J_n(V_D)$ and the barrier $\Delta\psi(V_D)$.

**Definition.** $V_{pt}$ is the smallest $V_D$ at which $J_n$ exceeds $10^{3}$ times its value at $V_D$ = 10 V. The barrier at $V_{pt}$ is reported.

## 2. Criteria

| ID | Criterion | Gating |
|----|-----------|--------|
| P1 | solver verification: at $V_D$ = 10 V the 1D DD current agrees between two meshes (h, h/2) to ≤ 1e-3 relative in log, and the $\Delta\psi$ from DD at a low $V_D$ matches the stage-12c zero-current $\Delta\psi$ at the same column to ≤ 1 $V_T$ | yes |
| P2 | every continuation step converged up to the highest $V_D$ reported; no "On entry to" BLAS lines. If the numerics stop below 150 V, the stage reports the reached range and **does not pass** | yes |
| P3 | $V_{pt}$ changes by ≤ 5 % between meshes h and h/2 | yes |
| R1 | $J_n(V_D)$ and $\Delta\psi(V_D)$ curves at both columns; $V_{pt}$; comparison with the scoping estimate (100–150 V) and with the avalanche BV of the same device | reported |
| R2 | **verdict**: confirmed (DD $V_{pt}$ below the avalanche BV) or refuted | reported |

Stage 19 passes if P1–P3 hold. Either verdict is a valid outcome.

**If confirmed**, the earlier result documents (12, 13–13d, 15b, 17, 17b) get a dated note. Their BVs are avalanche BVs of the electrostatic model and are not the device's blocking voltage, because punch-through comes first. A design fix (body depth or a deeper p-well) becomes the next stage, before any stage 17c.

**Not claimed:** a 2D punch-through voltage (1D is a column model); agreement with measured devices.
