# Stage 15: closing the SiC caveats (15b field-line BV, 15d local-field mesh adequacy, 15f v_sat source)

Status: **pre-registered** (2026-10-22). No stage-15 code exists. This addresses
the caveats listed in the README and in `stage13d_result.md`, in order of how
much each affects a number already reported.

## Triage (recorded before any code)

| Item | Caveat | Reported number it could move | Decision |
|---|---|---|---|
| 15b | BV from vertical-cut ionization integrals | 13d headline (true-min BV 517 → 578 V); the limiting cut x = 1.70 µm lies beside the curved body junction | **do now** |
| 15d | local fields (W1–W3, E_ox) rise 1–6 % per refinement, observed order ≤ 0 | 13d D6 and R1/R2 | **do now** |
| 15f | Khan & Cooper v_sat (⊥c, one doping) applied to every edge | R_on, V_T (expected small at V_DS = 0.1 V) | **do now**: the source is already registered (Ishikawa thesis, Table 3.1) |
| 15a | zero-current off-state (φn = V_D·w), not full DD | off-state fields and BV | **recast**: the 12c C2/C3 gates already bound the electron error of the interpolation. A full-DD off-state at hundreds of volts runs into the same underflow as 15c. No new run; the bound is restated in the result |
| 15c | V_DS ≤ 15 V (δe^{φp−ψ} underflow in hole rows) | none (R_on is at V_DS = 0.1 V) | **deferred**; needed only for output curves. A future scheme must be gated on a 1D reverse-biased diode up to ≥ 1 kV against the existing log-space solution |
| 15e | Hall, not channel, mobility | R_on (optimistic) | **blocked on sources** (paper request sent to the owner) |
| — | 3 MV/cm oxide guideline is unsourced | 13d R2 wording | **blocked on a source** |

## 15b. Field-line ionization integral (post-hoc re-evaluation)

**Definition.**
- At drain voltage $V$, solve the corrected (stage-12c) off-state electrostatics.
- From each seed point, trace the field line $d\mathbf r/ds=\pm\mathbf E/|\mathbf E|$ in the SiC region, in both directions:
  - fixed-step RK4 with step $\Delta s$, using bilinear interpolation of the nodal field reconstructed from $\psi$ on the tensor mesh;
  - stop on leaving the SiC region or where $|\mathbf E|<10^{-3}\max|\mathbf E|$.
- Along each line, evaluate the electron ionization integral, with the same Niwa coefficients and orientation convention as the vertical cuts:
  $I_n=\int_0^L\alpha_n\exp\!\left(-\int_0^{s}(\alpha_n-\alpha_p)\,ds'\right)ds$
  - trapezoidal rule on the RK4 nodes;
  - integrated from the p-side end (the hole-injecting end, as for the cuts).
- The field-line BV is the smallest $V$ at which $\max_{\text{lines}} I_n = 1$, using the same march/secant root-finding tolerances as `dense_bv`.
- Seeds:
  - $N_s$ points uniformly spaced along the arc of the metallurgical junction ($N_D-N_A=0$, linearly interpolated contour) inside the SiC;
  - plus the nodes of the dense vertical cuts at the junction crossing.
  - Default $N_s=400$.

**Scope.**
- Designs: the initial and final designs of 13b, 13c and 13d; re-evaluation only, no optimization.
- Meshes: $h$ and $h/2$.

| ID | Criterion | Gating |
|----|-----------|--------|
| F1 | 1D reduction: on the laterally uniform 2D p⁺n of 12c-C5 ($N_D=10^{16}$, Niwa), the field lines are vertical and the field-line BV matches the stage-12a 1D BV on the same mesh to ≤ 1e-3 relative | yes |
| F2 | step convergence: halving $\Delta s$ changes the field-line BV by ≤ 1e-4 relative (13d initial design, $h$) | yes |
| F3 | seed convergence: doubling $N_s$ changes the field-line BV by ≤ 2e-3 relative (13d initial and final designs, $h$) | yes |
| F4 | mesh adequacy: the field-line BV of the 13d initial design differs by ≤ 3 % between $h$ and $h/2$ | yes |
| F5 | every root-find converged; every line terminated by one of the two stopping rules (no step limit hit) | yes |
| R1 | field-line BV and limiting line (seed position and end points) vs dense vertical-cut BV and limiting cut, for every design on $h$ and $h/2$ | reported |
| R2 | **does the 13d improvement survive?** Field-line BV of the 13d final minus the 13d initial design, on $h$ and $h/2$. Also the same for 13c and 13b | reported |

No sign is pre-registered for field-line BV minus vertical-cut BV; either sign is a valid finding.

## 15d. Local-field mesh adequacy

**Diagnosis (reported).**
- For E_ox, W1, W2 and W3, at 400 V, on $h$, $h/2$ and $h/4$, for the 13d initial and final designs, report:
  - the position of the hard-max edge;
  - its distance to the nearest fixed feature: gate/oxide edge, oxide/SiC interface corner, n⁺ source junction, p-body junction, or window boundary.
- A measure is **pinned** if, on all three meshes, its argmax lies within 2 cells (of that mesh) of the same feature, or on the window boundary.

**Replacement measure (defined now, used for every measure, pinned or not).**
- $\bar E_k$ is the area-weighted mean of $|\mathbf E|$ over a fixed square of side $L=0.1$ µm (the rounding length $r$).
- The square is centred at the hard-max position of the **13d initial design on $h$**; the centre is fixed for all designs and meshes.
- For E_ox the square is clipped to the oxide; for W1–W3 it is clipped to the SiC.
- Edge fields are assigned to dual-face areas.

| ID | Criterion | Gating |
|----|-----------|--------|
| L1 | each $\bar E_k$ (13d initial design) on $h,h/2,h/4$: observed order ≥ 0.8 and the finest within 2 % of Richardson | yes |
| L2 | D6 re-evaluated with $\bar E_k$: does any measure increase from the 13d initial to the final design, on $h$ and $h/2$? | reported |
| R3 | the diagnosis table; for every pinned measure, which feature pins it | reported |

If L1 fails, the stage records that local field maxima in this device are not mesh-resolvable with the present tensor meshes, and states which measure failed.

## 15f. Saturation velocity from the Ishikawa thesis

**Source.** Ishikawa (Kyoto University doctoral thesis, 2025, doi:10.14989/doctor.k26017, already registered), Table 3.1 and §3.4.3.
- Conductance measurements on n-type SiC(11-20) epilayers, //c and ⊥c, at four doping levels.
- Recommended value: $v_{sat}=1.7\times10^7$ cm/s for **both** directions, averaged over $N_D<10^{16}$ cm⁻³ because of self-heating at higher doping.
- Fit exponents α (Eq. 3.3, the Caughey–Thomas form) per doping and direction.
- The thesis attributes the higher 2.2e7 of Khan & Cooper to measurement difficulties on heavily doped samples. It also notes that the Hall scattering factor may make its own velocities overestimates.

**Registry.**
- Add `vsat_n_par_c` and `vsat_n_perp_c_ishikawa` (1.7e7 cm/s), plus per-doping α entries, with locators.
- Each entry is independently re-verified (verified_by ≠ extracted_by), as the registry already requires.
- The Khan & Cooper entries stay; they are not replaced.

| ID | Criterion | Gating |
|----|-----------|--------|
| V1 | the new registry entries pass `check_entry`, including independent re-verification | yes |
| V2 | with the Ishikawa set (isotropic 1.7e7, β = α at the nearest tabulated doping ≤ 1e16), the DD solves converge for every evaluation listed in R4 | yes |
| R4 | R_on,sp, V_T and $I_D$ at $V_{DS}$ = 15 V for the 12b device and the 13d initial and final designs: Khan & Cooper vs Ishikawa | reported |

**Additional caveat recorded now.** The model applies the electron v_sat and β to holes as well (`sic.py`, `vsat_p = vsat_n`). There is no registered hole v_sat. Holes are minority carriers in this n-channel device, and R4 will show whether that matters.

## Stage pass

Stage 15 passes if F1–F5, L1, V1 and V2 hold. The R items are reported whatever the outcome.

**Not claimed:**
- a full 2D avalanche simulation (no impact-ionization generation in DD);
- agreement with measured devices;
- mesh convergence beyond $h/4$.
