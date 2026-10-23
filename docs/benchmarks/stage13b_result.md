# Stage 13b result: gate FAILED (U1), off-state model artifact identified (2026-10-22)

> **Note (2026-10-23, stage 19):** the BVs reported here are avalanche BVs of the zero-current off-state model, which cannot carry source-injected electrons. Stage 19 ([result](stage19_result.md)) found a 1D punch-through onset of about 90 to 130 V for this geometry (the body under the source is too thin), so these BVs are not the device's blocking voltage. Optimization results here (relative BV changes) remain statements about the avalanche objective only. The original text below is unchanged.

Raw data: `results/stage13b/` (commit `2a76a5c`, clean tree). Implemented by
a Sonnet subagent from the pre-registered spec, then reviewed. The official
run reproduces the development run: **U1 fails**, U2–U7 pass.

## U1 failure and its cause

The corner peak field does not converge. At the initial design's BV it is
4.02 → 4.48 → 4.77 MV/cm on $h$, $h/2$, $h/4$ (observed order 0.67), and
the finest mesh is 9.3% from Richardson against a 2% gate. The cause is
the **stage-12 zero-current off-state boundary model**, not the rounding:

- $\phi_n$ is set to $V_D$ or 0 by local doping type.
- Near the lateral junction this creates an artificial p-side electron
  reservoir (about 2e18 cm⁻³ at 400 V).
- The depletion voltage therefore drops across a single cell next to the
  classification contour.

The artifact is present for every rounding length tried (0.1, 0.15 and
0.2 µm). Caveats were added to the stage-12 and stage-13 results.

## Disclosed process deviations

- **Rounding length.** r = 0.1 µm was chosen after r = 0.2 µm gave a
  path-dependent on-state. A Gaussian source tail left a floating thin
  p layer, the Jacobian was near-singular, and $R_{on}$ varied by 1.7e-4
  between cold starts. The reason is technical, but the choice was made
  after seeing results.
- **Number of development runs.** More than one development run was made:
  one was killed for a 20× speed fix, and the r = 0.2 run was retained in
  scratch.

## What passed (meaningful away from the corner artifact)

| ID | Result | Value |
|---|---|---|
| U2 adjoint vs Richardson FD (initial design) | pass | $BV_{obj}$ 4.3e-9, $R_{on}$ 1.1e-9 |
| U3 final-design gradients | pass | |
| U4 feasibility | pass | |
| U5 dense true-minimum BV (49 independent cuts) | pass | +87.2 V on $h$ (395 → 483 V), +95.5 V on $h/2$; no new blind spot found |
| U6 global peak field at the initial BV | pass | reduced on both meshes |
| U7 SLSQP | pass | exit 0, 15 iterations |

Compared with stage 13: stage 13's optimum, evaluated on this rounded
device, has a *lower* true-minimum BV (379 V vs 395 V). The corner-aware
objective fixed the blind spot. The limiting location stays at the corner
($x\approx0.85$ µm).

**Not claimed:** any corner BV or peak-field value as physical (artifact
above), measured-device agreement, or a global optimum.

## Next (proposed)

Fix the off-state model, for example with stage 12c below. Then re-run
13b unchanged.

- Use full drift–diffusion in the off state, or a quasi-Fermi potential
  that varies with $\psi$ through the depletion region.
- Gate: the corner peak field converges (U1), and the planar-cut BV
  matches the current model where both are valid.
