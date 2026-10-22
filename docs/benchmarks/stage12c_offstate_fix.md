# Stage 12c: corrected off-state electrostatics, and stage 13c (13b re-run)

Status: **pre-registered** (2026-10-22). No stage-12c/13c code exists.

## 1. The flaw (stage 13b diagnosis)

The stage-12 zero-current off-state model sets $\phi_n=V_D$ on n-type
nodes and $\phi_n=0$ on p-type nodes. p-type nodes inside the drain-side
depletion region have $\psi\gg0$, so $n=\delta e^{\psi}$ becomes large and
spurious (about 2e18 cm⁻³). The depletion voltage then drops within one
cell, and the corner field never converges.

## 2. Corrected model

- $\phi_p=0$ everywhere, as before. The holes are physical, and they are
  negligible in the n regions under reverse bias.
- $\phi_n=V_D\,w(\mathbf r)$, with $w=0$ on the n⁺ source region (tied to
  the source) and $w=1$ on the drift/JFET/drain n-type region. On the
  remaining p-body nodes $w$ is a **smooth interpolation**: the discrete
  harmonic extension, i.e. Laplace's equation on the mesh graph with those
  Dirichlet values.
- **Justification.** Wherever $0<w<1$ (the p-body), $\psi$ is near the
  p-neutral value or inside a depletion region. Electrons there are
  negligible for any $w$ consistent with the endpoints, except next to
  the drift junction, where $w\to1$ makes them negligible too. This claim
  is checked by C2 and C3, not assumed.
- $w$ depends only on the fixed doping classification at $\theta=0$, as in
  stage 13. The off-state adjoint therefore keeps its form.

## 3. Criteria for stage 12c (corrected model)

| ID | Criterion | Gating |
|----|-----------|--------|
| C1 | the corner peak field at a fixed drain voltage (400 V), initial 13b rounded device, on $h,h/2,h/4$: observed order ≥ 0.8 and the finest mesh within 2% of Richardson | yes |
| C2 | insensitivity to the interpolation. Repeat with a second interpolation ($w$ = normalized graph-distance weighting). The peak field, planar-cut BV and dense-minimum BV each differ by ≤ 1e-6 relative | yes |
| C3 | self-consistency: on every node with $0.01<w<0.99$, $n\le10^{-6}\max(N_A,N_D)$ | yes |
| C4 | the planar-cut BV of the stage-12 12b device (Niwa) agrees with the stage-12 value 512.27 V to ≤ 0.5%. That cut was mesh-stable, so it should be unaffected | yes |
| C5 | a laterally uniform 2D p⁺n (no lateral junction, $N_D=10^{16}$) reproduces the stage-12 12a BV (Niwa) to ≤ 1e-6 relative. The old and new models coincide there | yes |
| C6 | analytic Jacobian of the corrected off-state system vs FD, per row, ≤ 1e-6 | yes |

## 4. Stage 13c

Stage 13c re-runs stage 13b **unchanged**: the same rounded device (r =
0.1 µm), the same K cuts, τ and dense set, and the same U1–U7 thresholds.
Only the off-state model changes, to the stage-12c one. Stage 13c passes if
U1–U7 hold.
