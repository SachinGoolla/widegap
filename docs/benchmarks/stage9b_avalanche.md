# Stage 9b: impact ionization and avalanche breakdown (1D)

Status: **pre-registered** (2026-10-21). No stage-9b code exists. The
Chynoweth coefficients here are **study parameters**. 4H-SiC values enter
later through the provenance registry. This replaces the stage-4
field-based proxy with an avalanche criterion. It does not replace the
proxy's verified numerics.

## 1. Models

- Chynoweth ionization coefficients $\alpha_{n,p}(E)=a_{n,p}\exp(-b_{n,p}/|E|)$,
  with $E$ the element field.
- **Avalanche breakdown, ionization integral.** For the reverse-bias
  zero-current field of stage 4:
  $I_n=\int_0^W\alpha_n\exp\!\big(-\int_0^x(\alpha_n-\alpha_p)dx'\big)dx$.
  $BV_{av}$ solves $I_n=1$; the electron- and hole-initiated forms give the
  same root. The integral uses the element-wise exact exponential
  quadrature of piecewise-constant $\alpha$.
- **Generation in DD.** $G=\alpha_n|j_n|+\alpha_p|j_p|$ (scaled) is added
  to both continuity equations, built on the stage-7.1 1D quasi-Fermi
  solver. In reverse bias the current directions are fixed, so $|j|$ is
  smooth on the solution branch. The implementation guards the sign.

## 2. Acceptance criteria

| ID | Criterion | Gating |
|----|-----------|--------|
| B1 | **Exact uniform-field reference.** p⁺–i–n⁺ structure: undoped i-layer of width $W$, heavy contacts, $\delta$ small. In the i-layer the field is uniform to $\le10^{-6}$ relative away from the edge tails, which is verified and reported. The closed form $\alpha_n(1-e^{-(\alpha_n-\alpha_p)W})/(\alpha_n-\alpha_p)=1$ (or $\alpha W=1$ if $\alpha_n=\alpha_p$) gives $E^*$ and hence $BV$. The numerical $BV_{av}$ matches within 1e-3 relative; the residual mismatch is attributable to the edge tails and is reported | yes |
| B2 | The electron- and hole-initiated ionization integrals give the same $BV_{av}$ to ≤ 1e-8 relative (an exact identity of the discrete scheme, which the implementation must ensure) | yes |
| B3 | $BV_{av}$ for a p⁺n diode on meshes $h,h/2,h/4$: observed order in [0.8, 2.5], and the finest mesh within 1% of Richardson | yes |
| B4 | **DD multiplication.** Reverse-biased diode with generation on, small injected current. The multiplication $M=I(\text{with }G)/I(\text{without }G)$ vs the ionization-integral formula $M=1/(1-I_n)$ on the computed field: ≤ 1e-2 relative for $M\le10$. The tolerance allows the self-consistent field change, which is negligible at small current, and the edge definition of $\alpha$ | yes |
| B5 | Analytic Jacobian (with $G$) vs FD, per row, ≤ 1e-6 | yes |
| R1 | $BV_{av}$ vs the stage-4 field-based proxy on the same structures | reported |

Stage 9b passes if B1–B5 hold.
