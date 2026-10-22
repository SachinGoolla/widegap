# Stage 13d: oxide-field-constrained, corner-aware SiC DMOSFET optimization

Status: **pre-registered** (2026-10-22). No stage-13d code exists. It is
stage 13c (corrected off-state model, rounded device r = 0.1 µm, the same
K cuts, τ and dense set) with two changes, both prompted by 13c:

1. **Gate-oxide field constraint.** Let $E_{ox}$ be a smooth maximum
   (log-sum-exp, fixed $\kappa$) of the field magnitude on the oxide edges
   adjacent to the SiC/SiO₂ interface. It is evaluated at a **fixed
   blocking voltage $V_{blk}=400$ V** (below the initial true BV of
   517 V). The constraint is $E_{ox}(\theta)\le E_{ox}(\theta_0)$: the
   design may not worsen oxide stress. The gradient comes from the adjoint
   of the corrected off-state electrostatics at fixed $V_{blk}$. The
   absolute value is **reported** against 3 MV/cm, a commonly used design
   guideline; that guideline is unsourced and is not gated.
2. **Defined local field measures** replace the ill-posed global-peak
   criteria. They are smooth maxima of the SiC edge-field magnitude in
   three fixed windows, recorded before running: (W1) JFET/body corner,
   (W2) under the source edge, (W3) the JFET-centre surface.

The problem is to maximize the soft-min BV (unchanged) subject to
$R_{on}\le R_{on,0}$ and the oxide constraint, with SLSQP, plus three
looser $R_{on}$ levels.

## Acceptance criteria

| ID | Criterion | Gating |
|----|-----------|--------|
| D1 | objective-mesh adequacy: the dense true-minimum BV (≥ 40 cuts) of the initial design differs by ≤ 3% between $h$ and $h/2$ | yes |
| D2 | adjoint vs Richardson FD (ε = 1e-3, ε/2) at the initial design: soft-min BV ≤ 1e-4, $R_{on}$ ≤ 1e-5, $E_{ox}$ ≤ 1e-4 relative | yes |
| D3 | final design: the same three checks in absolute form, $\le10^{-6}\max(1,\lVert\nabla\rVert)$, in units of the initial values | yes |
| D4 | feasibility: $R_{on}\le R_{on,0}(1+10^{-8})$, $E_{ox}\le E_{ox,0}(1+10^{-8})$, bounds respected | yes |
| D5 | the dense true-minimum BV improves on $h$ and $h/2$, with the $h/2$ improvement ≥ 50% of the $h$ improvement | yes |
| D6 | at $V_{blk}$, **none** of W1–W3 and $E_{ox}$ increases relative to the initial design, on $h$ and $h/2$ (a relative tolerance of 1e-6 absorbs round-off) | yes |
| D7 | SLSQP exit 0; every evaluation converged | yes |
| R1 | local window fields on $h,h/2,h/4$ with observed order and Richardson (mesh-adequacy information) | reported |
| R2 | BV, $R_{on}$, $E_{ox}$ (absolute, MV/cm) before and after, against the 3 MV/cm guideline; Pareto set; comparison with 13c | reported |

Stage 13d passes if D1–D7 hold. **Not claimed:** measured-device agreement,
oxide reliability, or a global optimum.
