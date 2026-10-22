# Stage 13: full-physics adjoint doping optimization of a 4H-SiC DMOSFET

Status: **pre-registered** (2026-10-21). No stage-13 code exists. This
combines the stage-10 method with the stage-12 sourced physics: incomplete
ionization, Caughey–Thomas field mobility, interface charge, avalanche BV
via the ionization integral, and registry parameters
($\varepsilon_{ox}$ = 3.9 owner-approved).

## 1. Problem

- **Device.** The stage-12 12b DMOSFET, in physical units.
- **Design variables.** Drift donor profile
  $N_D(y)=N_{D0}\exp(\sum_{k=0}^{3}\theta_kP_k(\eta))$ and the JFET donor
  level $\theta_J$ (log scale). Bounds: drift $N_D\in[10^{15},10^{17}]$
  cm⁻³, JFET $\in[10^{16},5\cdot10^{17}]$.
- **On-state objective.** Specific on-resistance $R_{on,sp}$ at
  $V_{GS}$ = 15 V, $V_{DS}$ = 0.1 V, from the full stage-9a DD with
  stage-12 parameters.
- **Off-state objective.** Avalanche BV (Niwa coefficients) by the
  ionization integral along a **fixed** vertical planar cut, at the
  stage-12 $x$ position under the body contact, on zero-current 2D
  electrostatics with incomplete ionization (stage-12 machinery).
- **Problem.** Maximize BV subject to $R_{on,sp}\le R_{on,sp,0}$ (the
  initial value). Solved with SLSQP on adjoint gradients, plus three
  looser $R_{on}$ levels for a Pareto set.

## 2. Gradients

- **$R_{on}$.** Adjoint of the full DD system,
  $dJ/d\theta=-\lambda^T\partial F/\partial\theta$. $\partial F/\partial\theta$
  includes doping entering (i) Poisson through the *ionized* fractions and
  (ii) the continuity fluxes through the doping-dependent mobility. It is
  analytic where implemented. A residual-only central difference in
  $\theta$ (no re-solve) is allowed for $\partial F/\partial\theta$,
  provided its step and accuracy are documented.
- **BV.** Implicit function on $I_n(V_D,C)=1$:
  $dBV/d\theta=-(dI_n/d\theta)/(dI_n/dV_D)$. The total derivatives use the
  adjoint of the electrostatic system, and $\partial I_n/\partial\psi$ along
  the cut is analytic.

## 3. Acceptance criteria

| ID | Criterion | Gating |
|----|-----------|--------|
| T1 | $dR_{on}/d\theta$ vs Richardson-extrapolated central FD ($\epsilon,\epsilon/2$; $\epsilon=10^{-3}$), 3 seeded directions, initial design: ≤ 1e-5 relative | yes |
| T2 | $dBV/d\theta$, same check: ≤ 1e-4 relative | yes |
| T3 | at the final design, both gradients vs Richardson FD: absolute ≤ $10^{-6}\max(1,\lVert\nabla\rVert)$ | yes |
| T4 | feasibility: $R_{on}\le R_{on,0}(1+10^{-8})$; bounds respected | yes |
| T5 | independent re-verification on $h$ and $h/2$ (fresh solves): BV improves on both meshes, with the $h/2$ improvement ≥ 50% of the $h$ improvement; $R_{on}$ at $h/2$ within 3% of $h$ | yes |
| T6 | nested family: final BV ≥ the best uniform-drift design (21-point grid in $\theta_0$, $\theta_J$ fixed) satisfying the constraint, minus $10^{-6}$ BV | yes |
| T7 | SLSQP exit 0 within 60 iterations; every evaluation converged | yes |
| R1 | Pareto set (BV, $R_{on,sp}$) at 4 constraint levels, plotted against the unipolar 4H-SiC limit computed from the same registry parameters (no fitting) | reported |
| R2 | 2D sensitivity maps $\gamma_{BV}$ and $\gamma_{R_{on}}$ (ZA19/BA26 sign convention) | reported |
| R3 | optimized profile vs the published trends (BA26), qualitatively only; all stage-12 flags carried over | reported |

Stage 13 passes if T1–T7 hold. **Not claimed:** agreement with measured
devices, fabrication feasibility, or a global optimum.
