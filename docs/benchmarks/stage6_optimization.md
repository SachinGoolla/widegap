# Stage 6: constrained doping optimization (1D, field-based BV at fixed drift resistance)

Status: **pre-registered** (2026-10-20). No stage-6 code exists at the time
of this commit. The FD error budget and the treatment of the
non-differentiable maximum were decided analytically beforehand.

## 1. Problem

The classic breakdown / on-resistance trade-off (Baliga, *Fundamentals of
Power Semiconductor Devices*, Springer 2008, the ZA19 reference [1])
applied to a punch-through diode:

- Device (scaled units, $\delta=10^{-3}$): p⁺ ($a=100$) on $[-3,0]$, n
  drift with design doping $b(x)$ on $[0,W]$, $W=12$, n⁺ ($b=100$) on
  $[W,W+3]$. Uniform mesh $h=0.025$. $E_{crit}=15$, a study constant as in
  stage 4.
- **Objective:** maximize the field-based BV (stage 4, zero-current
  electrostatics, no impact ionization). It is a **proxy**, not avalanche
  breakdown.
- **Constraint:** drift-resistance proxy $R[b]=\int_0^W dx/b(x)=R_0=W/b_0$,
  with $b_0=1$ (the initial uniform design). $R$ is the exact ohmic
  resistance of a neutral drift layer at constant mobility. It ignores
  depletion and injection, so it too is a proxy.
- **Parametrization (smoothness, positivity, exact constraint):**
  $g(\xi)=\sum_{k=1}^{5}\theta_kP_k(2\xi-1)$ with $\xi=x/W$, Legendre
  $P_k$, and $b(x)=e^{g}\int_0^We^{-g}dx\,/R_0$. Then $R[b]=R_0$ holds
  identically, and $\theta=0$ is the initial uniform design. The doping
  enters the forward model as control-volume averages over the true
  piecewise doping (the stage-2 convention, generalised to smooth $b$ by
  quadrature).
- **Bounds:** $0.05\le b\le20$ at every drift node. These are enforced by
  rejecting trial steps that violate them.

## 2. Smoothed breakdown index

$B_s=\kappa^{-1}\log\sum_e\exp(\kappa B_e)$ with $B_e=|E_e|/E_{crit}$ and
$\kappa=200$. It is smooth, and $B\le B_s\le B+\log(N_e)/\kappa$
($\approx B+0.033$ here), so $BV_s\le BV$ is conservative. $BV_s$ solves
$B_s=1$ (stage-4 bracketing plus Brent, `xtol` $10^{-12}$). Its gradient
is $dBV_s/dC=-(\partial B_s/\partial C)/(\partial B_s/\partial V_R)$ by the
stage-5 adjoint, and $dBV_s/d\theta$ follows by the chain rule through the
parametrization.

## 3. Optimizer

Gradient ascent with BFGS curvature. Armijo backtracking (factor 1/2,
$c_1=10^{-4}$, at most 30 halvings) accepts a step only if it increases
$BV_s$ **and** satisfies the bounds. Stop when the relative improvement over
the last 3 accepted iterations is below $10^{-8}$, when
$\lVert\nabla_\theta BV_s\rVert\le10^{-6}BV_s$, or when the line search
fails. At most 60 iterations.

## 4. Acceptance criteria

FD budget: $\epsilon=10^{-4}$ in $\theta$, Brent noise
$10^{-12}/\epsilon=10^{-8}$, truncation about $10^{-8}$. Threshold: 1e-5.

| ID | Criterion | Gating |
|----|-----------|--------|
| S6a | $dBV_s/d\theta$ (adjoint plus chain rule) vs central FD along 3 seeded random unit directions, at the initial and the final design: relative error ≤ 1e-5 | yes |
| S6b | $BV_s$ strictly increasing over accepted iterates | yes |
| S6c | at every accepted iterate, $|R/R_0-1|\le10^{-12}$ and $0.05\le b\le20$ at all drift nodes | yes |
| S6d | the optimizer stops by a rule of §3 within 60 iterations, and the stop reason is recorded | yes |
| S6e | **independent re-verification** with exact (non-smoothed) BV, fresh stage-4 sweeps and new `Diode1D` instances. On $h=0.025$ and $h=0.0125$: $BV(\text{final})>BV(\text{initial})$, and the improvement on $h=0.0125$ is at least 50% of that on $h=0.025$ | yes |
| S6f | nested-family check: $BV_s(\text{final})\ge\max$ over a 41-point grid of linear gradings ($\theta_1\in[-2,2]$, others 0, bounds respected) minus $10^{-6}BV_s$ | yes |
| R6a | optimized profile, field profiles, and $BV$, $BV_s$ before and after | reported |
| R6b | final $\lVert\nabla_\theta BV_s\rVert$ and $\theta$ | reported |

## 5. Gate and claims

Stage 6 passes if S6a–S6f hold. Even then, the optimized profile is **not**
a fabrication-ready design and **not** a claim about real breakdown. It
demonstrates verified constrained optimization on verified 1D proxies.
