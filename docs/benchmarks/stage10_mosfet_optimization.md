# Stage 10: 2D adjoint sensitivities and doping optimization of a vertical MOSFET half-cell

Status: **pre-registered** (2026-10-21). No stage-10 code exists. This is an
independent, open implementation of the *method* of Zhu & Andrei (2019) and
Brungi & Andrei (2026): the adjoint doping sensitivities of a field-based
breakdown index and of the on-current, used for doping optimization. It
does **not** reproduce their device results; their parameters are
unpublished. All device constants are study parameters unless drawn from the
stage-7 registry.

## 1. Device and models

- **Device.** Vertical DMOSFET half-cell, 2D, symmetric about $x=0$ (zero
  flux there):
  - n⁺ source and p-body at the top;
  - gate over the channel through oxide (stage-8.1 regions);
  - JFET region, n⁻ drift layer and n⁺ drain at the bottom (drain contact).

  The geometry and mesh are chosen by the implementation and recorded. The
  mesh is moderate so the official run finishes in under about 30 minutes.
- **Physics.** Stage-8.1 quasi-Fermi DD with regions, constant mobility,
  complete ionization, $R\equiv0$, $\delta=10^{-24}$. The stage-9a models
  are not used here, which keeps the adjoint scope verifiable; they are a
  later extension.
- **On-state objective.** $I_D$ at $V_{GS}=V_{on}$ and a small
  $V_{DS}=V_{lin}$.
- **Off-state objective.** Smoothed field-based breakdown voltage $BV_s$.
  It is found by solving the full QF DD at $V_{GS}=0$ with the drain swept,
  using $B_s=\kappa^{-1}\log\sum_e e^{\kappa|E_e|/E_{crit}}$ over
  semiconductor edges (stage-6 smoothing, conservative). $BV_s$ solves
  $B_s=1$.
- **Adjoint.** $F_X^T\lambda=\partial J/\partial X$. Doping enters the
  Poisson rows as $-V_{s,i}C_i$, so $dJ/dC_i=\lambda_{\psi,i}V_{s,i}$. For
  $BV_s$, the implicit-function form is used, with
  $\partial B_s/\partial V_D=-\lambda^T\partial F/\partial V_D$ (Dirichlet
  rows at the drain).
- **Design variables.** Drift-layer donor profile
  $N_D(y)=N_0\exp\big(\sum_{k=0}^{3}\theta_kP_k(\eta)\big)$ over the drift
  depth, plus the JFET-region donor level $\theta_J$ (log scale). Bounds
  are enforced by rejecting steps. Source, body and drain dopings are fixed.
- **Problem.** Maximize $BV_s$ subject to $I_D\ge I_{D,0}$, with
  $I_{D,0}$ = the initial design's $I_D$ (and three looser levels for a
  Pareto set). Solved with SLSQP on adjoint gradients, with every
  evaluation recorded.

## 2. Acceptance criteria

FD budget: $\epsilon=10^{-3}$ in $\theta$; $BV_s$ from Brent with
`xtol` $10^{-10}$; Newton on the per-row floor.

| ID | Criterion | Gating |
|----|-----------|--------|
| O1 | $dI_D/d\theta$ (adjoint) vs central FD along 3 seeded unit directions at the initial design: ≤ 1e-5 relative | yes |
| O2 | $dBV_s/d\theta$ (adjoint plus implicit function) vs central FD, same directions, initial design: ≤ 1e-4 relative | yes |
| O3 | At the final design (main problem): adjoint–FD absolute difference ≤ 1e-6·max(1, ‖∇‖) for both gradients (no relative errors near stationarity; stage-6 lesson) | yes |
| O4 | Final design feasible: $I_D\ge I_{D,0}(1-10^{-8})$ on the optimization mesh; bounds respected | yes |
| O5 | **Independent re-verification** on mesh $h$ and $h/2$, with fresh solves: exact (non-smoothed) $BV$ improves over the initial design on both meshes, with the $h/2$ improvement ≥ 50% of the $h$ improvement; $I_D$ at $h/2$ is within 2% of its $h$ value | yes |
| O6 | **Nested family:** final $BV_s$ ≥ the best uniform-drift-doping design ($\theta_{1..3}=0$, 21-point grid in $\theta_0$, $\theta_J$ fixed) satisfying $I_D\ge I_{D,0}$, minus $10^{-6}BV_s$ | yes |
| O7 | Optimizer terminates successfully (SLSQP exit 0) within 60 iterations | yes |
| R1 | Pareto points ($BV$, $I_D$) for 4 constraint levels | reported |
| R2 | 2D sensitivity maps $\gamma_{BV}(x,y)$ and $\gamma_{I_D}(x,y)$ in the ZA19/BA26 sign convention ($D=N_A-N_D$) for the initial design: the analogue of the papers' signature figures | reported |
| R3 | Qualitative comparison with the published trends (e.g. lower drift doping near the junction raises BV; higher JFET doping raises $I_D$), with no quantitative claim | reported |

Stage 10 passes if O1–O7 hold.

## 3. Amendments and disclosed deviations (before the official run, after a development run)

**O3 (restated).** The development run showed that central FD at
$\epsilon=10^{-3}$ carries an $O(\epsilon^2)$ truncation error of up to
0.30 at the final design. The FD error fell as $\epsilon^2$: 0.30, 0.075
and 0.0014 at $\epsilon=10^{-3}$, $5\times10^{-4}$ and $10^{-4}$. At the
minimax optimum two field maxima nearly tie (the junction and the JFET
surface), which makes $B_s$ strongly curved. The original O3 tolerance,
$10^{-6}\lVert\nabla\rVert$, is unattainable with that step. **O3′:** the
adjoint is compared with the Richardson extrapolation
$(4\,FD(\epsilon/2)-FD(\epsilon))/3$, $\epsilon=10^{-3}$, with the same
tolerance $|\cdot|\le10^{-6}\max(1,\lVert\nabla\rVert)$ for both gradients.
In the development run the adjoint matched it to 3.8e-5 ($BV_s$) and
1.1e-10 ($I_D$), and matched FD at $\epsilon=10^{-5}$ to 4e-5.

**Deviation 1: breakdown index domain.** $B_s$ excludes a box around the
n⁺ source / p-body junction (1609 of 1934 semiconductor edges are kept).
At $\delta=10^{-24}$ that junction's built-in field (about 108
$V_T/L_D$) exceeds the study $E_{crit}=75$ at zero drain bias and does
not depend on $V_D$. Raising $E_{crit}$ above it pushes BV to about 600
$V_T$ or more, where Newton fails. This is an artifact of the study
scaling; the device's breakdown-relevant fields are the drain-dependent
ones, which are kept.

**Deviation 2: bounds.** These are enforced as SLSQP box and linear
constraints, not by step rejection. The maximum violation was 5.6e-15.
