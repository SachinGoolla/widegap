# Stage 17: channel-aware re-optimization of the 4H-SiC DMOSFET

Status: **pre-registered** (2026-10-23). No stage-17 code exists. This is the
13d problem with **R_on from the stage-16 channel-mobility model**, using an
exact adjoint through the nonlocal mobility. It covers the 5 design variables
of 13d only. A p-body doping variable (17b) needs a V_T constraint and a
punch-through guard, and will be pre-registered separately.

## 0. Basis and disclosures

- **Stage 16 did not pass.** It failed H5 (R_on changed 3.15 % from h to h/2, then 1.27 % from h/2 to h/4). The remedy here is the mesh choice:
  - the on-state (R_on, its adjoint, the R_on constraint) is evaluated on **h/2**;
  - the off-state (soft-min BV and E_ox) stays on **h**, as in 13d, where D1 showed h is adequate for these;
  - the on- and off-state are separate solves, so using two meshes is consistent;
  - R_on is re-checked independently on **h/4** (J6).
- **Exploratory diagnostics were seen before writing this.** They are scratchpad only, labelled exploratory, and used the stub-verified Noguchi entries before they were merged.
  - At the 13d initial and final designs on h/2, the channel model leaves the absolute ∂R_on/∂θ for the four drift coefficients within about 2 % of the bulk model.
  - It multiplies the JFET-level component by 2.2–2.4. The mechanism is not yet identified.
  - A non-negative least-squares stationarity check at the 13d final, with the R_on constraint from the channel model, leaves a 29 % residual, against 8 % with the bulk model.
  - These diagnostics set the scope (the 5-variable problem is not trivially at its optimum). They did not set any threshold.
  - A bracket for the missing accumulation-layer mobility on the n-type surface under the gate (fixed 20–50 cm²/Vs) changes R_on by +2.0 % to +6.2 %. This is reported as a caveat, since there is no sourced model.

## 1. Problem

Maximize the 13d soft-min BV (same K cuts and τ, on h) subject to:
- $R_{on,ch}(\theta)\le R_{on,ch}(0)$ on h/2;
- $E_{ox}(\theta)\le E_{ox}(0)$ at $V_{blk}$ = 400 V on h;
- the 13d bounds.

The solver is SLSQP. Pareto runs at $R/R_0$ = 1.05, 1.1 and 1.2 follow, as in 13d.

## 2. Exact adjoint of the channel-model R_on

- The forward solve is the stage-16 Picard loop (unchanged). At the converged fixed point, $F(X,\mu(X,\theta),\theta)=0$.
- The adjoint is $(F_X+F_\mu\,\mu_X)^T\lambda=J_X$, and $dJ/d\theta=-\lambda^T(F_\theta+F_\mu\,\mu_\theta)$.
- $\mu_X$ is sparse: each channel edge's $\mu$ depends on $N_S$ of its two columns, and so on $n$ at the nodes within $d_{ch}$ in those columns. It is analytic, through the Noguchi formulas, the harmonic edge mean and the trapezoidal $N_S$.
- $\mu_\theta$ comes in through doping-dependent bulk terms only.
  - The 5 variables change $N_D$ in the drift and JFET regions.
  - $N_{dep}$ uses the surface $N_A$, which the 5 variables do not change.
  - The p-surface column classification is frozen at θ = 0 (J4).
- The lagged gradient (without $F_\mu\mu_X$) is computed and its error **reported**.

## 3. Acceptance criteria

| ID | Criterion | Gating |
|----|-----------|--------|
| J1 | exact adjoint $dR_{on,ch}/d\theta$ vs Richardson FD (ε = 1e-3 and 5e-4, all 5 unit directions) at θ = 0 on h/2: $\lVert g_{adj}-g_{FD}\rVert\le10^{-5}\lVert g_{FD}\rVert$. The diagnostics put the FD noise floor at ≤ 1e-6 relative. The soft-min BV and E_ox gradients are re-checked as in 13d D2 | yes |
| J2 | at the final design, the same three gradients in absolute form, $\le10^{-6}\max(1,\lVert\nabla\rVert)$ in units of the initial values | yes |
| J3 | feasibility: $R_{on,ch}\le R_{0}(1+10^{-8})$ (h/2), $E_{ox}\le E_{ox,0}(1+10^{-8})$, bounds respected | yes |
| J4 | smoothness guards at every evaluated θ: no surface column changes the sign of $N_D-N_A$ relative to θ = 0; no gated channel column has $N_S$ within a factor 2 of the $10^{11}$ clamp or $E_{eff}$ within a factor 2 of the $10^{-2}$ MV/cm clamp | yes |
| J5 | the dense true-minimum BV (≥ 40 cuts) improves on h and h/2, with the h/2 improvement ≥ 50 % of the h improvement | yes |
| J6 | independent on-state re-check on **h/4**: $R_{on,ch}$(final) ≤ $R_{on,ch}$(initial)·(1 + 1e-3) | yes |
| J7 | SLSQP exit 0; every Picard and Newton solve converged; `widegap.linalg.nonfinite_calls` reported. Non-finite calls are allowed only inside rejected Newton trial steps, and the log has no "On entry to" BLAS lines | yes |
| R1 | final vs the 13d final: BV gain, θ (especially the JFET level), the channel-model R_on of the 13d final, and the Pareto set | reported |
| R2 | lagged-gradient error; the JFET-level sensitivity (which edges carry it: map of $\lambda$-weighted $\partial F/\partial\theta_5$) | reported |
| R3 | W1–W3 and E_ox at 400 V on h and h/2 (W measures are not mesh-resolvable per 15d, so same-mesh comparison only); field-line BV (15b method) of the final design | reported |

Stage 17 passes if J1–J7 hold.

**Caveats carried over:**
- all stage-16 caveats: Hall factor 1, NO anneal, no $N_{it}$ dependence, 300 K;
- no accumulation-layer mobility (bracket above);
- the h/2 on-state mesh uncertainty of about 1.3 %;
- the stage-12c zero-current off-state.

**Not claimed:** measured-device agreement, global optimality.
