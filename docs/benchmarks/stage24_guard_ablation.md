# Stage 24 (experiment A): prospective test — does an optimizer without a punch-through guard trade blocking capability for avalanche BV?

Status: **pre-registered** (2026-10-23). No stage-24 code exists.

**Motivation.**
- Stage 23 found, **post hoc**, that the optimized 0.8 µm designs 13d and 17b raised avalanche BV by 61–56 V while their blocking voltage $\min(BV,V_{pt})$ **fell** from 94.8 V to 83.2 and 80.9 V.
- This stage tests the underlying claim **prospectively**: an adjoint optimizer whose model cannot see a failure mode moves toward it. Both the prediction and the thresholds are fixed here.

## 1. Device and arms

- **Device.** The stage-20 design at $d_{body}$ = 1.2 µm. Its punch-through margin is thin: $V_{pt,es}/BV$ = 1.057 (h) and 1.074 (h/2) at θ = 0, per stage 20. The stage-20 mesh rule applies.
- **Common to both arms.** The stage-21 problem:
  - 5 variables;
  - maximize the soft-min BV (h);
  - $R_{on,ch}\le R_0$ (h/2), $E_{ox}\le E_{ox,0}$ (400 V, h), and the 17b classification constraints;
  - SLSQP with the 13d settings.
- **Arm U (unguarded):** no punch-through information.
- **Arm G (guarded):** one additional SLSQP constraint, $B_\tau(\theta;V^*)\ge10\,V_T$ at the fixed voltage $V^*$ = 650 V.
  - $B_\tau$ is a smooth minimum over the stage-20 barrier columns. It nests log-sum-exp soft-minima with temperature τ = 0.5 $V_T$: one over $y$ (for $\min_y\psi$) and one over columns.
  - The barrier is $\psi_{surface}-\min_y\psi$, evaluated on the 12c off-state at $V_D$ = 650 V on h.
  - The gradient comes from the off-state adjoint at fixed $V_D$, the same structure as the 13d $E_{ox}$ constraint.
  - $V^*$ = 650 V is chosen above every BV this design family reached (the stage-21 optimum was 608 V).
  - At θ = 0 the guard must hold: $V_{pt,es}$ is 578–586 V < 650 V, so it may be **infeasible at the start**. If it is, Arm G starts from the feasible point found by a short pre-registered restoration: SLSQP minimizing $-B_\tau$ subject to the common constraints. The start point is recorded.

## 2. Pre-registered predictions

1. **P24a (Arm U):** at the Arm-U final design, $V_{pt,es}$(h) **< the dense true-minimum BV** (h). The optimized design is punch-through-limited, so its blocking voltage is below its avalanche BV.
2. **P24b (Arm U vs initial):** the Arm-U final blocking voltage $\min(BV_{dense},V_{pt,es})$ (h) is **less than 1.05×** the initial design's. That means less than 5 % of real blocking gain, whatever the avalanche-BV gain.
3. **P24c (Arm G):** the Arm-G final has $V_{pt,es}$(h) ≥ $BV_{dense}$(h). Its dense-BV gain is smaller than Arm U's by **≥ 5 V**, and its blocking voltage is **higher** than Arm U's.

Each prediction is reported as confirmed or refuted. Refutation is a valid result. It would mean that, on this geometry, the optimizer does not drift toward punch-through.

## 3. Criteria (methodology)

| ID | Criterion | Gating |
|----|-----------|--------|
| X1 | the $B_\tau$ adjoint gradient vs Richardson FD (ε = 1e-3 and 5e-4) at θ = 0 and at the Arm-G final: ≤ 1e-4 relative norm. The other gradients are checked as in stage 21 S1 | yes |
| X2 | feasibility of each arm's own constraints at its final design, as in stage 21 S2; plus $B_\tau\ge10V_T(1-10^{-8})$ for Arm G | yes |
| X3 | SLSQP exit 0 for both arms (and for the restoration, if used); all solves converged; no "On entry to" BLAS lines | yes |
| X4 | at each final design, $V_{pt,es}$ changes by ≤ 5 % and the dense BV by ≤ 3 % between h and h/2. P24a–c are also evaluated on h/2 and reported | yes |
| R1 | per arm: θ, dense BV (h, h/2), field-line BV, $V_{pt,es}$ (h, h/2), the hard-min barrier at 650 V, the blocking voltage, R_on,ch, E_ox and the active constraints | reported |
| R2 | P24a–c outcomes | reported |

Stage 24 passes if X1–X4 hold. The predictions are reported whatever the gate outcome.

**Not claimed:**
- generality beyond this geometry and model family;
- the 1D DD onset of the finals is reported in R1 if numerics allow, but it is not gated.
