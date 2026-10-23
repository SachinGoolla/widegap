# Stage 21: guarded re-optimization on the punch-through-safe geometry

Status: **pre-registered** (2026-10-23). No stage-21 code exists. Every
optimized design from stages 13–17b sits on the 0.8 µm-body geometry, which
punches through at about 90–130 V (stage 19). Stage 20 found that $d_{body}$ =
1.4 µm is punch-through-safe ($V_{pt,es}/BV\approx1.6$) at +1.4 % R_on. This
stage repeats the stage-17b problem on that geometry, with every guard active.
It is the first design result in this project that does not rest on a
geometry that punches through.

## 1. Problem

- **Device:** the stage-20 design at $d_{body}$ = 1.4 µm, using the stage-20 mesh rule. Everything else is the 13d/17b initial design.
- **Variables and bounds:** the 13d design variables and bounds (4 drift Legendre coefficients and the JFET level).
- **Objective:** maximize the 13d soft-min BV (K cuts, τ; on h). The cuts are those of 13d, mapped to the same x positions.
- **Constraints:**
  - $R_{on,ch}(\theta)\le R_{on,ch}(0)$ (h/2, stage-16 model);
  - $E_{ox}(\theta)\le E_{ox}(0)$ at 400 V (h);
  - the 17b classification constraints ($m$ = 0.02, on h/2 and h/4);
  - **new:** the punch-through guard $V_{pt,es}(\theta)\ge1.2\,BV_{obj}(\theta)$, checked at the final design and at every SLSQP iterate (callback). It is not an SLSQP constraint. If it is violated at any iterate, the stage records that and fails P-guard (below).
- **Solver and Pareto levels:** SLSQP with the 13d settings. Pareto at $R/R_0$ = 1.05, 1.1 and 1.2.

## 2. Criteria

| ID | Criterion | Gating |
|----|-----------|--------|
| S1 | gradient checks at θ = 0 and at the final design, as 17b J1/J2 (same thresholds) | yes |
| S2 | feasibility (R_on,ch, E_ox, classification) as 17b J3/J4b | yes |
| S3 | the dense true-minimum BV improves on h and h/2, with the h/2 improvement ≥ 50 % of the h improvement | yes |
| S4 | **P-guard:** $V_{pt,es}\ge1.2\,BV$ at every iterate, the final design and every Pareto design (h); also checked on h/2 for the final design | yes |
| S5 | field-line BV (15b method) of the final design within 3 % of its dense-cut BV, on h and h/2 | yes |
| S6 | SLSQP exit 0; all solves converged; no "On entry to" BLAS lines | yes |
| R1 | final vs initial: BV (dense, field-line), R_on,ch on h/2 and h/4, the accumulation bracket (18b L′) at initial and final, V_T, E_ox, and the $V_{pt,es}$ margin | reported |
| R2 | comparison with the 17b result on the 0.8 µm geometry (BV gain, θ), labelled clearly as a different geometry | reported |

Stage 21 passes if S1–S6 hold.

**Caveats carried over:**
- the stage-16 channel-model caveats;
- the h/2 on-state mesh uncertainty of about 1.3 % (16 H5) and the 18b B3 narrow failure;
- punch-through safety above 400 V rests on the electrostatic barrier criterion (stage 20).
