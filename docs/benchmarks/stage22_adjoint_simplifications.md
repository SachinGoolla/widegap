# Stage 22: how much do common adjoint simplifications cost?

Status: **pre-registered** (2026-10-23). No stage-22 code exists.

**Motivation.** The closest prior SiC adjoint optimization (Brungi & Andrei 2026, read via a fetched summary; to be re-checked) appears to drop how the residual depends on doping through **incomplete ionization and mobility** when forming $\partial F/\partial\theta$. It reports no adjoint-vs-FD check. Stage 17 found that dropping the channel-mobility term $F_\mu\mu_X$ alone gives a gradient error of about 5 % (8 % on the JFET variable). This stage measures each simplification separately, on the punch-through-safe geometry.

## 1. Setup

- **Device and design points.** The stage-20 device with $d_{body}$ = 1.4 µm (stage-20 mesh rule), at θ = 0 and at two seeded random interior points: θ ~ U(−0.1, 0.1)⁵, seeds 0 and 1, clipped to the bounds.
- **Figures of merit:** $R_{on,ch}$ (h/2, stage-16 model), the soft-min BV (h) and $E_{ox}$ at 400 V (h). All use the 13d/17 definitions.
- **Gradient variants:**
  - **E** (exact): the verified adjoint, reference.
  - **A:** drop $F_\mu\mu_X$ (the lagged channel mobility). Applies to R_on only.
  - **B:** drop the doping dependence of the bulk mobility in $\partial F/\partial\theta$ ($\partial\mu_{bulk}/\partial N=0$). Applies to R_on only.
  - **C:** replace the incomplete-ionization doping derivative with complete ionization ($\partial N^+/\partial N_D:=1$ and $\partial N^-/\partial N_A:=1$) in $\partial F/\partial\theta$. Applies to all three figures of merit.
  - **D:** A + B + C combined.

## 2. Pre-registered prediction (falsifiable)

1. **P22a:** for at least one figure of merit, the D gradient differs from E by ≥ 1 % in relative norm at θ = 0.
2. **P22b:** an SLSQP run of the 17b-type problem on this geometry, started from θ = 0, using **D gradients** with exact function values, ends at a design that differs from the run with E gradients by at least one of:
   - (i) dense true-minimum BV differing by ≥ 1 V;
   - (ii) a constraint violated by > 1e-6 relative at its final design;
   - (iii) a KKT residual, evaluated with E gradients, ≥ 10 % of $\lVert\nabla BV\rVert$;
   - (iv) SLSQP exit status ≠ 0.

   The E-gradient run is stage 21's main run, reused if stage 21 is recorded first; otherwise it is run here with identical settings.

Either outcome is reported. A refuted prediction is a valid finding: it would show that these simplifications do not matter here.

## 3. Criteria (methodology)

| ID | Criterion | Gating |
|----|-----------|--------|
| Q1 | E at each design point matches Richardson FD (ε = 1e-3 and 5e-4) to ≤ 1e-5 relative norm for every figure of merit | yes |
| Q2 | every variant gradient is finite, deterministic (identical on a re-run) and assembled from the same converged state as E | yes |
| Q3 | every solve converged; no "On entry to" BLAS lines | yes |
| R1 | per variant, figure of merit and design point: relative norm error, cosine to E, per-component error | reported |
| R2 | P22a and P22b outcomes, with the D-run trajectory, final θ, BV, constraint values, KKT residual and exit status | reported |

Stage 22 passes if Q1–Q3 hold.

**Not claimed:** that Brungi & Andrei's results are wrong. Their device, model and objective differ, and their derivation is to be re-read from the PDF.
