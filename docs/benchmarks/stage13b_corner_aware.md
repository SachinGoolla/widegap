# Stage 13b: corner-aware BV objective with a rounded body junction

Status: **pre-registered** (2026-10-22). No stage-13b code exists. This
fixes the stage-13 specification flaw: the fixed planar-cut BV was gamed
while the JFET/body corner BV fell (see `stage13_result.md`).

## 1. Changes from stage 13 (everything else unchanged)

- **Rounded body junction.** The p-body acceptor profile has smooth lateral
  and vertical edges: a product of erfc profiles with a lateral
  straggle/rounding length $r$ (study value, recorded; for example 0.1 µm).
  The same treatment applies to the n⁺ source. This removes the
  geometric field singularity, so the corner field is mesh-resolvable.
- **BV objective.**
  $BV_{obj}=-\tau\log\sum_{k=1}^{K}e^{-BV_k/\tau}$, a smooth minimum over
  $K$ fixed vertical cuts. The cuts span from the symmetry line through the
  JFET/body corner region to the planar region under the body. $K$, the
  positions and $\tau$ (in volts) are chosen and recorded before running.
  $BV_k$ is the avalanche BV (Niwa ionization integral) along cut $k$.
  Gradients are the stage-13 implicit-function adjoint per cut, combined
  with soft-min weights.
- **Problem.** Maximize $BV_{obj}$ subject to $R_{on,sp}\le R_{on,0}$
  (SLSQP), plus three looser levels.

## 2. Acceptance criteria

| ID | Criterion | Gating |
|----|-----------|--------|
| U1 | the corner peak field converges with mesh refinement ($h,h/2,h/4$, initial design): finest within 2% of Richardson | yes |
| U2 | $dBV_{obj}/d\theta$ and $dR_{on}/d\theta$ vs Richardson FD (ε = 1e-3), initial design: ≤ 1e-4 and ≤ 1e-5 relative | yes |
| U3 | final design: the same checks in absolute form, $\le10^{-6}\max(1,\lVert\nabla\rVert)$, with $R_{on}$ in units of $R_{on,0}$ | yes |
| U4 | feasibility, as stage-13 T4 | yes |
| U5 | **no blind spot**: the *true* minimum BV over a dense set of $\ge40$ vertical cuts (independent of the objective's cuts, covering the whole semiconductor width), from fresh solves on $h$ and $h/2$, improves over the initial design on both meshes. The $h/2$ improvement must be ≥ 50% of the $h$ improvement | yes |
| U6 | the global peak field of the final design is not higher than the initial design's at the same drain voltage (initial BV), on $h$ and $h/2$ | yes |
| U7 | SLSQP exit 0, every evaluation converged | yes |
| R1 | per-cut BV profiles before and after; which cut limits BV; Pareto set vs the unipolar limit | reported |
| R2 | sensitivity maps; optimized profile; comparison with stage 13 | reported |

Stage 13b passes if U1–U7 hold. **Not claimed:** the full 2D
ionization integral along field lines, since vertical cuts approximate it,
especially near the curved junction; agreement with measured devices; a
global optimum.
