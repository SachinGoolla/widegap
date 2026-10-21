# Stage 4b: re-registered criteria with an exact first-integral reference (2026-10-20)

**Disclosures.**
1. These criteria were written after the v1 run failed
   (`docs/benchmarks/stage4_result.md`), so this is not a blind test.
2. The 4b evaluation code (`benchmarks/stage4/run.py --criteria 4b`) was
   committed *before* this document, together with the solver (commit
   "feat(fv): reverse-bias electrostatics…"). That breaks the project's
   "criteria before code" rule. The thresholds below are the ones in that
   code. They are stated here before the 4b run is executed.

The substitute reference is exact rather than approximate. For 1D
zero-current electrostatics with Boltzmann carriers, the first integral
$E(x_0)^2/2=\int_{\psi(x_0)}^{\psi_{cathode}}\rho(\psi)\,d\psi$ holds
exactly. It is evaluated with the true piecewise doping per element and the
numerical nodal $\psi$, whose error is $O(h^2)$.

| ID | Criterion | Gating |
|----|-----------|--------|
| G1 | unchanged: $\lvert B(BV)-1\rvert\le10^{-9}$ on all meshes | yes |
| G2′ | meshes $h\in\{0.05,\dots,0.003125\}$ (halving). The last two observed BV orders lie in [0.8, 1.2] (first order is expected from the element-midpoint maximum). The finest-mesh BV is within 1% of the first-order Richardson extrapolation | yes |
| G3′ | at $V_R\in\{20,50,80\}$, B1 and B2: the junction-element field, Richardson-extrapolated from the two finest meshes, matches the exact first-integral $\lvert E(0)\rvert$ (finest mesh) to ≤ 2e-3 relative | yes |
| G6′ | B1, $V_R=50$: slope of $\lvert E\rvert$ vs $x$ over element midpoints with $x\ge1$ and $n<0.01$ equals −1 within 1% | yes |
| G7, G8, G9 | unchanged | yes |
| — | v1 references (depletion BV, $E^2$ and $E$ slopes) | reported only |

Stage 4 passes if G1, G2′, G3′, G6′, G7, G8 and G9 hold.
