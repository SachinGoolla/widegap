# Stage 7b: re-registered M5(b) geometry (2026-10-21)

**Disclosure.** Written after run 2 failed M5(b)
(`docs/benchmarks/stage7_result.md`). Only the M5 geometry and the bulk
check points change. All parameters, models, criteria and thresholds
(including 1e-10 for M5(b)) are unchanged, and everything is re-run fresh.

The geometry is sized from the screening analysis in the diagnosis. The
p-side screening length is about 0.82 $L_D$, with depletion depth about
1.5 $L_D$. The n side has $L_D$ = 1, with depletion depth about 15. The
domain is $[-60,80]\,L_D$. The bulk is checked at $x=-30$, which is 28.5
$L_D$ beyond the p depletion edge (tail about $e^{-35}$), and at $x=45$,
30 beyond the n edge (tail about $e^{-30}$). Both are far from the
contacts. Mesh $h\in\{0.025,0.0125\}$ as before.

Stage 7 passes if M1–M6 hold with this geometry.
