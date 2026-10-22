# Stage 10 result: gate PASSED (2026-10-21)

Raw data: `results/stage10/` (commit `a151495`, clean tree). Implemented by a
Sonnet subagent from the pre-registered spec, then reviewed. O3 was amended
before the official run and two deviations were disclosed (spec §3). The
official numbers are identical to the development run.

| ID | Result | Value |
|---|---|---|
| O1 $dI_D/d\theta$ vs FD (initial design) | pass | 6.0e-8 relative |
| O2 $dBV_s/d\theta$ vs FD (initial design) | pass | 3.8e-7 relative |
| O3′ final design vs Richardson FD | pass | $I_D$ 1.1e-10, $BV_s$ 3.8e-5 (tolerance 3.3e-4). Plain FD at ε = 1e-3 differs by 0.30, which is truncation error |
| O4 feasibility | pass | $I_D/I_{D,0}-1=3.3\times10^{-12}$; bounds respected |
| O5 independent re-verification | pass | exact BV +84.6 $V_T$ (233.4 → 318.0) on $h$, +70.3 on $h/2$ (ratio 0.83); $I_D$ differs 0.14% between meshes |
| O6 nested family | pass | 308.4 vs best feasible uniform design 227.1 |
| O7 SLSQP | pass | exit 0, 22 iterations |

**Outcome (study parameters, field-based BV proxy as in ZA19/BA26).** At an
unchanged on-current, the optimized drift and JFET doping raises the
field-based breakdown voltage by **36%** (exact, non-smoothed, mesh $h$).

Pareto set (R1), $BV_s$ vs $I_D/I_{D,0}$:

| $I_D/I_{D,0}$ | 1.00 | 0.98 | 0.95 | 0.90 |
|---|---|---|---|---|
| $BV_s$ | 308.4 | 315.2 | 325.6 | 344.9 |

**Sensitivity maps (R2, `figures/stage10_sensitivity_maps.png`).**
$\gamma_{BV}>0$ in the drift directly below the p-body junction, so
lowering $N_D$ there raises BV. It is negative in the body and the JFET
column. $\gamma_{I_D}<0$ along the JFET and channel, so donors raise the
current.

**Comparison with the published trends (R3, qualitative only).** The
optimizer lowers the drift doping at the body junction (16 → 5.1) and adds
a doping peak about 2 $L_D$ below it. The JFET doping rises only slightly
(18 → 19.5), because the JFET surface field becomes the competing maximum.
The published direction is reproduced only qualitatively and weakly. There
is **no quantitative comparison** with Brungi & Andrei (2026), whose
parameters are unpublished.

**Not claimed:**
- a 4H-SiC MOSFET result: oxide, $E_{crit}$ and doping levels are study
  parameters, and the BV is the field-based proxy, not avalanche;
- that the optimum is global;
- fabrication feasibility.
