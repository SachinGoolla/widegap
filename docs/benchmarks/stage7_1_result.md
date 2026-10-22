# Stage 7.1 result: gate PASSED on the first official run (2026-10-21)

Raw data: `results/stage7_1/` (commit `f1c09ce`, clean tree). Implemented by a
Sonnet subagent from the pre-registered spec. The flux and identity formulas
were reviewed by the lead.

| ID | Result | Value |
|---|---|---|
| T1 equivalence with stage 3b (δ=1e-3) | pass | max\|Δψ\| 3.0e-14; current 1.6e-11 relative |
| T2 MMS | pass | orders 2.000 ± 0.001 (ψ, n, p) |
| T3 Jacobian vs FD | pass | ≤ 4.7e-7 per-row relative (FD-limited) |
| T4 convergence at δ=1e-24 | pass | all 111 biases, ≤ 6 Newton iterations |
| T5 equilibrium | pass | ψ vs stage 2 1.6e-11; \|φ\| ≤ 1.2e-26 |
| T6 minority currents vs exact SG identity | pass | ≤ 3.9e-14 relative, at every bias (currents down to ~1e-45) |
| T7 total-current constancy, \|j\| ≥ 1e-6 | pass | 3.6e-6 (10 biases, V ≥ 101) |
| T8 d ln j/dV, V=95–102 | pass | 0.996–0.997 |

**R1, the old (ψ,n,p) solver at the same δ.** It reports "convergence" at
low bias while its total current is wrong by up to 9×10¹¹ relative. At
V = 100 it produces **negative hole densities** (≈ −1e-35). This confirms
the stage-3 variables are unusable for realistic SiC.

**R2 (SRH, τ = 10³).** Converges everywhere. d ln j/dV rises from 0.56
to 0.74 over V = 95–102, the expected recombination-dominated regime.

**R3 (illustration only).** Complete ionization, μ_n = μ_p = 1101 cm²/V s.
j ≈ 7e-3 A/cm² at 100 V_T (2.59 V), 1 A/cm² at 105 V_T (2.71 V), and
114 A/cm² at 110 V_T (2.84 V). These are plausible orders of magnitude for a
4H-SiC pn diode's turn-on near 2.7 V, but they are **not** validated
against measurements.

**Resolvability limit (by design, §3 of the spec).** Majority-carrier current
is resolved only for |j| ≳ 1e-6 scaled (about 1e-2 A/cm²). Below that, total
current is reported but is not meaningful. Minority currents are resolved at
all biases.

**Not claimed:** incomplete ionization inside drift–diffusion (stage 7 covers
equilibrium only); field-dependent mobility; impact ionization; 2D.
