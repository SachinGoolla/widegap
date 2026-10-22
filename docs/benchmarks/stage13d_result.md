# Stage 13d result: gate PASSED (D1-D7), with a large absolute oxide field disclosed (2026-10-22)

Raw data: `results/stage13d/` (code commit `8769780`; see provenance). Figures:
`figures/stage13d_pareto.png`, `stage13d_oxide_field.png`, `stage13d_fields.png`,
`stage13d_per_cut_bv.png`. Implemented by a Sonnet subagent from the
pre-registered spec (`stage13d_oxide_constrained.md`), then reviewed.
Label: **new** (an oxide-field-constrained, corner-aware variant of stage 13c).
The adjoint, solver and optimizer are an **independent implementation**
(NumPy/SciPy, no external TCAD reference).

## Verdict

**All of D1-D7 pass (`gate_passed: true`).** Headline: at $R_{on}=R_{on,0}$ and
$E_{ox}(400\,\mathrm{V})=E_{ox,0}$, the dense true-minimum BV improves from
**516.9 to 577.9 V** on $h$ and **518.1 to 580.4 V** on $h/2$ (+61.1 V and +62.3 V).
SLSQP exit status 0 after 12 iterations, 0 failed evaluations.

| ID | Criterion | Measured | Threshold | Result |
|---|---|---|---|---|
| D1 | dense true-min BV of the initial design, $h$ vs $h/2$ (49 cuts) | 516.87 vs 518.09 V, rel. diff 0.236% | <= 3% | pass |
| D2 | adjoint vs Richardson FD at the initial design, worst relative | BV 1.7e-8, $R_{on}$ 1.1e-9, $E_{ox}$ 1.3e-12 (windows W1/W2/W3: 3.3e-11, 1.6e-10, 1.5e-11) | 1e-4 / 1e-5 / 1e-4 | pass |
| D3 | final-design checks, worst absolute (normalized) | BV 6.8e-11, $R_{on}$ 4.2e-11, $E_{ox}$ 4.1e-11 | 1e-6 each | pass |
| D4 | feasibility | $R/R_0-1$ = 2.1e-12; $E_{ox}/E_{ox,0}-1$ = 4.5e-12; bound violation 3.3e-16 | <= 1e-8 | pass |
| D5 | dense true-min BV improves on $h$ and $h/2$; $h/2$ gain >= 50% of $h$ gain | +61.07 V ($h$), +62.32 V ($h/2$), ratio 1.02 | ratio >= 0.5 | pass |
| D6 | none of W1-W3, $E_{ox}$ increases at 400 V ($h$, $h/2$) | see below; all `ok` | rel. tol 1e-6 | pass |
| D7 | SLSQP exit 0, all evaluations converged | status 0, nit 12, 0 failed evals; all three Pareto runs status 0, 0 failures | exit 0 | pass |

D3 also passes in literal units: BV 3.5e-8 V against a 2.8e-4 V tolerance;
$E_{ox}$ 2.1e-10 against 3.2e-6 MV/cm.

D6 detail at $V_{blk}=400$ V, MV/cm, initial to final (relative change):

| Measure | $h$ | $h/2$ |
|---|---|---|
| $E_{ox}$ | 5.2525 to 5.2525 (+4e-12, active constraint) | 5.2867 to 5.2819 (-0.09%) |
| W1 (JFET/body corner) | 2.388 to 2.290 (-4.1%) | 2.447 to 2.345 (-4.1%) |
| W2 (under source edge) | 2.651 to 2.524 (-4.8%) | 2.727 to 2.596 (-4.8%) |
| W3 (JFET-centre surface) | 2.114 to 2.113 (-0.08%) | 2.179 to 2.175 (-0.16%) |

So, unlike 13c (whose U6 failed on the JFET-centre surface field), the
constrained optimum lowers the corner and source-edge fields and leaves W3
essentially unchanged. The limiting cut location is unchanged
(x = 1.70 um on $h$, 1.77 um on $h/2$), and the per-cut BV does not drop at
any dense cut: the smallest per-cut BV gain is +53.8 V on L0 and +54.3 V on
L1 (`worst_cut_BV_drop_V` = -53.8 / -54.3 in the JSON, i.e. a negative drop).

## R1 (reported, not gated): mesh study of local fields

Values at $V_{blk}=400$ V, MV/cm, on $h$, $h/2$, $h/4$ (initial design):

| Measure | $h$ | $h/2$ | $h/4$ | change $h \to h/4$ | observed order | Richardson |
|---|---|---|---|---|---|---|
| $E_{ox}$ | 5.2525 | 5.2867 | 5.3160 | +1.2% | 0.22 | 5.49 |
| W1 | 2.388 | 2.447 | 2.534 | +6.1% | -0.57 | 2.27 |
| W2 | 2.651 | 2.727 | 2.800 | +5.6% | 0.07 | 4.25 |
| W3 | 2.114 | 2.179 | 2.244 | +6.1% | -0.005 | -17.96 |

**The local window fields do not converge; they creep upward with
refinement, monotonically, by roughly 6% from $h$ to $h/4$** (E_ox by 1.2%).
Observed orders are near zero or negative, and the Richardson extrapolations
are not meaningful (W3 extrapolates to a negative field). The final design
behaves the same way (e.g. $E_{ox}$ 5.2525, 5.2819, 5.3116; W1 2.290, 2.345,
2.432). All eight convergence flags in the JSON are `pass: false`. This is
reported only and does not enter the gate, but it means the absolute window
fields are mesh-dependent at the 5-6% level and should be read as relative
comparisons on a fixed mesh (which is how D6 uses them). The
mesh-to-mesh differences for the initial vs final design are consistent
(D6 holds on $h$ and $h/2$).

## R2 (reported, not gated): absolute oxide field, Pareto set, comparison with 13c

**Absolute $E_{ox}$ is far above the design guideline.** At $V_{blk}=400$ V,
$E_{ox}$ = 5.25 MV/cm on $h$ (5.29 on $h/2$, 5.32 on $h/4$), about **1.75x**
the 3 MV/cm design guideline, both before and after. The guideline is
unsourced and is not gated; the constraint only prevents $E_{ox}$ from
getting worse than the initial design.

**At fixed 400 V the constraint holds, but the oxide field at rated BV
rises.** Evaluated at each design's own true-minimum BV on $h$:
$E_{ox}$ = 5.90 MV/cm for the initial design at its own BV (516.9 V),
and 6.20 MV/cm for the final design at its own BV (577.9 V). The
final design at the initial design's BV is 5.87 MV/cm, slightly lower than
the initial 5.90 (JSON: `E_ox_at_initial_true_BV_initial/final`). The
higher rated BV therefore costs a higher oxide field when the device is
actually operated at its new limit; the 400 V constraint does not cover that.
The hard maximum over oxide edges (not the smooth maximum) is 5.2409 to
5.2411 MV/cm, a 0.003% increase; the largest single-edge increase is
0.034 MV/cm at x = 1.07 um.

**Pareto set (dense true-min BV on $h$, $E_{ox}=E_{ox,0}$ in all cases):**

| $R_{on}/R_{on,0}$ | 13d BV (V) | 13c BV (V) | 13d minus 13c (V) | SLSQP nit |
|---|---|---|---|---|
| 1.00 (initial design) | 516.9 | 516.9 | | |
| 1.00 | 577.9 | 618.8 | -40.9 | 12 |
| 1.05 | 624.4 | 642.2 | -17.8 | 10 |
| 1.10 | 651.9 | 660.2 | -8.3 | 15 |
| 1.20 | 684.0 | 685.4 | -1.4 | 16 |

All Pareto runs exit with status 0 and no failed evaluations.
$R_{on}/R_{on,0}$ and $E_{ox}/E_{ox,0}$ are met to ~1e-11. W1 falls monotonically
along the set (2.29, 2.15, 2.06, 1.96 MV/cm) and W3 stays near 2.11.

**Comparison with 13c.** The 13c optimum (same problem without the oxide
constraint) has true-min BV 618.8 V at $R/R_0=1$ but $E_{ox}$ = 7.71 MV/cm
on $h$ (7.75 on $h/2$), **+46.7%** over the initial 5.25 MV/cm, and its
JFET-centre window W3 rises to 3.03 MV/cm. Its all-SiC maximum field also rises
(2.54 to 2.94 MV/cm). The oxide constraint costs **about 41 V of BV at
$R/R_0=1$** (618.8 to 577.9 V); the penalty shrinks as $R_{on}$ is relaxed
and is 1.4 V at $R/R_0=1.2$. This quantifies what the 13c breakdown-only
objective was buying at the expense of the oxide.

**Reported blind-spot scan** (final vs initial design at 400 V on $h$, all SiC
edges, from `blind_spot`): the largest SiC field increase anywhere is
+0.32 MV/cm, located at x = 0, y = -3.29 um (deep in the drift region),
with 666 edges outside the windows showing some increase. The largest
SiC field overall falls (2.54 to 2.40 MV/cm). This is reported, not gated.

## Provenance

- Command `benchmarks/stage13d/run.py`, 2026-10-22T19:07:46 UTC.
- Manifest git commit `8769780`, `git_dirty: true`. The only uncommitted
  content at run time was untracked stage-11 files
  (`src/widegap/pinn_study.py`, `tests/test_pinn_study.py`,
  `benchmarks/stage11/`), which stage 13d does not import. They were
  committed afterward in `28cc75a`; the results were then committed in `f19ccce`.
- Python 3.12.3, widegap 0.2.0, NumPy 2.5.4, SciPy 1.18.1, torch 2.14.1+cpu
  (1 torch thread), AMD Ryzen 7 5700U, 16 logical CPUs, Linux.
- Manifest seeds: 0, 1, 2. The manifest does not record how the seeds are used;
  the pipeline here is deterministic gradient-based optimization from a fixed
  initial design.
- Total runtime 250.6 s (SLSQP main 42.5 s, Pareto runs 107.9 s in parallel,
  adjoint checks 38.6 s initial and 24.2 s final).

## Not claimed

- No agreement with measured devices; no external TCAD reference.
- No oxide **reliability** claim. $E_{ox}$ is 1.75x the unsourced 3 MV/cm
  guideline before and after, and rises to about 6.2 MV/cm at the new rated BV.
  The constraint prevents worsening at 400 V; it does not make the oxide safe.
- Not a global optimum: SLSQP finds a local solution from one starting
  design in a 5-parameter family.
- BV uses a vertical-cut ionization-integral approximation of a 2D device
  at K cuts; the dense true minimum (>= 40 cuts) is the check, but it is
  still not a full 2D avalanche solution.
- The reference voltage $V_{blk}=400$ V is a choice (below the initial true BV
  of 517 V); a different reference would give a different constraint and
  different Pareto numbers.
- Local window fields and $E_{ox}$ are not mesh-converged (R1): they rise
  by ~1-6% from $h$ to $h/4$, so absolute values are mesh-dependent.
