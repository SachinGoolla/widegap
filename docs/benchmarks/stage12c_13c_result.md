# Stages 12c and 13c results (2026-10-22): artifact fixed; gates FAILED on criterion definitions and on a newly found design issue

Raw data: `results/stage12c/` and `results/stage13c/` (commit `8f4fa5a`).
The 13c manifest shows `git_dirty: true` only because of the freshly written,
untracked 12c results folder. Implemented by a Sonnet subagent from the
pre-registered spec, then reviewed.

## Stage 12c: the corrected off-state model

The model is $\phi_n=V_D\,w$, with $w$ the harmonic extension across the
p-body.

| ID | Result | Value |
|---|---|---|
| C1 "corner peak field" convergence | **fail** | global peak 2.536 / 2.545 / 2.554 MV/cm (variation 0.7%, against 19% before). The order (−0.06) is ill-posed, see below |
| C2 insensitivity to the interpolation | pass | 0 relative (bit-identical; n < 1e-50 where the two $w$ fields differ) |
| C3 carrier self-consistency | pass | n/N ≤ 2e-19 |
| C4 planar BV vs stage 12 (512.27 V) | **fail** | 556.46 V (+8.6%), see below |
| C5 1D reduction | pass | 1.3e-16 |
| C6 Jacobian | pass | 2.4e-10 |

- **The artifact is removed.** The corner field no longer runs away with
  refinement, and a local corner-window measure converges (order 1.07).
  That measure is 2.5% from Richardson, so it would still miss a 2% gate,
  which is why C1 is recorded as a failure.
- **C1 is ill-posed as written.** After the fix the global peak moved to a
  broad maximum under the source edge. A maximum over edge midpoints
  jitters there, so an observed order is meaningless. The criterion
  should have defined a local measure.
- **C4 shows the stage-12 planar BV was also biased.** The same artifact
  (an artificial electron reservoir at the p side of the junction) sat at
  the *planar* junction too: the old model reproduces 512.27 V exactly.
  The corrected 556.46 V is within 3.4% of an independent 1D
  flat-quasi-Fermi bound through the same column (576.16 V); the rest is
  2D fringing. The old value was 11.1% below that bound. "Mesh-stable" did
  not mean "correct". The earlier caveat in `stage12_result.md` is
  corrected accordingly.

## Stage 13c: 13b re-run, unchanged except the off-state model

| ID | Result | Value |
|---|---|---|
| U1 | pass | 1.26% from Richardson; the observed order is not meaningful (same issue as C1) |
| U2, U3 adjoint vs Richardson FD | pass | 1.7e-8 / 1.1e-9; final-design checks pass |
| U4 feasibility | pass | $R/R_0-1=2\times10^{-12}$ |
| U5 dense true-minimum BV | pass | **516.9 → 618.8 V** on $h$, 518.1 → 622.5 V on $h/2$, at unchanged $R_{on}$ |
| U6 no higher global peak field | **fail** | 2.94 → 3.18 MV/cm on $h$ (2.95 → 3.19 on $h/2$), at the JFET-centre surface under the gate oxide |
| U7 SLSQP | pass | exit 0, 17 iterations |

**The new blind spot is a real design issue.** The optimizer raised the
SiC surface field at the JFET centre, directly under the gate oxide. The
BV there (about 760 V) stays far above the device minimum, so the soft-min
objective cannot see it. In SiC MOSFETs, though, the **gate-oxide field**
(about $\varepsilon_{SiC}/\varepsilon_{ox}\approx2.6\times$ the SiC surface
field) is a known reliability limit. U6 caught a problem that a breakdown-only
objective misses.

**Comparison with 13b** (artifact model): the initial true-minimum BV is
517 V (was 395 V), and the limiting location moved from the corner
artifact to under the source edge. The stage-13 "optimum" evaluated on
the corrected rounded device has a true-minimum BV of only 181 V, which
confirms stage 13 was gamed.

## Lessons

1. A mesh-stable result can still be biased; independent reference
   bounds are needed.
2. A "peak" criterion must define a local, well-posed measure.
3. Device optimization must constrain every failure mode, not only
   avalanche breakdown (here, the oxide field).

## Proposed next (13d)

Add a gate-oxide-field constraint, $E_{ox}\le E_{ox,max}$ at the off-state
operating voltage, with $E_{ox,max}$ either sourced or a declared design
value, to the 13c problem. Replace the global-peak criteria with defined
local measures. Re-run.
