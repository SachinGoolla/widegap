# Stage 23: artifact budget — claimed vs. real gains under each model correction

Status: **pre-registered** (2026-10-23). No stage-23 code exists. This turns
the project's sequence of corrections into one like-for-like table. Each entry
compares what an optimization **claimed** with what holds **after one
correction**, on the **same device, mesh and design**. Definitions are fixed
here, before anything is computed.

## 1. Entries

| ID | Correction | "Claimed" | "Real" | Same device / mesh / design |
|---|---|---|---|---|
| AB1 | objective blind spot (planar cut → dense true minimum) | stage-13 planar-cut BV gain, initial → final, under the stage-12 off-state | dense true-minimum BV change (≥ 40 cuts) of the same two designs, same off-state model, mesh h. **Computed here** | yes |
| AB2 | frozen surface classification (stage 17) | stage-17 BV gain at $R_{on,ch}$(frozen) = R0 | **two readings, both reported, neither preferred:** (a) the 17b gain (classification kept valid), and (b) the stage-17 final under dynamic reclassification (R_on 4.1 % *lower*, so feasible). The sign of the artifact is model-dependent | yes (h/2 on-state, h off-state) |
| AB3 | missing punch-through physics | the avalanche BV (dense true minimum, h) of the 13, 13c, 13d and 17b final designs and their initial design (0.8 µm body) | the blocking voltage $\min(BV, V_{pt})$, with $V_{pt}$ from the stage-19 1D tied-body DD onset ($J_n$ = 1 A/cm²) at the stage-20 worst barrier column of **each** design, and $V_{pt,es}$ (stage-20 barrier criterion) reported alongside. **Computed here** | yes |
| AB4 | bulk → channel mobility | $R_{on}$ (bulk Hall) of the 13d initial and final designs | $R_{on,ch}$ (stage 16) of the same designs, h/2. Already recorded | yes |
| AB5 | exact vs simplified adjoint | stage-22 D-gradient optimum | stage-22 E-gradient optimum: the BV difference and the iteration ratio. Already recorded | yes |

**Metric.** For each entry: claimed, real, and the artifact fraction $1-\text{real}/\text{claimed}$, applied to the claimed quantity. The quantity is a gain for AB1, AB2 and AB5 and an absolute value for AB3 and AB4. Where "real" has no unique definition (AB2), the table says so.

## 2. Criteria

| ID | Criterion | Gating |
|----|-----------|--------|
| T1 | AB1: the recomputed claimed planar-cut gain reproduces the stage-13 recorded value to ≤ 1e-6 relative (same code path) | yes |
| T2 | AB3: the 1D DD reaches the onset for every design (stage-19 numerics), or the reached voltage is reported and that design's entry is marked "lower bound only"; $V_{pt,es}$ changes by ≤ 5 % from h to h/2 | yes |
| T3 | every solve converged; no "On entry to" BLAS lines | yes |
| R1 | the artifact-budget table AB1–AB5, with sources (result file, key) for every number | reported |

Stage 23 passes if T1–T3 hold.

**Not claimed:**
- that any one artifact fraction generalizes beyond this device and model family;
- a ranking of the corrections' "importance" across devices.
