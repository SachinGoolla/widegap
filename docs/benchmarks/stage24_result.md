# Stage 24 result: gate NOT passed (X1, X4 pass; X2, X3 FAIL, spec error); predictions P24a to P24c reported (2026-10-24)

Raw data: `results/stage24/` (`summary.json`, `results.json`, `manifest.json`). Figures:
`figures/stage24_theta.png`, `figures/stage24_blocking.png`, `figures/stage24_barrier_iterates.png`.
Pre-registration: `stage24_guard_ablation.md`. Context: `stage23_result.md`, `stage21_result.md`,
`stage20_result.md`, `stage19_result.md`. All numbers below are quoted from the JSON files above
(`summary.json` unless another file is named); anything not recorded there is marked "derived".

## Verdict

`criteria`: **X1 true, X2 false, X3 false, X4 true**; `gate_passed: false`, `dev: false`.
**Stage 24 does not pass.**

| ID | Criterion | Measured | Result |
|---|---|---|---|
| X1 | B_tau adjoint vs Richardson FD <= 1e-4 at theta = 0 and at the Arm-G final; other gradients as in stage 21 S1 | `values/X1/B_tau` `rel_err_norm` 2.64e-8 (theta 0), 7.22e-10 (G final), 3.84e-11 (U final, extra); `S1_theta0` R 8.2e-9, BV 1.6e-8, E_ox 1.5e-12; `S1_finals` U and G `ok` true | pass |
| X2 | each arm feasible at its own final design; Arm G also B_tau >= 10 V_T (1 - 1e-8) | `values/X2/U/ok` true; `values/X2/G/ok` **false**: `B_tau_650_VT` 1.696 against the required 10 | **FAIL** |
| X3 | SLSQP exit 0 for both arms and the restoration; all solves converged; no "On entry to" lines | `values/X3/status` U 0, **G 8** ("Positive directional derivative for linesearch"); restoration `status` 0 but `reached_feasible_B` false; `on_entry_to_lines` 0, `job_errors` [] | **FAIL** |
| X4 | V_pt,es change <= 5 %, dense BV change <= 3 % between h and h/2 at each final | `values/X4` U `dV_pt_rel` 1.329 %, `dBV_rel` 0.103 %; G `dV_pt_rel` 1.341 %, `dBV_rel` 0.099 % | pass |

X2 and X3 fail for the same reason, and the reason is the pre-registered guard level, not the numerics: see
the next section.

## Why X2 and X3 fail: the guard level was unreachable (a SPEC ERROR by the planner)

The pre-registered guard is B_tau(theta; V* = 650 V) >= 10 V_T. On this design family that level cannot be
reached:

- At theta = 0, `values/R1/G_B_tau_at_theta0` = **1.4738 V_T** (`B_tau_650_VT`; hard-min barrier at 650 V
  `hard_barrier_650_VT` 1.1746 V_T).
- The restoration (SLSQP maximizing B_tau subject to the common constraints, `values/R1/G_restoration`)
  exited 0 after 18 iterations at B_tau = **1.6962 V_T** (hard 1.5053 V_T), `reached_feasible_B` false. The
  Arm-G final reaches the same value (`values/X2/G/B_tau_650_VT` 1.6962). The Arm-U final has 1.6520 V_T
  without any guard. The restoration maximum is therefore about **1.70 V_T**, i.e. 17 % of the required 10 V_T.
- Arm G started from the restoration end point and its constraint B_tau >= 10 V_T was infeasible throughout.
  SLSQP exited with status 8 after 21 iterations (`values/X3/status/G`, `nit` G 21). X2 and X3 then fail by
  construction.

**Cause.** V* = 650 V was chosen above every BV this design family reached (the stage-21 optimum was 608 V),
but it also lies above V_pt,es of the initial design (578.24 V h, 586.38 V h/2; `values/R1/blocking/initial`
`V_pt_es_h`, `V_pt_es_h2`). The spec assumed the restoration would reach 10 V_T at 650 V; it does not, so the
barrier at V* is under-depleted for every design in the family reachable under the common constraints. This is
a **spec error by the planner**, not an optimizer or numerics failure. The recorded `min_constraint_value`
(-6.1e-14) for Arm G is round-off level; the X2 failure is the B_tau shortfall (derived reading: `ok` is a
single flag and the B_tau clause is the only quantity recorded far from its bound).

**No re-registration with a different V* is made.** The Arm-U outcome is already known. Choosing a new guard
level after seeing it would tune the guard to the observed result and would no longer be a prospective test.
The failed gate stands as recorded.

Not gated, noted: `nonfinite_calls_total` 112, all `nonfinite_in_failed_newton`, 0 in the adjoint and in
converged Newton solves (`values/X3`); `failed_picard`, `failed_evals_on`, `failed_evals_off` all 0.

## R1: per-arm table

Source: `values/R1/blocking/{initial,U,G}` and `values/R1/theta`. theta is the 5-vector, start from 0 for U.
Blocking voltage = min(dense BV, V_pt,es). R_on,ch is on h/2 (`R_on_ch_h2`, ohm cm2); E_ox at 400 V on h
(units as recorded).

| Quantity | Initial (theta = 0) | Arm U | Arm G |
|---|---|---|---|
| theta | 0 | [-0.15994, -0.03663, 0.85309, 0.03663, 0.008567] | [-0.08719, 0.12214, 0.70434, -0.04615, 0.008567] |
| Dense BV h (V) | 547.126 | 612.701 | 602.910 |
| Dense BV h/2 (V) | 545.974 | 612.068 | 602.311 |
| Dense BV gain over initial, h / h/2 (V, derived) | 0 | 65.576 / 66.094 | 55.785 / 56.336 |
| Field-line BV h / h/2 (V) | 545.174 / 541.573 | 608.566 / 605.863 | 599.759 / 597.233 |
| V_pt,es h / h/2 (V) | 578.238 / 586.383 | 604.085 / 612.114 | 609.048 / 617.216 |
| Blocking h / h/2 (V) | 547.126 / 545.974 | 604.085 / 612.068 | 602.910 / 602.311 |
| Hard-min barrier at 650 V (V_T) | 1.1746 | 1.4448 | 1.5053 |
| B_tau at 650 V (V_T) | 1.4738 | 1.6520 | 1.6962 |
| R_on,ch h/2 (ohm cm2) | 2.45948e-3 | 2.45948e-3 | 2.45948e-3 |
| E_ox | 4.2774 | 4.2246 | 4.1938 |
| Active constraints (`active_constraints`) | n/a | R_on,ch; classification index 12 | R_on,ch; classification index 12 |
| Not active | n/a | E_ox, B_tau, bounds | E_ox, B_tau, bounds |
| SLSQP | n/a | exit 0, 10 iterations, 10 evaluations | exit 8, 21 iterations, 129 evaluations (from the restoration point) |

The initial design is not punch-through limited on either mesh (V_pt,es/dense BV 1.057 on h, 1.074 on h/2,
as in stage 20). R_on,ch is equal to its start value in both finals (the R constraint is active).
`B_tau` is not active in Arm G because the constraint is infeasible, not because it is slack. The 1D DD
onset was not run (no DD job in `results.json` `diagnostics`; not gated).

## R2: predictions (thresholds unchanged from the spec)

| Prediction | h | h/2 | Outcome |
|---|---|---|---|
| P24a: Arm-U V_pt,es < dense BV | 604.085 < 612.701 V: true | 612.114 vs 612.068 V: **false** (V_pt,es is 0.046 V above BV) | **CONFIRMED on h, NOT on h/2** |
| P24b: Arm-U blocking < 1.05x the initial blocking | ratio **1.1041** (604.085 / 547.126 V) | ratio **1.1211** (612.068 / 545.974 V) | **REFUTED** |
| P24c: Arm G has V_pt,es >= BV, dense-BV gain at least 5 V smaller than U, blocking higher than U | parts: G not punch-through limited true; gain gap 9.79 V >= 5 true; blocking G 602.910 > U 604.085 **false** | same parts: true, gap 9.76 V, blocking G 602.311 > U 612.068 false | **REFUTED** (weak evidence, see below) |

(`values/R2/{h,h2}/P24a`, `P24b`, `P24c`.)

- **P24a.** On h the prediction holds with 8.6 V of margin (derived: 612.701 - 604.085). On h/2 the two quantities
  differ by 0.046 V, which is 0.04 V in the rounded figures: the difference is inside mesh resolution (the
  recorded h-to-h/2 changes are 1.3 % in V_pt,es and 0.10 % in BV, `X4`), so the h/2 result is a tie, not
  evidence for either side. The P24a result is mesh-resolved only on h.
- **P24b.** The Arm-U blocking voltage rose by 56.96 V on h (+10.4 %) and 66.09 V on h/2 (+12.1 %, derived), about
  twice the 5 % threshold. V_pt,es itself rose from 578.24 V to 604.08 V (h), so the unguarded optimizer moved
  the punch-through voltage up, not down.
- **P24c.** It is evaluated against an Arm G that is infeasible (B_tau 1.70 V_T against 10 V_T required,
  exit 8) and whose optimization was therefore not the pre-registered guarded problem; it is weak evidence in
  either direction. It fails on a single component (blocking G below blocking U by 1.17 V on h and 9.76 V on
  h/2, derived), and the other two components hold. Because Arm G is not a valid guard arm, the result says
  nothing about whether a working guard would help.

## Key conclusion

The prospective test does **not** support the claim that an unguarded optimizer drifts toward the unmodelled
failure mode on this geometry. In Arm U, the optimizer raised both the avalanche BV (547.1 to 612.7 V, h) and
V_pt,es (578.2 to 604.1 V, h); the blocking voltage rose by 10.4 % (h) and 12.1 % (h/2). The stage-23 post hoc
pattern (the 0.8 um designs 13d and 17b gained avalanche BV while their blocking voltage fell) **does not
replicate** on the 1.2 um body. This is the single prospective case registered; it refutes the general
statement only for this geometry and model family.

The supported, narrower statement: avalanche-BV objectives do not track blocking capability in a fixed
direction. The coupling between them is geometry-dependent (down on the 0.8 um body in the post hoc data, up on
the 1.2 um body in this prospective test), so every failure mode has to be checked for each design rather than
assumed to be improved or degraded by the objective.

## Interpretation (not a measured result)

On the 1.2 um body, the design change that raised the avalanche BV also raised V_pt,es; the mechanism is not
tested here. Arm G ends with a lower BV and a higher V_pt,es than Arm U (609.0 against 604.1 V on h), but it is
the end point of a restoration toward a higher barrier, not of a valid guarded problem. Anything beyond this
requires a corrected guard level and is not claimed.

## Caveats

- One geometry (d_body = 1.2 um), one coefficient set, one starting point per arm (the manifest lists `seeds` [0, 1, 2]; I did not check whether they affect the optimization).
- Arm G is not a valid guarded problem (infeasible constraint); P24c is weak evidence.
- Blocking uses the electrostatic V_pt,es criterion of stage 20, not a 2D or DD current-carrying breakdown; no
  DD onset was computed for the finals.
- The soft-min `B_tau` exceeds the hard-min barrier in every record (e.g. 1.4738 against 1.1746 V_T at theta
  = 0, 1.6962 against 1.5053 V_T at the G final). A smooth minimum of the nested log-sum-exp form is normally
  not above the hard minimum. Explanation (planner, from the pre-registered definition): B_tau is nested, psi_surface minus a soft-min over y of psi. A log-sum-exp soft-min lies BELOW the true minimum, so the inner term raises each column's barrier above its hard value. The outer soft-min over columns lowers it again by at most tau*ln(N_columns). The net B_tau can therefore exceed the hard-min barrier. That makes the guard slightly LENIENT (it over-states the barrier), not strict. It does not affect the verdict, since even the lenient guard was unreachable.
  It does not change the conclusion (both are far below 10 V_T).
- Both mesh-resolution statements (P24a on h/2) rest on the recorded h versus h/2 changes only.

## Manifest

Code commit `453cf8adc23d04a4be501925812093b98e2172e6` (`453cf8a`), `git_dirty: false`,
2026-10-24T01:02:49 UTC, `benchmarks/stage24/run.py --workers 12`, runtime 1103.4 s (`runtime_seconds`).
`widegap 0.6.0`. Python 3.12.3, NumPy 2.5.4, SciPy 1.18.1, torch 2.14.1+cpu (1 thread), AMD Ryzen 7 5700U,
16 logical CPUs, Linux 7.0.0-38-generic. `seeds` [0, 1, 2] (recorded in the manifest); `dev: false`.

## Not claimed

- That an unguarded optimizer drifts toward punch-through in general (refuted prospectively on this geometry).
- That a punch-through guard improves or degrades blocking capability: the guarded arm was infeasible.
- A different guard level (V*) as a repaired experiment: no such experiment is registered or run.
- A 2D punch-through voltage, generality beyond this geometry and model family, or agreement with measured
  devices.
