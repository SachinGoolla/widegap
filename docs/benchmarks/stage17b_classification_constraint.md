# Stage 17b: stage 17 with the surface-classification guard enforced as a constraint

Status: **pre-registered** (2026-10-23), written after the stage-17 dev run and
before the stage-17 official result. This is a disclosed re-registration
prompted by the stage-17 J4 finding. Stage 17 is run and recorded unchanged.

## 1. The finding (stage-17 dev run)

- The channel model classifies each surface column once, at θ = 0, as p-type (inversion, Noguchi) or n-type (bulk). The adjoint is exact for that frozen map.
- One gated column at the rounded JFET/body corner (x ≈ 0.733 µm on h/2) starts p-type with $N_D/N_A$ = 0.93. It becomes n-type at θ_J ≈ 0.074.
- The stage-17 optimum (θ_J = 0.147) and every Pareto design pass that point. Re-classifying the column at the stage-17 final lowers R_on,ch by 4.1 %.
- So the binding R_on constraint was partly evaluated in a regime the frozen model does not represent. J4 (a guard, not a constraint) failed on 113 of 173 evaluations.

## 2. Change

Stage 17b adds one family of SLSQP inequality constraints:

$$c_k(\theta)=\frac{N_{A,k}-N_{D,k}(\theta)}{N_{A,k}}-m\ \ge 0,\qquad m=0.02,$$

- **Which columns.** Every surface column that is p-type at θ = 0 and lies under the gate, on **both** the h/2 (on-state) and h/4 (J6 re-check) meshes. Doping is evaluated at each column's surface node.
- **Gradient.** These are explicit, smooth functions of θ, so the gradients are analytic and checked by FD.
- **Effect.** The design is restricted to the region where the frozen classification is valid. This is a model-validity restriction, not a physical design rule.
- **Everything else** is stage 17 unchanged: the problem, the meshes, the adjoint, the Pareto levels, and the J1, J2, J3, J5, J6 and J7 criteria with the same thresholds.

## 3. Criteria (changes only)

| ID | Criterion | Gating |
|----|-----------|--------|
| J4b | at every SLSQP iterate (callback) and at the final and Pareto designs, no gated surface column changes sign relative to θ = 0, on h/2 or h/4. Clamp margins are as in J4. Line-search trial points that cross are counted and **reported** | yes |
| J8 | classification-constraint gradients vs FD ≤ 1e-8 relative | yes |
| R4 | active classification constraints at the optimum; BV cost of the constraint (17b vs 17 final dense BV); re-classified R_on,ch of the stage-17 final for comparison | reported |

Stage 17b passes if J1, J2, J3, J4b, J5, J6, J7 and J8 hold.

**Not claimed:**
- that the margin m = 0.02 is physical;
- a model of the compensated corner. Dynamic reclassification or a smooth inversion/accumulation blend would need a sourced model, which we do not have.
