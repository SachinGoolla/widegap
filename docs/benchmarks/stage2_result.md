# Stage 2 result: gate FAILED (S2), 2026-10-20

Raw data: `results/stage2/` (manifest, results, summary).

| Criterion | Result | Value |
|---|---|---|
| S1 MMS order (E1) | pass | 2.001, 2.000, 2.000, 2.000 |
| S2 quadratic Newton (E1, N=256) | **fail: not evaluable** | only 2 unlimited iterations (ratios 3e-10, 1.1e-3) |
| S2 quadratic Newton (E2, N=1024) | **fail** | final three unlimited ratios 4.9, 4.6, **1.8e3** (> 1e3) |
| S3 bulk neutrality (E2) | pass | 5.0e-7 (≤ 1e-6) |
| S4 total charge (E2) | pass | 0.0 to round-off |
| S5 peak field vs depletion approx (E2) | pass | 0.66 % (≤ 10 %) |
| S6 (E3, asymmetric 100:1, reported) | n/a | neutrality 9.2e-6 at \|x\|≥20; peak field 18 % above the depletion approximation |

## Diagnosis

The solver converges in every case, to $\lVert F\rVert_\infty\le10^{-12}$, in
3 (E1), 12 (E2) and 24 (E3) iterations. Every failing S2 ratio is the
**last** Newton step, which lands on the double-precision floor of the
residual. The flux terms are of order $|\psi|/h\approx 2\times10^2$, so the
floor is about $10^{-13}$. A step into round-off cannot show quadratic
contraction: E3's last ratio is 6.8e10 for the same reason. Before that final
step the ratios are O(1) (E2: 4.9, 4.6; E3: 0.12, 0.11), which is consistent
with quadratic convergence.

E1 converges so quickly that "the final 3 unlimited iterations" do not exist.

**Conclusion:** the evidence points to a flaw in the criterion's wording,
not in the solver. Under the project rules, a failed gate stays failed. The
criterion is not reinterpreted after the results are known, and stage 3 does
not start.

## Proposed path (requires a new pre-registration, not applied here)

Register S2′ in a new dated document *before* a fresh run:

- quadratic contraction is evaluated only on steps whose output residual
  exceeds $10^3\,\epsilon_{mach}\max_i(\text{flux and source term magnitudes})$;
- at least 2 such steps are required.

Then run the stage-2 benchmark again from a clean commit.

## E3 note (reported only)

For the 100:1 junction the depletion approximation underestimates the
numerical peak field by 18 %. The neutrality bound at $|x|\ge20$ is not met
on the lightly doped side, whose depletion region and Debye tail are wider.
Both are expected qualitatively, and neither is gating. The E3 mesh
($h=0.078\,L_D$, against a heavy-side Debye length of $0.1\,L_D$) is
marginal, and a refinement study belongs to S2′'s rerun.
