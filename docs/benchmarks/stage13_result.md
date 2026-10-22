# Stage 13 result: gates PASSED, but the optimized design is NOT a better device (2026-10-22)

## Read this first: the objective was gamed

The pre-registered objective measured avalanche BV along a **fixed planar
cut** under the body contact. The cut was chosen to avoid the mesh-dependent
field singularity at the JFET/body corner, and the optimizer exploited that
blind spot:

| Quantity | Initial | Optimized |
|---|---|---|
| BV on the planar cut (the objective) | 512 V | **633 V** (+121 V, mesh-confirmed) |
| Peak field at the JFET/body corner | 4.54 MV/cm | **6.32 MV/cm** |
| BV on the corner column | 417 V | **338 V** |

A real device would break down at the corner, at a **lower** voltage than
before. The +121 V is real for the objective as specified and physically
misleading for the device. This is a **specification flaw** by the lead,
not an implementation error, and stage 13b addresses it. Nothing below
should be read as "a better SiC MOSFET".

## What *is* verified

Raw data: `results/stage13/` (commit `28c590d`, clean tree). Implemented by
a Sonnet subagent from the pre-registered spec, then reviewed. The official
numbers are identical to the development run.

| ID | Result | Value |
|---|---|---|
| T1 $dR_{on}/d\theta$ vs Richardson FD | pass | 6.7e-9 relative |
| T2 $dBV/d\theta$ vs Richardson FD | pass | 2.0e-8 relative |
| T3 final-design gradients | pass | BV 2.8e-7 V absolute (tolerance 1.1e-4) |
| T4 feasibility | pass | $R/R_0-1=2\times10^{-12}$ |
| T5 re-verification on $h$ and $h/2$ | pass | +120.9 V / +121.4 V; $R_{on}$ differs 0.7% between meshes |
| T6 nested family | pass (weak) | the initial design is itself the best feasible uniform design |
| T7 SLSQP | pass | exit 0, 15 iterations |

The **full-physics adjoint machinery** is verified for these objectives.
That covers doping entering through the ionized fractions and the
doping-dependent mobility (analytic $\partial F/\partial\theta$, matching
residual FD to about 1e-9), the full-DD $R_{on}$ adjoint, and the
implicit-function avalanche-BV gradient.

**Disclosed deviation (T3).** The spec's tolerance
$10^{-6}\max(1,\lVert\nabla\rVert)$ depends on units and is vacuous for
$R_{on}$ in Ω·cm². It was applied to $R/R_0$ instead. The literal check in
Ω·cm² also passes and is reported.

**Reported (R1–R3).** The Pareto set lies 4.0–6.6× above the 1D unipolar
4H-SiC limit computed from the same registry parameters, with no fitting.
The optimized drift profile is non-monotone (6.6e16 → 1.65e16 → 1e17
cm⁻³), and the JFET sits at its upper bound. All of this is *for the flawed
objective*.

## Lesson

An optimizer will find any region the objective cannot see. Breakdown
objectives must cover **every** location that can break down first, either
as a smooth minimum over candidate paths or as the full 2D ionization
integral, and the corner singularity must be resolved physically, for
example with a rounded junction.

## Post-hoc caveat (2026-10-22, found in stage 13b)

The zero-current off-state electrostatics introduced in stage 12 sets
$\phi_n=V_D$ on nodes with $N_D>N_A$ and $\phi_n=0$ elsewhere. Near the
**lateral** body/JFET junction this creates an artificial electron
reservoir on the p side (about 2e18 cm⁻³ at 400 V), and the depletion
voltage then drops within a single cell. As a result, **peak fields and BVs
near the JFET/body corner are mesh-dependent artifacts**: in stage-13b
diagnostics the corner peak field went 4.0 → 4.5 → 4.8 MV/cm on
successive refinements, at observed order 0.67.

- The **planar-cut** BV (stage 12: 512 V) was mesh-stable and appears
  unaffected.
- **Corner-column** BV values (stage 12: 375–417 V; stage 13: 417 → 338 V)
  must not be used quantitatively.
- Stage 13's qualitative conclusion, that the planar-cut objective was
  exploited, stands.

A model fix is planned (a quasi-Fermi level that varies through the
depletion region, or full drift–diffusion in the off state).
