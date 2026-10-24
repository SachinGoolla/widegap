# Stage 26: comparison with a fabricated 4H-SiC DMOSFET (Saha & Cooper 2007)

Status: **pre-registered** (2026-10-23). No stage-26 code exists. This is a
**consistency check against one fabricated device with several unsourced
inputs**, not a validation. Each unsourced input gets a range fixed here,
taken from sources or stated plainly as a guess. Nothing is tuned to the
measurement.

**Source.** A. Saha and J. A. Cooper, IEEE TED 54, 2786 (2007), doi:10.1109/TED.2007.904577 ("S07").
- Context: S06 = Saha & Cooper, MSF 527–529, 1269 (2006); M04 = Matin, Saha & Cooper, IEEE TED 51, 1721 (2004).
- Every value below was extracted by one agent and **independently verified on the rendered pages by a second agent** (2026-10-23).

## 0. Verified facts used

- **Structure (S07):**
  - drift epi 6 µm at 1e16 cm⁻³; a 1 µm epi at 1e17 forms both the CSL and the JFET (p. 2789);
  - Fig. 1 depths: P⁺ and N⁺ source 0.25 µm, P-base 0.35 µm more (base total 0.6 µm), CSL 0.3 µm (p. 2787);
  - fabricated channel length about 0.3 µm (pp. 2789–2790);
  - fabricated 1 µm-JFET cell: gate length 3 µm, gate–source separation 0.5 µm, pitch 10 µm (Fig. 8);
  - source contact 3 µm (Figs. 7–8).
- **Oxide (S07):** 50 nm pyrogenic at 1150 °C, NO anneal at 1175 °C for 2 h, C–V thickness 52 nm (p. 2789).
- **Measured (S07):**
  - R_on,sp = **6.95 mΩ·cm²** (1 µm JFET, type-a contacts, area 0.32e-4 cm²). The bias is not stated; V_G = 20 V is only the top-curve label of Fig. 9. Temperature and substrate subtraction are not stated.
  - The 1.5 µm JFET gives 7.5 and the 2 µm JFET 7.9 (type a).
  - BV = 1050 V, avalanche-limited, with a trench termination.
  - TLM contact resistivity 1e-3 Ω·cm².
  - V_T is **not given**.
- **Simulation (S07, MEDICI):**
  - μ = 20 cm²/Vs, field-independent; L_ch = 0.5 µm; ρc = 5e-4;
  - Table I: 3.89 (ρc 5e-4, μ 20), 5.38 (1e-3, 20), 4.53 (5e-4, Fig. 12 fit), 6.02 (1e-3, fit);
  - the 3.89 row is the 8 µm-pitch optimized design (LS = 2 µm). The pitch of the other rows is **inferred** to be the same;
  - S06: simulated VB = 1230 V.
- **Mobility (S07 Fig. 12):**
  - lateral MOSFETs in implanted p-wells, constant-current technique;
  - the x-axis is the **"gate electric field" (MV/cm)**, not $E_{eff}$;
  - fitted curve (approximate, read from the figure): 11, 13.6, 13.4, 11.5, 9 and 7 cm²/Vs at 0.5, 1, 1.5, 2, 3 and 4 MV/cm.
- **Body (M04, same group, earlier process):** Al retrograde with surface concentration 5e17. The M04 simulation used 2.5e17 in the channel and 1e18 near the bottom. No implant schedule is given in any of these papers.

## 1. Model extensions (built and verified first, as stage 26a)

1. **CSL layer.** A separate n layer of thickness $t_{csl}$ under the body, with its own doping, between the body bottom and the drift. It defaults to absent, which keeps every earlier result bit-identical.
2. **Retrograde body.** $N_A(y)$ is log-linear from $N_{A,s}$ at the surface to $N_{A,b}$ at the body bottom, times the existing rounding factors. The acceptor ionization levels and $V_{FB}$ use $N_{A,s}$; this is a documented approximation.
3. **Gate.** $V_{GS}$ = 20 V. The oxide is 52 nm. The gate can extend over the source by $L_{ov}$.
4. **Series terms**, added to R_on,sp after the solve:
   - contact: $\Delta R_c=\rho_c\,x_{max}/L_c$, where $x_{max}$ is the half-pitch and $L_c$ is the modelled source-contact length;
   - substrate: $R_{sub}=\rho_{sub}t_{sub}$.

**Gates for 26a:**

| ID | Criterion | Gating |
|----|-----------|--------|
| E1 | with every extension at its default, the stage-21 initial on-state and off-state values are reproduced to ≤ 1e-12 relative | yes |
| E2 | the CSL and retrograde profiles reproduce their analytic doping at the nodes to ≤ 1e-12 | yes |
| E3 | the series terms match their formulas exactly; the R_on Jacobian/adjoint paths are not required | yes |

## 2. Comparisons (stage 26b)

**Mesh:** h/2, with h/4 reported. $d_{ch}$ = 20 nm.

**(a) Cross-simulator (fewest unknowns).**
- Inputs: S07 simulation conditions; constant μ = 20 cm²/Vs on the channel edges (replacing Noguchi); L_ch = 0.5 µm; Fig. 1 depths; N_J = N_CSL = 1e17; W_J = 1 µm (**full** width, so a half-cell $x_{jfet}$ of 0.5 µm); pitch 8 µm with LS = 2 µm; t_ox 50 nm; V_G = 20 V; no substrate; ρc of 5e-4 and 1e-3.
- Unknown (bracketed): the S07 simulation's body doping, as a uniform $N_A\in\{2.5\times10^{17},5\times10^{17},10^{18}\}$.
- Compare with Table I: 3.89 and 5.38.
- **Criterion A** (reported, with a prediction): Table I lies within our bracket, and the central value ($N_A$ = 5e17) is within ±20 %.

**(b) Mobility model vs measurement.**
- Convert the S07 "gate field" to the SiC surface field, $E_s=(\varepsilon_{ox}/\varepsilon_{SiC})E_{gate}$ ($\varepsilon_{ox}$ = 3.9 owner-approved, $\varepsilon_{SiC}$ = 10.2 registry). Then $N_S=\varepsilon_{SiC}E_s/q-N_{dep}$ and $E_{eff}=(q/\varepsilon_{SiC})(N_{dep}+\eta N_S)$, with η = 1/3 and $N_{dep}$ from the p-well $N_A\in[2.5\times10^{17},5\times10^{17}]$ (M04).
- Evaluate the Noguchi $\mu_{ch}$ and compare with the Fig. 12 fit at gate fields of 1, 2, 3 and 4 MV/cm. Points with $N_S$ ≤ 0 are excluded and reported.
- **Prediction P26b:** within a factor of 1.5 at every included point, for at least one $N_A$ in the range.
- Caveats: Hall vs constant-current mobility, implanted vs epitaxial p-well, and $N_A$ above the Noguchi fit range.

**(c) Measured R_on,sp.**
- **Central inputs (sourced where possible):**
  - fabricated geometry: pitch 10 µm, W_J 1 µm (full), L_ch 0.3 µm, gate length 3 µm, source contact 3 µm;
  - Fig. 1 depths, N_J = N_CSL = 1e17, drift 6 µm at 1e16;
  - retrograde body 5e17 → 1e18 (M04), t_ox 52 nm, V_G = 20 V, ρc = 1e-3 (TLM);
  - two mobility arms: **Noguchi** (main) and the **S07 Fig. 12 fit** (sensitivity arm, using the (b) conversion with the central $N_A$).
- **Ranges (one-at-a-time and corner bracket):**
  - $N_{A,s}\in[2.5\times10^{17},5\times10^{17}]$ (M04 sim vs fab);
  - $N_{A,b}\in[5\times10^{17},10^{18}]$;
  - W_J effective ∈ [0.8, 1.0] µm. S07 attributes part of its own miss to p-well straggle; 0.8 is a **guess**.
  - L_ch ∈ [0.3, 0.4] µm (S07);
  - substrate $R_{sub}\in[0, 0.7]$ mΩ·cm². The upper value is an **unsourced guess**: 350 µm at 0.02 Ω·cm. 0 means the substrate was subtracted.
  - ρc ∈ [5e-4, 1e-3] (S07 sim vs TLM).
- **Prediction P26c:** 6.95 mΩ·cm² lies within the Noguchi-arm bracket, and the Noguchi central value is within ±20 % of 6.95. This is anchored to S07's own 13 % miss (6.02 vs 6.95). The bracket width is reported. If it spans about 2× or more, the comparison is labelled weak.
- **Also reported:**
  - the predicted/measured ratio for both arms, whatever it is;
  - the channel share of R_on (vs S07's qualitative "almost equal" components);
  - the punch-through check (stage-19/20 method) of the central design;
  - the cell avalanche BV vs S06's simulated 1230 V (context; measured 1050 V includes termination and is not gated);
  - the R_on vs W_J trend (1, 1.5, 2 µm) against the measured 6.95, 7.5 and 7.9.

**A disclosed one-at-a-time sensitivity scan** may run before the official run, to rank which unknowns dominate. It may not change any range above.

## 3. Criteria for 26b (methodology)

| ID | Criterion | Gating |
|----|-----------|--------|
| C1 | every solve converged; no "On entry to" BLAS lines | yes |
| C2 | central-design R_on (both arms) changes by ≤ 3 % from h/2 to h/4 | yes |
| C3 | Picard convergence as in stage 16 | yes |
| R1 | (a), (b) and (c) tables and the P26b and P26c outcomes | reported |

Stage 26 passes if E1–E3 and C1–C3 hold. The predictions are reported whatever the outcome.

**Not claimed:** validation of the model in general; agreement with any device but this one; that the unsourced substrate value is correct.

## Note (2026-10-23, after a dev pre-run, before the official run; no criterion, range or prediction changed)

- **Channel-length mapping.**
  - The literal mapping ch_len = L_ch (50 % points, r = 0.1 µm) is degenerate: the 1e20 n⁺ source tail removes all net-p surface at $N_{A,s}$ = 2.5e17.
  - In a first dev pre-run, 5 of 16 bracket corners had no channel columns and failed.
  - The runner therefore implements $L_{ch}$ as the **electrical channel length** (the extent of the net-p gated surface), calibrating ch_len per configuration.
  - The literal mapping is kept as the reported scan row `nominal_l_ch`.
  - This choice was made after seeing the pre-run failure, and is disclosed as such.
- **Other choices,** made where the spec is underdetermined:
  - $t_{csl}$ = 0.3 µm (Fig. 1);
  - central $R_{sub}$ = 0.35 mΩ·cm², the range midpoint, since the spec gave no central value;
  - $V_{DS}$ = 0.1 V;
  - the W_J trend uses pitch 9 + W_J;
  - W_J,eff changes $x_{jfet}$ only.
- **Dev outcome seen:** E1–E3 and C1–C3 passed. Criterion A, P26b and P26c failed. The official run decides.
