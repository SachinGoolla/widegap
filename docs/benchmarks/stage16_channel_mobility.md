# Stage 16: sourced inversion-layer (channel) mobility for the 4H-SiC DMOSFET

Status: **pre-registered** (2026-10-23). No stage-16 code exists. This closes
caveat 15e: R_on is optimistic because the channel uses bulk Hall mobility.
Stage 16 is a **forward re-evaluation only**. Re-optimizing with this model
needs an exact adjoint through a nonlocal mobility and is left to a later
stage.

## 1. Model (Noguchi et al.)

**Sources (full text read):**
- Noguchi et al., IEDM 2017, doi:10.1109/IEDM.2017.8268358 (N17)
- Noguchi et al., JJAP accepted manuscript 2020, doi:10.35848/1347-4065/ab8b3c (N20)

Both are Si-face 4H-SiC with NO-nitrided gate oxide, measured by Hall effect at room temperature. The total channel mobility follows Matthiessen's rule (N20 Eq. 4):

$$\mu_{ch}^{-1}=\mu_{ph}^{-1}+\mu_{SR}^{-1}+\mu_{C}^{-1}$$

| Term | Formula | Source |
|---|---|---|
| phonon | $\mu_{ph}=66.5\,E_{eff}^{-0.39}$ | N17 Eq. (1); N20 Eq. (1) |
| roughness | $\mu_{SR}=76.7\,E_{eff}^{-2}$ | N17 Eq. (2) |
| Coulomb | $\mu_C=\mu_{C,0}(N_S/N_{S,0})^{\gamma}$, $\mu_{C,0}=A(N_{dep}/N_{dep,0})^{B}$, $\gamma=C\log_{10}(N_{dep}/N_{dep,0})+D$ | N20 Eqs. (6)–(8) |
| effective field | $E_{eff}=(q/\varepsilon_{SiC})(N_{dep}+\eta N_S)$, $\eta=1/3$ | N20 Eq. (2) |
| depletion charge | $N_{dep}=\sqrt{(2\varepsilon_{SiC}N_A/q)\,[(2kT/q)\ln(N_A/n_i)-V_B]}$, $V_B=0$ | N20 Eq. (3) |

$\mu$ is in cm²/Vs, $E_{eff}$ in MV/cm, and $N$ in cm⁻². The Coulomb constants are $A=59$ cm²/Vs, $B=-0.939$, $C=0.428$, $D=0.334$, $N_{S,0}=1.5\times10^{12}$ and $N_{dep,0}=10^{12}$ cm⁻² (N20). The N17 Fig. 13 Coulomb set uses a different normalization and is recorded as superseded, not used.

**Pre-registration checks (done before writing this document; reported again in the result):**
- **Units.** $E_{eff}$ in MV/cm comes from the N17 Fig. 6 axis and caption: data on the line $66.5E_{eff}^{-0.39}$, about 163 cm²/Vs at 0.1 MV/cm, which matches the plotted line. In V/cm the formula would give about 0.75 cm²/Vs.
- **Log base in N20 Eq. (8).** "log" is log10. Using registry ε (10.2) and n_i (1.0e-8 cm⁻³) at $V_B=0$:
  - log10 gives γ = 0.229, 0.333, 0.447, 0.578 for N_A = 1e16, 3e16, 1e17, 4e17;
  - N20 measured γ = 0.20, 0.34, 0.46, 0.54;
  - ln gives 0.092, 0.331, 0.593, 0.895.
- **Mesh.** At level 0 the top SiC cells are 3.4 and 5.7 nm, with 3 / 4 / 5 nodes within 10 / 20 / 40 nm of the interface. Levels 1 and 2 have 5 / 7 / 9 and 9 / 13 / 17.

## 2. Discretization (our choices, recorded now)

- **Sheet evaluation per channel column.**
  - For each mesh column $x_i$ whose surface node is p-type ($N_A>N_D$), $N_S(x_i)=\int_{-d_{ch}}^{0}n\,dy$ (mobile electrons only, trapezoidal on nodes).
  - $N_{dep}(x_i)$ uses that column's surface $N_A$; the p-body is laterally graded.
  - $\mu_{ch}(x_i)$ follows from §1.
- **Where it applies.** $\mu_{ch}$ replaces the bulk mobility on **lateral electron edges** with both nodes within $d_{ch}=20$ nm of the interface, in p-type surface columns. Elsewhere the bulk mobility is unchanged: vertical edges, deeper edges, and the n-type JFET or accumulation surface, which Noguchi does not cover.
  - The Ishikawa bulk mobility is **not** added inside the channel. Noguchi's terms are fitted to the total measured inversion-layer mobility, so adding it would double-count lattice scattering.
  - Velocity saturation (Caughey–Thomas, stage 9a) still multiplies the low-field mobility on those edges.
- **Validity clamp.** $N_S\leftarrow\max(N_S,10^{11})$ and $E_{eff}\leftarrow\max(E_{eff},10^{-2})$ MV/cm. These are the lower ends of the plotted fit ranges (N17 Figs. 6 and 8); they are model choices, not paper values. Clamped columns are counted, as are columns with $N_S>10^{13}$ (an extrapolation).
- **Nonlocality.** $\mu_{ch}$ depends on $N_S$, an integral over the column, so it is lagged in an outer Picard loop:
  - solve the DD with fixed $\mu$;
  - update $\mu_{ch}$;
  - repeat until $\max|\Delta\mu_{ch}/\mu_{ch}|<10^{-10}$, at most 60 iterations.
- **Hall factor** taken as 1, as Noguchi does.

## 3. Acceptance criteria

| ID | Criterion | Gating |
|----|-----------|--------|
| H1 | the registry entries for §1 pass `check_entry`, including independent re-verification against the **rendered page images** (signs and exponents) | yes |
| H2 | code-path limit: with $\mu_{ch}$ replaced by the bulk mobility on the channel edges, stage-12 R_on,sp of the 12b device is reproduced to ≤ 1e-12 relative | yes |
| H3 | discretization check on a long-channel lateral n-MOSFET (uniform $N_A\in\{10^{16},10^{17},3\times10^{17}\}$, L = 10 µm, $V_{DS}$ = 0.05 V, $V_{GS}$ giving $N_S\approx3\times10^{12}$ at mid-channel): $\mu_{eff}=I_DL/(qWN_SV_{DS})$ matches $\mu_{ch}$ at the simulated mid-channel $(E_{eff},N_S)$ to ≤ 2 % | yes |
| H4 | Picard convergence (§2) for every evaluation, with every inner Newton converged | yes |
| H5 | mesh: R_on,sp of the 12b device with the channel model differs by ≤ 3 % between $h$ and $h/2$ | yes |
| H6 | $d_{ch}$ insensitivity: R_on,sp changes by ≤ 1 % between $d_{ch}$ = 20 and 40 nm (12b, $h/2$) | yes |
| R1 | R_on,sp, V_T and the channel share of R_on for 12b and the 13d initial and final designs: bulk-Hall vs channel model; also at $h/4$ | reported |
| R2 | clamp counts; columns above 1e13; $d_{ch}$ = 10 nm result; $\mu_{ch}$ along the channel | reported |
| R3 | a γ(N_A) and $\mu_{ph}$ reproduction table (§1 checks) | reported |

Stage 16 passes if H1–H6 hold.

**Caveats recorded now:**
- Hall factor 1.
- NO-annealed devices only.
- The Coulomb term does **not** depend on $N_{it}$ or $Q_f$. The registered $D_{it}$ (Saks/Dhar) and $Q_f$ still set V_T and $N_S$, but they do not enter the mobility.
- 300 K only.
- No mobility model on the n-type accumulation surface.
- The 13d optimum was found with the bulk model, so it is not re-optimized here.

**Not claimed:** agreement with measured R_on of any commercial device.

## 4. Oxide-field reference (doc update, not a criterion)

- Shi et al., Mater. Sci. Semicond. Process. 174, 108194 (2024), doi:10.1016/j.mssp.2024.108194, p. 6, state an "industry standard of 20 years at 4 MV/cm and 150 °C" for gate-oxide lifetime under **positive gate bias** (TDDB). The paper gives no citation for it.
- The 13d E_ox (5.25 MV/cm at 400 V) is an **off-state, drain-induced** field, which is a different stress polarity.
- The README and result wording become: "5.25 MV/cm vs the 4 MV/cm gate-bias TDDB convention stated by Shi et al. (2024), 1.3×, different stress polarity". The unsourced 3 MV/cm guideline is dropped.
- The paper's 11 MV/cm is an assumption it makes to estimate oxide thickness, and is not cited.
