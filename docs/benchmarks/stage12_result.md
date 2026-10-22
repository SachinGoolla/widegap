# Stage 12 result: gate PASSED after one disclosed re-registration (2026-10-21)

These are the first widegap simulations that use **sourced 4H-SiC
parameters** (87-entry registry). The single exception is the
owner-approved convention $\varepsilon_{SiO_2}=3.9$. Implemented by a Sonnet
subagent from the pre-registered spec, then reviewed.

**Runs:**
- **Run 1** (`results/stage12_run1/`, commit `b8001d0`): passed S1 and
  S3–S6, **failed S2′** in one case, Niwa at $N_D=10^{15}$ (order 2.90,
  superconvergence; finest mesh within 4.6e-7 of Richardson).
- **Run 2** (`results/stage12b/`, commit `48fe4b4`, clean tree), under
  S2″ with no order upper bound: **all S1–S6 pass**.
- Amendments are disclosed in spec §4–5.

| ID | Result | Value |
|---|---|---|
| S1 traceability | pass | 63 registry rows tabulated; non-registry design numbers listed separately |
| S2″ BV mesh convergence (10 cases) | pass | orders 1.24–2.90 (or converged to the floor); finest mesh ≤ 3.5e-6 from Richardson |
| S3 electron/hole ionization-integral identity | pass | ≤ 1.7e-10 |
| S4′ KCL, $|I_D|\ge10^{-5}$ | pass | ≤ 1e-8 relative |
| S5 on-current, mesh $h$ vs $h/2$ | pass | 2.36% (gate 3%; margin thin) |
| S6 Jacobian, all models | pass | 7.2e-9 |

**12a: avalanche BV of 4H-SiC p⁺n diodes** (non-punch-through, ∥c, 300 K):

| $N_D$ (cm⁻³) | BV, Niwa 2015 | BV, Hatakeyama 2004 | Ratio |
|---|---|---|---|
| 1e15 | 11.2 kV | 14.2 kV | 0.79 |
| 3e15 | 4.62 kV | 5.58 kV | 0.83 |
| 1e16 | 1.79 kV | 2.03 kV | 0.88 |
| 3e16 | 773 V | 822 V | 0.94 |
| 1e17 | 322 V | 310 V | 1.04 |

Log-log slopes are −0.77 (Niwa) and −0.83 (Hatakeyama), against the
textbook ~$N_D^{-3/4}$; no fitting was done. At $10^{17}$ the Niwa peak
field (3.4 MV/cm) exceeds its stated 1–3 MV/cm validity, which is flagged.
The two published coefficient sets differ by up to 21% in BV. That spread
is a real modelling uncertainty for SiC power-device design.

**12b: vertical 4H-SiC DMOSFET half-cell.** The device: 3 µm drift at
5e16, JFET 1e17, Al body 3e17, 50 nm gate oxide, n⁺ poly gate.
- **Threshold.** $V_T$ = 7.9 V with the sourced flat band ($V_{FB}$ =
  −2.61 V from Afanas'ev offsets) and the Dhar $D_{it}$.
- **$D_{it}$ sensitivity.** ×0 gives $V_T$ 7.6 V; ×10 gives 10.6 V.
- **Fixed-charge sensitivity.** $Q_f=\pm1.5\times10^{12}$ cm⁻² (the
  Afanas'ev bound) moves $V_T$ to 4.3 or 11.3 V. The sign is unknown from
  the source.
- **On-state.** $R_{on,sp}\approx2.9\times10^{-4}$ Ω·cm², using bulk
  Hall mobility with no interface-scattering model, so it is optimistic.
  Velocity saturation is visible in $I_D$–$V_{DS}$.
- **Breakdown.** Avalanche BV is **512 V** on the planar cut under the
  body (Niwa; 522 V Hatakeyama), mesh-stable to 0.03%. The corner cut
  (375–417 V, mesh-dependent) is a geometric singularity and is reported
  only. BV is insensitive to $D_{it}$ and $Q_f$ (< 0.06 V).

**Known limitations and flags** (from `R3_flags` in the results):
- Hall, not drift, mobility; no interface-scattering (channel) mobility
  model.
- The Khan & Cooper $v_{sat}$ (⊥c, 23 °C) is applied to all edges.
- $D_{it}$ is taken uniform from an order-of-magnitude text value.
- The owner-approved $\varepsilon_{ox}$.
- Avalanche BV uses a 1D cut, not the full 2D ionization integral.
- $V_{DS}$ is limited to 15 V by double-precision underflow of
  $\delta e^{\phi_p-\psi}$ in hole rows.
- The Niwa field-range excursion at $10^{17}$.

**Not claimed:** agreement with measured devices. Validation against
experimental 4H-SiC MOSFET data is out of scope here, and these numbers are
model outputs from traceable parameters.

## Post-hoc caveat (2026-10-22, found in stage 13b)

The zero-current off-state electrostatics introduced in stage 12 sets
$\phi_n=V_D$ on nodes with $N_D>N_A$ and $\phi_n=0$ elsewhere. Near the
**lateral** body/JFET junction this creates an artificial electron
reservoir on the p side (about 2e18 cm⁻³ at 400 V), and the depletion
voltage then drops within a single cell. As a result, **peak fields and BVs
near the JFET/body corner are mesh-dependent artifacts**: in stage-13b
diagnostics the corner peak field went 4.0 → 4.5 → 4.8 MV/cm on
successive refinements, at observed order 0.67.

- ~~The planar-cut BV (512 V) appears unaffected~~ **Corrected by stage 12c:** the planar BV was also biased, about 8% low. The corrected model gives 556.46 V, within 3.4% of an independent 1D bound of 576.16 V. See `stage12c_13c_result.md`.
- **Corner-column** BV values (stage 12: 375–417 V; stage 13: 417 → 338 V)
  must not be used quantitatively.
- Stage 13's qualitative conclusion, that the planar-cut objective was
  exploited, stands.

The model fix was implemented in stage 12c (harmonic interpolation of
$\phi_n$); see `stage12c_13c_result.md`.
