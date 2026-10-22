# Stage 8.1: MOS structures in 2D (oxide, gate, and a MOSFET half-cell)

Status: **pre-registered** (2026-10-21). No stage-8.1 code exists at the
time of this commit.

## 1. Model additions to stage 8

- **Regions.** Each primal cell is semiconductor or oxide. The Poisson
  equation holds everywhere, with cell-wise $\varepsilon$. Carrier
  equations hold only on nodes that touch a semiconductor cell. Carrier
  edge fluxes use only semiconductor cells' dual-face weights, and carrier
  volumes use only the semiconductor part of each control volume. The
  space charge acts only in the semiconductor part. Nodes entirely inside
  oxide carry identity rows for $\phi_n,\phi_p$. The semiconductor/oxide
  interface therefore has zero carrier flux automatically.
- **Gate contact.** Dirichlet $\psi=V_G-V_{FB}$ on the gate boundary
  (scaled units).
- **Study parameters, not sourced material values:**
  $\varepsilon_{ox}/\varepsilon_s=0.38$, and $V_{FB}$ given per case.
  Physical oxide, work-function and interface-charge values require
  sourced registry entries in a later stage. No interface traps or fixed
  charge in this stage.
- **Carriers.** Quasi-Fermi DD (stage 7.1/8): complete ionization,
  constant mobility, $R\equiv0$ unless stated.

## 2. Exact reference for the MOS capacitor (1D)

A p-type substrate (net acceptor $a$, scaled) under an oxide of thickness
$t_{ox}$, in equilibrium. Integrating Poisson's equation once from the
neutral bulk gives the semiconductor charge as a closed function of the
surface potential $\psi_s$ (measured from the bulk):

$$
Q_s(\psi_s)=-\operatorname{sign}(\psi_s)\sqrt{2\big[a(e^{-\psi_s}+\psi_s-1)\big]+2\big[\delta e^{\psi_b}(e^{\psi_s}-\psi_s-1)\big]}
$$

with Boltzmann carriers, and $\psi_b$ the bulk potential. Gauss's law
across a charge-free oxide gives
$V_G-V_{FB}-\psi_{b}=\psi_s-\tfrac{\varepsilon_s}{\varepsilon_{ox}}t_{ox}Q_s(\psi_s)$
(in scaled units, $\varepsilon_s=1$). The implementation must derive the
exact expressions in its docstring, consistent with stage-2/7 signs; the
form above is indicative. The resulting $\psi_s(V_G)$ is solved by a 1D
root-find independent of the FV code.

## 3. Cases ($\delta=10^{-24}$ realistic-SiC scale unless stated)

- **C1 (1D MOS capacitor):** p-substrate $a=10$, $t_{ox}=0.5$ (in units
  of $L_D$, a study value), semiconductor depth 30, $V_{FB}=0$. Equilibrium
  ($\phi=0$). $V_G$ from accumulation through depletion to inversion
  (range chosen so that $\psi_s$ spans $[-5,\,2\psi_F+10]$).
- **C2:** C1 extruded in 2D (zero flux on the sides).
- **C3 (MOSFET half-cell):** p-body $a=10$; n⁺ source and drain wells
  ($b=100$) at both ends; gate over the channel through oxide $t_{ox}$;
  source, drain and body contacts. Symmetric geometry with channel length
  $L$ and source/drain spacing as set by the implementation, all recorded.
  Bias $V_{DS}$ small and $V_{GS}$ swept through threshold.

## 4. Acceptance criteria

| ID | Case | Criterion | Gating |
|----|------|-----------|--------|
| G1 | C1 | numerical $\psi_s(V_G)$ (node at the interface), Richardson-extrapolated over two meshes, matches the exact reference to ≤ 2e-3 relative (or ≤ 1e-4 absolute when $|\psi_s|<0.05$) at every $V_G$ | yes |
| G2 | C1 | the gate-side oxide field times $\varepsilon_{ox}$ (Gauss) equals $-Q_s$ from the semiconductor charge integral of the numerical solution, to ≤ 1e-10 relative (discrete conservation) | yes |
| G3 | C2 | every 2D column equals the 1D C1 solution to ≤ 1e-10 | yes |
| G4 | C3 | Newton converges at every bias (per-row tolerance) | yes |
| G5 | C3, $|I|\ge10^{-6}$ | source current + drain current + body current = 0 to ≤ 1e-8 relative to $|I_D|$ | yes |
| G6 | C3, symmetric device | $I_D(V_{DS})=-I_D(-V_{DS})$ with source and drain exchanged, to ≤ 1e-8 relative, where $|I|\ge10^{-6}$ | yes |
| G7 | C3 | $I_D$ at one on-state bias on three meshes ($h,h/2,h/4$): observed order in [0.8, 2.5], and the finest mesh within 2% of the Richardson extrapolation | yes |
| G8 | C3 | analytic Jacobian (with region-masked assembly) vs FD, per row, ≤ 1e-6 | yes |
| R1 | C3 | transfer curve $I_D(V_{GS})$, threshold by linear extrapolation, and gradual-channel comparison $I_D\approx\mu\frac{W}{L}C_{ox}(V_{GS}-V_T)V_{DS}$ | reported |
| R2 | C3 | figure: potential and electron density maps in the on state | reported |

## 5. Gate and claims

Stage 8.1 passes if G1–G8 hold. A pass verifies the MOS **machinery**
(regions, gate, interface, 2D MOSFET currents). It is **not** a 4H-SiC
MOSFET prediction. Sourced oxide, interface-trap and channel-mobility
physics are needed for that, and they are later stages.

## 6. Amendments

### 2026-10-21: G2 round-off floor and G6 wording (before the official run, after a development run)

**Disclosure.** Both points were raised by the implementing agent after a
development run, which is not an official result.

- **G2.** A relative error on $Q_s$ is ill-posed at flat band, where
  $Q_s=0$. This is the same flaw as stage 6's gradient check at the
  optimum. Restated: $|\text{err}|\le10^{-10}|Q_s|+10^{2}\epsilon_{mach}S$,
  where $S$ is the magnitude scale of the charge and Gauss terms. That is
  the explicit round-off floor. At every non-flat-band point the
  $10^{-10}$ relative bound must hold by itself.
- **G6, precise wording.** With fixed $V_G$, the case $V_D=-V_{DS}$ is not
  the mirror of $V_D=+V_{DS}$, because $V_{GD}$ differs. The exact symmetry
  is the **mirror test**. In configuration A, source = 0 and drain = $v$;
  in configuration B, source = $v$ and drain = 0, with the same $V_G$ and
  the body at 0. The test requires $I_{\text{right}}(A)=-I_{\text{right}}(B)$
  to ≤ 1e-8 relative. "With source and drain exchanged" in §4 meant
  exactly this.

Thresholds are otherwise unchanged.
