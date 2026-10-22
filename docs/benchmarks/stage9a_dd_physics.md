# Stage 9a: additional drift–diffusion physics (incomplete ionization, field-dependent mobility, interface charge)

Status: **pre-registered** (2026-10-21). No stage-9a code exists. All model
constants here are **study parameters** chosen to exercise each model.
4H-SiC values enter later only through the provenance registry (stage 7
rules). Lessons applied: relative errors only away from zeros (explicit
round-off floors otherwise), exact references, and resolvability limits.

## 1. Models (added to the stage-8/8.1 2D quasi-Fermi core; 1D is the special case)

1. **Incomplete ionization out of equilibrium.** Donor occupancy follows
   the electron quasi-Fermi level, and acceptor occupancy the hole level:
   $N_D^+=N_D/(1+e^{\psi-\phi_n+\ell_D})$ and
   $N_A^-=N_A/(1+e^{\phi_p-\psi+\ell_A})$, with $\ell_{D,A}$ as in stage 7.
2. **Field-dependent mobility** (Caughey–Thomas). On each edge,
   $\mu=\mu_0/[1+(\mu_0F/v_{sat})^{\beta}]^{1/\beta}$. The driving force
   $F=|\phi_{b}-\phi_{a}|/h_e$ is the quasi-Fermi gradient along the edge,
   separately for electrons and holes. Newton includes the exact
   derivative. At $F=0$ the derivative is taken as the analytic limit,
   which is smooth for $\beta\ge1$.
3. **Interface charge** at semiconductor/oxide interface nodes, as a sheet
   density per unit interface length entering the Poisson row:
   - fixed charge $Q_f$;
   - trap charge $Q_{it}$ for a uniform $D_{it}$ over an energy window
     $[E_1,E_2]$ relative to midgap, acceptor-like above midgap and
     donor-like below. Occupancy is Fermi–Dirac with respect to the local
     electron quasi-Fermi level. Its energy integral has a closed form
     (log terms), which the implementation derives and documents.

## 2. Acceptance criteria

| ID | Criterion | Gating |
|----|-----------|--------|
| A1 | Equilibrium ($\phi=0$) 1D junction with incomplete ionization through the DD core equals the stage-7 `IonizedJunction` solution: $\max|\Delta\psi|\le10^{-10}$ | yes |
| A2 | MMS with incomplete ionization, quasi-Fermi DD, 1D and 2D (smooth manufactured $\psi,\phi_n,\phi_p$, sources derived symbolically): observed orders, last three, in [1.8, 2.2] | yes |
| A3 | Field-dependent mobility, $v_{sat}\to\infty$ ($10^{30}$): currents equal the constant-mobility solution to ≤ 1e-12 relative where $|j|\ge10^{-6}$ | yes |
| A4 | Uniform n-type bar ($N_D$ uniform, $\delta=10^{-3}$, ohmic contacts, length $L$), $V$ up to $20\,v_{sat}L/\mu_0$: $j=qn\mu(F)F$ with $F=V/L$ and $n=N_D$ in the bulk. The interior element current matches the closed form to ≤ 1e-6 relative (contact layers excluded, $\ge5L_D$ from contacts); saturation $j\to qnv_{sat}$ is reproduced | yes |
| A5 | MOS capacitor with $Q_f$ and $Q_{it}$: numerical $\psi_s(V_G)$, Richardson-extrapolated over two meshes, matches the exact reference (stage-8.1 Gauss relation with $Q_f+Q_{it}(\psi_s)$ added; independent root-find) to ≤ 2e-3 relative where $|\psi_s|\ge0.05$, and ≤ 1e-4 absolute otherwise | yes |
| A6 | Analytic Jacobian vs central FD with all three models active, per row, ≤ 1e-6 ($\epsilon=10^{-6}$; FD step chosen for the expm1/saturation scales), at random states including $\delta=10^{-24}$ | yes |
| A7 | 2D MOSFET (stage-8.1 C3 geometry) with all models on, at study values: converges; source + drain + body current balance ≤ 1e-8 relative where $|I_D|\ge10^{-6}$ | yes |
| R1 | Transfer curve with and without $Q_{it}$ and with field-dependent mobility (threshold shift, saturation) | reported |

Stage 9a passes if A1–A7 hold.
