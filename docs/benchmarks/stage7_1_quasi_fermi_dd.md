# Stage 7.1: drift-diffusion at realistic 4H-SiC $n_i/N_0$ (quasi-Fermi variables)

Status: **pre-registered** (2026-10-21). No stage-7.1 code exists at the
time of this commit. Tolerances and the regimes in which each check gates
were derived from round-off analysis before any code (§3).

## 1. Why stage 3 cannot be used as-is

With $\delta=n_i/N_0\approx10^{-24}$ (stage 7), minority densities in
neutral regions are about $\delta^2/N\sim10^{-48}$. The stage-3 unknowns
are $(\psi,n,p)$ and its tolerance takes a block-wise maximum over rows,
so minority rows are judged on the majority-carrier scale. Those rows
"converge" without resolving any minority physics. This is reported as R1
below.

## 2. Formulation (same physics and discretization as stage 3)

The unknowns per node are $X=(\psi,\phi_n,\phi_p)$, with
$n=\delta e^{\psi-\phi_n}$ and $p=\delta e^{\phi_p-\psi}$. The scaled
equations are those of stage 3: Poisson, $j_n'=R$, $j_p'=-R$. In these
variables the Scharfetter–Gummel fluxes are algebraically identical to
stage 3, written without cancellation:

$$
j_n=\tfrac{\mu_n}{h}B(\Delta)\,\delta e^{\psi_{b}-\phi_{n,a}}\operatorname{expm1}\!\big(-(\phi_{n,b}-\phi_{n,a})\big),\qquad
j_p=-\tfrac{\mu_p}{h}B(-\Delta)\,\delta e^{\phi_{p,a}-\psi_{b}}\operatorname{expm1}\!\big(\phi_{p,b}-\phi_{p,a}\big),
$$

where $\Delta=\psi_b-\psi_a$ on element $(a,b)$. SRH uses
$np-\delta^2=\delta^2\operatorname{expm1}(\phi_p-\phi_n)$. Newton uses the
analytic Jacobian, the stage-2 step limiter on $\psi$, and bias
continuation in steps of $1\,V_T$. Ohmic contacts:
$\phi_n=\phi_p=V_c$ and $\psi=\psi_{eq}+V_c$.

**Per-row tolerance.** Converged when, for every interior row $i$,
$|F_i|\le10^2\epsilon_{mach}S_i$ with
$S_i=\sum_j|J_{ij}|\max(|X_j|,1)+V_i|\text{volume terms}_i|$, and the
final $\lVert\Delta X\rVert_\infty\le10^{-8}$. The difference from stage 3b
is that this test is per row, not a block-wise maximum.

## 3. Round-off analysis (what is resolvable)

- **Minority currents.** Across a minority region $\phi$ varies by
  $O(V)$, so $\Delta\phi$ per element is about $10^{-2}$ to $10^{-1}$. The
  relative resolution is about $\epsilon|\phi|/|\Delta\phi|\sim10^{-12}$,
  which is well resolved at any bias.
- **Majority currents.** In neutral majority regions
  $\Delta\phi\approx jh/(\mu N)$. The absolute resolution of $j$ is about
  $\epsilon|\phi|\mu N/h\sim10^{-13}$ per element, and with the $10^2$
  tolerance multiplier about $3\times10^{-11}$. Total-current constancy is
  therefore meaningful only for $|j|\gtrsim10^{-6}$ (scaled; about
  $10^{-2}$ A/cm² in 4H-SiC units). Practical forward conduction is far
  above this.
- **Exact discrete identity** ($R\equiv0$, 1D). Summing the SG electron
  flux telescopes, giving
  $e^{-\phi_{n,R}}-e^{-\phi_{n,L}}=j_n\sum_e h_e/(\mu_n\delta B(\Delta_e)e^{\psi_{b,e}})$.
  Given the computed $\psi$, this yields the exact discrete $j_n$. The hole
  identity is analogous. This checks *global* consistency of the minority
  current at any bias, including tiny currents.

## 4. Cases

- **Q1 (equivalence):** the stage-3 D3 diode ($\delta=10^{-3}$, $a=b=1$,
  $R\equiv0$, $[-40,40]$, $N=1024$), $V=0\ldots10\,V_T$.
- **Q2 (MMS):** the stage-3 D1 manufactured solution, $N=32\ldots512$.
- **Q3 (realistic):** $\delta=1.0\times10^{-24}$ (stage-7 $n_i/N_0$ at
  300 K, $N_0=10^{16}$ cm⁻³). $a=10$ (p), $b=1$ (n). Complete ionization,
  which isolates the numerics; incomplete-ionization DD is a later step.
  $R\equiv0$. Domain $[-30,60]$, $h=0.025$. $V=0,1,\dots,110\,V_T$, where
  $V_{bi}\approx112.8\,V_T$.
- **Q4:** Q3 with SRH, $\tau_n=\tau_p=10^3$. Reported.

## 5. Acceptance criteria

| ID | Case | Criterion | Gating |
|----|------|-----------|--------|
| T1 | Q1 | $\max|\psi-\psi_{\text{stage 3b}}|\le10^{-9}$ for all $V$; total current vs stage 3b within $10^{-6}$ relative for $V\ge3$ | yes |
| T2 | Q2 | observed orders of $\psi,n,p$ (last three) in [1.8, 2.2] | yes |
| T3 | Q1, Q3 | analytic Jacobian vs central FD ($\epsilon=10^{-7}$, 5 seeded directions, per-row relative to $\sum_j|J_{ij}||v_j|$) ≤ 1e-6, at an initial guess and a solution | yes |
| T4 | Q3 | Newton converges (per-row tolerance) at every bias within 60 iterations | yes |
| T5 | Q3 | at $V=0$: $\max|\psi-\psi_{\text{stage 2}}|\le10^{-9}$ and $\max(|\phi_n|,|\phi_p|)\le10^{-12}$ | yes |
| T6 | Q3, all $V\ge1$ | anode-side minority $j_n$ (first p-region element) vs the exact discrete identity, relative ≤ 1e-8. The same for cathode-side $j_p$ vs the hole identity | yes |
| T7 | Q3, $|\bar j|\ge10^{-6}$ | total-current constancy $\max_e|j_e-\bar j|/|\bar j|\le10^{-4}$ | yes |
| T8 | Q3, $V\in[95,102]$ | $d\ln j/dV\in[0.97,1.03]$ (low injection; $j$ resolvable) | yes |
| R1 | Q3 | stage-3 $(\psi,n,p)$ solver attempted at the same $\delta$; outcome documented | reported |
| R2 | Q4 | convergence and $d\ln j/dV$ with SRH | reported |
| R3 | Q3 | current in A/cm² using the stage-7 parameters ($\mu$ from the Ishikawa fit at $N_D$), as an illustration only | reported |

## 6. Gate

Stage 7.1 passes if T1–T8 hold. A pass means the quasi-Fermi formulation
is a verified, numerically sound replacement for stage 3 at realistic
4H-SiC $n_i/N_0$, within the resolvability limits of §3.
