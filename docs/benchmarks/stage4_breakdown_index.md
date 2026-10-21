# Stage 4: field-based breakdown index in 1D (non-punch-through and punch-through diodes)

Status: **pre-registered** (2026-10-20). No stage-4 code exists at the time
of this commit. Tolerances were sized analytically against round-off floors
and known approximation errors *before* any code, after the lessons of
stages 2 and 3.

## 1. What is computed, and what is not claimed

The breakdown-voltage index of Zhu & Andrei (2019), eq. (1)
([doi:10.1063/1.5113764](https://doi.org/10.1063/1.5113764)):
$B=\max_\Omega|\nabla\varphi|/E_{crit}$. The field-based "BV" is the reverse
bias at which $B=1$. In this stage, $E_{crit}$ is a **constant chosen for the
numerical study** ($E_{crit}=15\,V_T/L_D$ in scaled units). It is not a SiC
or Si material value.

**Not claimed:** avalanche breakdown. No impact ionization is modelled. The
field-based BV is a numerical proxy, verified here as a *numerical quantity*
and against depletion-approximation physics. Validation against avalanche
breakdown requires an impact-ionization model and is out of scope.

## 2. Model and discretization

Reverse-bias electrostatics in the zero-current approximation. The
quasi-Fermi potentials are flat, $\phi_n=0$ (cathode) and $\phi_p=-V_R$
(anode), so

$$
-\psi''=\delta\,(e^{\phi_p-\psi}-e^{\psi-\phi_n})+C .
$$

Anode $\psi=\psi_{eq}(C_{anode})-V_R$; cathode $\psi=\psi_{eq}(C_{cathode})$.
The space charge is exact up to minority-carrier terms of size
$\lesssim\delta^2/N$, which are negligible in reverse bias. This is the
stage-2 operator extended with flat quasi-Fermi levels. With
$\phi_n=\phi_p=0$ it is bit-identical to stage 2.

Newton with the analytic Jacobian, the stage-2 step limiter, and bias
continuation in steps of $1\,V_T$. Convergence uses the stage-3b
round-off-relative tolerance:
$\lVert F\rVert_\infty\le10^{2}\epsilon_{mach}\max_i(\sum_j|J_{ij}\psi_j|+V_i|\rho_i|)$,
together with a final $\lVert\Delta\psi\rVert_\infty\le10^{-8}$.

The element field is $E_e=-(\psi_{e+1}-\psi_e)/h_e$, and
$B=\max_e|E_e|/E_{crit}$. BV is found by bracketing between integer
$V_R$ and then Brent's method (`xtol` $10^{-10}$). Each evaluation is warm
started from the bracketing solution.

## 3. Cases ($\delta=10^{-3}$, all units scaled)

- **B1 (non-punch-through):** p⁺ ($a=100$) on $[-5,0]$, n ($b=1$) on
  $[0,60]$. The depletion edge at BV is about 15, well inside.
- **B2 (punch-through):** p⁺ ($a=100$) on $[-5,0]$, n ($b=1$) on $[0,8]$,
  n⁺ ($b=100$) on $[8,13]$. Punch-through is expected near
  $V_R\approx16$.
- Meshes: uniform, $h\in\{0.05,0.025,0.0125\}$; $h=0.025$ is the primary
  mesh. $V_R=0,1,\dots,120$.

Depletion-approximation references (S. M. Sze, 2nd ed., ch. 2; derivation
to be added to `docs/equations.md` §1.4). B1: $E_{max}^2=\frac{2ab}{a+b}(V_{bi}+V_R-2)$,
$V_{bi}=\psi_{eq}(b)-\psi_{eq}(-a)$, giving
$BV_{dep}=E_{crit}^2\frac{a+b}{2ab}-V_{bi}+2$. B2 after punch-through:
$E_{max}=(V_{bi}^{p^+n^+}+V_R)/W+bW/2$ with $W=8$.

## 4. Acceptance criteria

| ID | Case | Criterion | Gating |
|----|------|-----------|--------|
| G1 | B1, B2 | BV extraction: $|B(BV)-1|\le10^{-9}$ | yes |
| G2 | B1, B2 | mesh convergence: $|BV_{h=0.025}-BV_{h=0.0125}|/BV\le1\%$ | yes |
| G3 | B1, B2 | $BV$ within 10% of the depletion-approximation BV | yes |
| G4 | B1 | least-squares slope of $E_{max}^2$ vs $V_R$ on $V_R\in[20,80]$ within 3% of $2ab/(a+b)$ | yes |
| G5 | B2 | least-squares slope of $E_{max}$ vs $V_R$ on $V_R\in[40,100]$ within 10% of $1/W$ | yes |
| G6 | B1 at $V_R=50$ | least-squares slope of $|E_e|$ vs $x$ for element midpoints in $[2,x_n-3]$ equals $-b$ within 1%, where $x_n$ is the first $x>0$ with $\lvert E\rvert<0.01E_{max}$ | yes |
| G7 | B1, B2 at $V_R=50$ and $V_R=0$ initial guess | analytic Jacobian vs central differences ($\epsilon=10^{-7}$, 5 seeded directions): relative error ≤ 1e-6 | yes |
| G8 | all solves | Newton converges within 60 iterations | yes |
| G9 | B1 | punch-through absent: $|E|$ at the cathode-side element below $10^{-3}E_{max}$ at BV | yes |
| R1 | B1 at $V_R\in\{2,5\}$ | $\max|\psi_{QF}-\psi_{DD}|$ against the stage-3 drift-diffusion solver | reported |
| R2 | B1, B2 | BV on three meshes, and observed order | reported |
| R3 | B1, B2 | $BV/BV_{dep}$ | reported |

## 5. Gate

Stage 4 passes if G1–G9 hold. Otherwise stage 5 (adjoint sensitivities)
does not start.
