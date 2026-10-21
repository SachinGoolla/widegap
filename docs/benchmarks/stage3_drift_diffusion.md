# Stage 3: 1D steady-state drift-diffusion (Scharfetter–Gummel)

Status: **pre-registered** (2026-10-20). No stage-3 code exists at the time
of this commit. This document applies the lesson of stage 2: the analytic
Jacobian is verified directly (D2, gating), and the asymptotic convergence
rate is only reported.

## 1. Model (scaled units of `docs/equations.md` §1.2)

Unknowns at nodes: $\psi$, $n$, $p$ (densities in units of $N_0$).

$$
-\psi''=p-n+C,\qquad j_n'=R,\qquad j_p'=-R,
$$
$$
j_n=\mu_n\,(n'-n\psi'),\qquad j_p=-\mu_p\,(p'+p\psi'),\qquad
R=\frac{np-\delta^2}{\tau_p(n+\delta)+\tau_n(p+\delta)} .
$$

Currents are in units of $q\mu_0V_TN_0/L_D$, time in units of
$L_D^2/(\mu_0V_T)$, and mobilities in units of $\mu_0$. The signs follow
from $J_n=q\mu_nnE+qD_n\nabla n$, $J_p=q\mu_ppE-qD_p\nabla p$,
$\nabla\cdot J_n=qR$ and $\nabla\cdot J_p=-qR$, with the Einstein relation.
Steady state, isothermal, Boltzmann statistics, complete ionization,
constant mobilities. "No recombination" means $R\equiv0$.

**Ohmic contacts** (Dirichlet): $n$ and $p$ take their charge-neutral
equilibrium values, and $\psi=\psi_{eq}+V_c$, where $V_c$ is the contact
voltage. The anode (left, p side) is at $V$ and the cathode is at 0.

## 2. Discretization

- Box method, as in stages 1–2.
- Scharfetter–Gummel fluxes (SG69,
  [doi:10.1109/T-ED.1969.16566](https://doi.org/10.1109/T-ED.1969.16566)).
  With $B(t)=t/(e^t-1)$ and $\Delta=\psi_{i+1}-\psi_i$:
  $j_n=\frac{\mu_n}{h}[B(\Delta)n_{i+1}-B(-\Delta)n_i]$ and
  $j_p=-\frac{\mu_p}{h}[B(-\Delta)p_{i+1}-B(\Delta)p_i]$.
- $R$ is lumped at nodes.
- Fully coupled Newton with the analytic Jacobian. The max-norm step limiter
  of stage 2 applies to $\psi$; $n$ and $p$ updates are scaled by the same
  factor.
- Bias continuation from the stage-2 equilibrium solution, in steps of at
  most $1\,V_T$.
- Convergence: $\lVert F_\psi\rVert_\infty\le10^{-12}$ and
  $\lVert F_{n}\rVert_\infty,\lVert F_{p}\rVert_\infty\le10^{-14}$
  (continuity rows carry currents of order $\delta^2$).

## 3. Cases

- **D1 (coupled MMS)** on $[0,1]$ with $\mu_n=\mu_p=\tau_n=\tau_p=1$,
  $\delta=0.1$, $C=0$, and exact
  $\psi=2x+0.3\sin(2\pi x)$, $n=\exp(0.5\sin\pi x)$, $p=1/(1+x)$. Sources
  for all three equations are derived symbolically (sympy). Uniform
  $N\in\{32,64,128,256,512\}$.
- **D3 (symmetric pn diode)**: $a=b=1$, $\delta=10^{-3}$, $\mu_n=\mu_p=1$,
  $R\equiv0$, domain $[-40,40]$, $N=1024$, junction at 0. Biases
  $V\in\{-5,-2,0,1,2,\dots,10\}\,V_T$.
- **D4 (with SRH)**: D3 with $\tau_n=\tau_p=10^{3}$. Reported only.

## 4. Acceptance criteria

| ID | Case | Criterion | Gating |
|----|------|-----------|--------|
| D1 | D1 | max-norm observed orders of $\psi$, $n$ and $p$: last three each in [1.8, 2.2] | yes |
| D2 | D1 at $N=128$; D3 at $V=4$ | analytic Jacobian vs central differences ($\epsilon=10^{-7}$, 5 seeded directions, at the solution and at the initial guess): relative error ≤ 1e-6 | yes |
| D3a | D3, $V=0$ | $\max|\psi-\psi_{\text{stage 2}}|\le10^{-10}$ and $\max(|j_n|,|j_p|)\le10^{-12}$ | yes |
| D3b | D3, all $V\ne0$ | total current constant across elements: $\max_e|j_e-\bar j|/|\bar j|\le10^{-8}$ | yes |
| D3c | D3, $V=4\ldots8$ | ideality $d\ln j/dV$ (successive biases) within [0.97, 1.03] | yes |
| D3d | D3, $V=2\ldots8$ | $j/j_{\text{short}}\in[0.8,1.2]$, with the short-diode formula $j_{\text{short}}=\delta^2(e^V-1)\big(\mu_n/(aW_p)+\mu_p/(bW_n)\big)$ and neutral widths from the depletion approximation with $(V_{bi}-V-2)$ | yes |
| D5 | all runs | $n>0$ and $p>0$ at every node | yes |
| D6 | all runs | Newton converges within 60 iterations at every continuation step | yes |
| D3e | D3, $V<0$ | $j(V)/(-j_{\text{short}}(V))$ | reported |
| D4 | D4 | ideality for $V=1\ldots4$ (expected to rise toward 2 from depletion-region recombination) | reported |
| D7 | D3 | same diode in DEVSIM (Apache-2.0, S22), in physical units chosen so that the scaled problem is identical: ratio of currents for $V=2\ldots8$ | reported. DEVSIM is not a CI dependency |
| D8 | all | quadratic-rate ratios above the round-off floor (stage-2b definition) | reported |

Wall time and Newton iteration counts are recorded.

## 5. Gate

Stage 3 passes if D1, D2, D3a–D3d, D5 and D6 hold. Otherwise stage 4 does
not start.

## 6. Amendments

### 2026-10-20: tolerances vs the round-off floor (before any stage-3 code ran)

Analysis before implementation shows two thresholds in §2 and §4 are below
the double-precision floor.

In the neutral regions the SG flux is the small difference of two terms of
size $\mu B(\pm\Delta)\,c/h$. That is about 13 for D3 ($h=0.078$, $c\approx1$)
and up to about 800 for D1 at $N=512$. The absolute floor of a continuity row
is therefore about $\epsilon_{mach}$ times these, between $6\times10^{-15}$
and $2\times10^{-13}$, so a fixed $10^{-14}$ tolerance is not attainable on
fine meshes. Element-to-element variation of the total current equals the
summed continuity residuals, so D3b at $10^{-8}$ relative is not attainable
where $j\lesssim10^{-6}$ (low bias, $\delta=10^{-3}$).

Restated before any code exists:

- **Convergence (§2, and D6):** for each block $b\in\{\psi,n,p\}$,
  $\lVert F_b\rVert_\infty\le\tau_b=10^{2}\,\epsilon_{mach}\,S_b$. Here $S_b$
  is the maximum over rows of the sum of absolute values of the terms
  entering that row (both element fluxes and the volume term), evaluated at
  the current iterate. The final Newton update must also satisfy
  $\lVert\Delta\psi\rVert_\infty\le10^{-8}$.
- **D3b:** $\max_e|j_e-\bar j|/|\bar j|\le10^{-6}$, evaluated for
  $|V|\ge3\,V_T$, where $|\bar j|$ exceeds the round-off scale by more than
  $10^{6}$.

All other criteria are unchanged.
