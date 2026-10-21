# Stage 5: discrete adjoint doping sensitivities

Status: **pre-registered** (2026-10-20). No stage-5 code exists at the time
of this commit. The finite-difference (FD) error budget below was worked out
analytically before any code, following the lessons of stages 2–4.

## 1. Method (independent implementation after Zhu & Andrei 2019, §III)

Each verified forward solver is written as $F(X,C)=0$, where $X$ is the
nodal state and $C$ the nodal net doping in the scaled units of
`docs/equations.md`. **Design variables** are the interior nodal values
$C_i$. The contact values stay fixed, because they set the boundary
conditions. For an objective $J(X)$:

$$
F_X^{T}\lambda=\Big(\frac{\partial J}{\partial X}\Big)^{T},\qquad
\frac{dJ}{dC_i}=-\lambda^{T}\frac{\partial F}{\partial C_i}=\lambda_{\psi,i}\,V_i ,
$$

because $C$ enters only the Poisson rows, as $-V_iC_i$. $F_X$ is the
analytic Jacobian already verified in stages 2–4. One transposed linear
solve gives the sensitivity to every node.

**Sign convention.** Zhu & Andrei define sensitivities with respect to
$D=N_A-N_D=-C$, so $\gamma^{(D)}=-\gamma^{(C)}$. Both are exported.

## 2. Cases and objectives

| Case | Forward problem | Objective $J$ |
|---|---|---|
| A1 | stage 2, E2 symmetric junction, $N=1024$ | $\int\psi\,dx$ (trapezoid) |
| A2 | stage 3, D3 diode ($R\equiv0$), $N=1024$, $V=6\,V_T$ | mean element total current $\bar j$ |
| A2s | as A2 with SRH, $\tau=10^3$ | $\bar j$ (reported) |
| A3 | stage 4, B1 at $h=0.025$, $V_R=50$ | $B=\lvert E_{e^*}\rvert/E_{crit}$, $e^*$ = argmax element |
| A4 | stage 4, B1 at $h=0.025$ | field-based BV via the implicit function theorem: $\frac{dBV}{dC_i}=-\frac{\partial B/\partial C_i}{\partial B/\partial V_R}$ (ZA19 eq. 14), with $\partial B/\partial V_R=-\lambda^TF_{V}$ |

$B$ is differentiable only while the argmax element $e^*$ is unique and
does not change. Every FD evaluation in A3 and A4 must report the same
$e^*$ as the base state.

## 3. Finite-difference verification

**Directions** (seeded, 5 per case): $\delta C_i=\epsilon\,\max(|C_i|,1)\,m(x_i)$
for interior $i$. Here $m(x)=\sum_{k=1}^{4}a_k\sin(k\pi\xi)$ with
$\xi\in[0,1]$ the normalised position, $a_k\sim\mathcal N(0,1)$, rescaled
so that $\lVert m\rVert_\infty=1$. Smooth directions make the FD signal
add up over many nodes.

**FD forward solves:** the verified stage solver warm started from the base
solution, followed by **3 additional full Newton steps** ("polish") to reach
the round-off floor. BV in A4 uses Brent with `xtol` $10^{-12}$.

**Error budget.** For a central difference, truncation error is about
$\epsilon^2|J'''|$ and noise is about (state round-off error)/$\epsilon$.

- A1 and A3: $|J|\sim10^{1}$ and residual floors about $10^{-13}$. With
  $\epsilon=10^{-4}$, noise is about $10^{-9}$ and truncation about
  $10^{-8}$ relative.
- A2: $\bar j\approx2\times10^{-5}$ at $V=6$, and the current round-off
  error is about $10^{-13}$. With $\epsilon=10^{-3}$ and a directional
  signal of about $10^{-1}\epsilon\bar j$, noise is about $5\times10^{-5}$.
  Truncation is about $10^{-6}$.
- A4: Brent noise $10^{-12}$ against a signal of about $10^{-3}\cdot10$.
  With $\epsilon=10^{-3}$, noise is about $10^{-8}$ and truncation about
  $10^{-6}$.

| ID | Case | Criterion (relative error of the directional derivative, all 5 directions) | Gating |
|----|------|------|--------|
| S5a | A1 | ≤ 1e-6 at $\epsilon=10^{-4}$ | yes |
| S5b | A2 | ≤ 1e-3 at $\epsilon=10^{-3}$ | yes |
| S5c | A3 | ≤ 1e-6 at $\epsilon=10^{-4}$, and $e^*$ unchanged in every FD solve | yes |
| S5d | A4 | ≤ 1e-4 at $\epsilon=10^{-3}$, and $e^*$ unchanged | yes |
| S5e | A2 | cost: median wall time of one adjoint gradient (Jacobian assembly + transposed solve + product) ≤ 3× the median time of one Newton iteration (assembly + solve), over 20 repeats | yes |
| S5f | A1–A4 | single-node central FD at 5 seeded interior nodes vs adjoint | reported |
| S5g | A2s | directional check as S5b | reported |
| S5h | A2 | FD-gradient cost estimate ($2\times$ number of design variables × forward-solve time) vs adjoint time | reported |

## 4. Gate

Stage 5 passes if S5a–S5e hold. Otherwise stage 6 (optimization) does not
start.
