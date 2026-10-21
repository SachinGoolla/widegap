# Benchmark P: one-dimensional Poisson problem

Status: **pre-registered**. This document fixes the problem definitions, the
error measures and the acceptance criteria *before* any solver or PINN result
exists. It is committed before the solver code so that the Git history shows
the thresholds were not tuned to the results. Changing a criterion afterwards
requires a new, dated section below explaining why. The original text stays.

## 1. Why this benchmark

Every drift–diffusion device simulator (including the adjoint-optimization
framework of Zhu & Andrei, AIP Adv. 9, 095301 (2019),
[doi:10.1063/1.5113764](https://doi.org/10.1063/1.5113764)) solves a Poisson
equation for the electrostatic potential. Before any coupled or nonlinear
device problem is attempted, both the classical discretization and the PINN
have to solve the *linear* Poisson operator with a heterogeneous permittivity
correctly. The benchmark uses the method of manufactured solutions
(Roache, J. Fluids Eng. 124, 4 (2002),
[doi:10.1115/1.1436090](https://doi.org/10.1115/1.1436090)). The exact
solution is chosen first, and the source term follows from it analytically.

This benchmark is **not** a device simulation. It contains no carriers, no
doping, no SiC material parameters and no physical units. Its only claim is
the correctness of the numerical operators.

## 2. Problem statement (dimensionless)

Find $u:[0,1]\to\mathbb{R}$ such that

$$
-\frac{d}{dx}\!\left(\varepsilon(x)\,\frac{du}{dx}\right) = f(x),\qquad 0<x<1,
\qquad u(0)=g_0,\quad u(1)=g_1 .
$$

The variables relate to device electrostatics as follows: $u$ is the potential
scaled by a reference voltage, $\varepsilon$ is the relative permittivity
scaled by a reference value, and $f$ is the scaled space charge. The scaling
is not used in this benchmark.

### Case P1: smooth, constant permittivity (**gating**)

$$
\varepsilon(x)=1,\qquad u(x)=\sin(\pi x)+x,\qquad f(x)=\pi^2\sin(\pi x),\qquad g_0=0,\ g_1=1 .
$$

The non-homogeneous Dirichlet data check that the boundary conditions are
imposed correctly, not just a zero solution.

### Case P2: smooth, variable permittivity (**gating for the reference solver**, reported for the PINN)

$$
\varepsilon(x)=1+x^2,\qquad u(x)=\sin(\pi x)+x,
$$
$$
f(x)=-\big(\varepsilon u'\big)' = -2x\,\big(\pi\cos(\pi x)+1\big)+(1+x^2)\,\pi^2\sin(\pi x),\qquad g_0=0,\ g_1=1 .
$$

### Case P3: piecewise-constant permittivity, charge-free (**gating for the reference solver**)

$$
\varepsilon(x)=\begin{cases}1, & x<x_s\\ \kappa, & x>x_s\end{cases},\quad x_s=\tfrac12,\ \kappa=\tfrac52,\qquad f=0,\qquad g_0=0,\ g_1=1 .
$$

The exact solution is piecewise linear, with $u$ and the flux
$\varepsilon u'$ continuous at $x_s$. With slopes $a_1$ (left) and
$a_2=a_1/\kappa$ (right), and $a_1x_s + a_2(1-x_s) = 1$:

$$
a_1=\frac{1}{x_s+(1-x_s)/\kappa},\qquad
u(x)=\begin{cases}a_1x, & x\le x_s\\ a_1x_s + a_2(x-x_s), & x\ge x_s.\end{cases}
$$

$\kappa = 5/2$ is a dimensionless contrast. It is not a claimed material
value. This case models a dielectric interface. A vertex-centred finite-volume
scheme that uses the length-weighted *harmonic* mean of $\varepsilon$ over
each element reproduces this solution to round-off, **whether or not the
interface falls on a mesh node**. An arithmetic mean does not. The case
therefore tests the flux discretization directly.

### Case P4: steep internal layer (**stress test, reported, not gating**)

$$
\varepsilon(x)=1,\qquad u(x)=\tanh\!\Big(\frac{x-\tfrac12}{\delta}\Big),\ \delta=0.02,
\qquad f(x)=\frac{2}{\delta^{2}}\,\operatorname{sech}^2(s)\tanh(s),\ s=\frac{x-\tfrac12}{\delta},
$$

with $g_0=u(0)$ and $g_1=u(1)$. This mimics the sharp potential variation at a
junction, which is the regime where PINNs are known to have difficulties
(e.g. Wang, Teng & Perdikaris, SIAM J. Sci. Comput. 43, A3055 (2021),
[doi:10.1137/20M1318043](https://doi.org/10.1137/20M1318043)). It is reported
to characterize behaviour, not to pass or fail the project.

## 3. Reference solver (independent of the PINN)

- Vertex-centred finite volume (the box method used in drift–diffusion
  solvers) on a mesh of $N$ elements with nodes $x_0=0<\dots<x_N=1$.
- Element permittivity: the length-weighted harmonic mean
  $\bar\varepsilon_e = h_e \big/ \int_e \varepsilon^{-1}\,dx$. For smooth
  $\varepsilon$ the integral uses Gauss–Legendre quadrature. For piecewise
  constant $\varepsilon$ it is computed exactly from the sub-segment lengths.
- Source: $f$ at the node times the control-volume length (half of each
  adjacent element).
- Dirichlet values are imposed directly at $x_0$ and $x_N$. The resulting
  tridiagonal system is solved directly (banded LU).
- The solver code shares **only the problem specification** (the functions
  $\varepsilon$, $f$, $u$, $g_0$, $g_1$) with the PINN. It does not import any
  PINN code, and the PINN does not import any solver code.

## 4. Error measures

Errors are measured at the evaluation points: the mesh nodes for the solver,
and a fixed uniform grid of 2001 points that is separate from the training
points for the PINN.

- $e_\infty = \max_i |u_h(x_i)-u(x_i)|$
- relative discrete $L^2$: $e_{2,\mathrm{rel}} = \|u_h-u\|_2 / \|u\|_2$, with
  trapezoidal weights.
- observed order between successive uniform refinements:
  $p = \log(e_\infty(N)/e_\infty(2N))/\log 2$.

## 5. Acceptance criteria

### Reference solver (uniform meshes $N\in\{32,64,128,256,512\}$)

| ID | Case | Criterion |
|----|------|-----------|
| R1 | P1, P2 | each of the last three observed orders $p$ (max-norm) lies in $[1.9, 2.1]$ |
| R2 | P1 | $e_\infty \le 1\times10^{-5}$ at $N=512$ |
| R3 | P2 | $e_\infty \le 5\times10^{-5}$ at $N=512$ |
| R4 | P3 | $e_\infty \le 1\times10^{-12}$ for $N=64$ (interface on a node) **and** $N=63$ (interface inside an element) |
| R5 | — | a linear exact solution ($f=0$, $\varepsilon=1$) is reproduced to $\le 10^{-12}$; the discrete flux balance equals $\int_0^1 f\,dx$ to $\le 10^{-10}$ relative for P1 at $N=512$ |

R2 was derived before implementation. On a uniform mesh with $\varepsilon=1$
the scheme equals the 3-point central difference. Its truncation error is
$\tau=\tfrac{h^2}{12}u''''$, which gives a nodal error of about
$\tfrac{h^2\pi^2}{12}\approx 3.1\times10^{-6}$ at $h=1/512$. The thresholds
leave a factor of about 3 (P1) and about 15 (P2) of margin for constants that
were not estimated.

### PINN

Fixed configuration, chosen before any run:

- Network: fully connected, 3 hidden layers of 32 units, $\tanh$, float64.
- Dirichlet conditions enforced exactly by construction:
  $u_\theta(x)=g_0(1-x)+g_1x+x(1-x)\,N_\theta(x)$.
- Loss: mean squared strong-form residual
  $\big(-(\varepsilon u_\theta')' - f\big)^2$ at collocation points, with
  derivatives from automatic differentiation.
- Collocation: 256 points drawn uniformly at random in $(0,1)$ and redrawn
  every 100 Adam steps. The random stream is seeded.
- Optimization: Adam with learning rate $10^{-3}$ for 5000 steps, then
  L-BFGS (strong-Wolfe line search) for at most 1000 iterations on a fixed
  set of 512 collocation points.
- Seeds: 0, 1, 2, 3, 4. **Every seed is reported.** None is discarded or
  selected.

| ID | Case | Criterion |
|----|------|-----------|
| N1 | P1 | $e_{2,\mathrm{rel}} \le 1\times10^{-3}$ against the exact solution, **for all five seeds** (gating) |
| N2 | P1 | error against the reference solver ($N=512$) is reported alongside N1. Because the reference error is about $10^{-6}$, the two comparisons must agree to within that level |
| N3 | P2, P4 | $e_{2,\mathrm{rel}}$ reported for all seeds. Not gating |

Wall-clock time, CPU model, thread count and peak resident memory are recorded
for every run. Runtime is not an acceptance criterion.

### Gate

The first benchmark **passes** only if R1–R5 and N1 all hold. If it fails, the
project stops and the failure is reported. It does not move on to coupled
drift–diffusion or device problems.

## 6. Reproduction

```bash
python benchmarks/poisson1d/run.py --out results/poisson1d
```

This writes the raw per-case and per-seed JSON results and a reproducibility
manifest (Git commit, package versions, CPU, threads, seeds, command,
timestamp).
