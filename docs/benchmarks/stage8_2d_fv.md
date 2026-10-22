# Stage 8: 2D finite volume on tensor-product meshes

Status: **pre-registered** (2026-10-21). No stage-8 code exists at the time
of this commit.

## 1. Discretization

The mesh is a rectangular tensor product with node coordinates
$x_0<\dots<x_{N_x}$ and $y_0<\dots<y_{N_y}$, which may be non-uniform. The
method is the vertex-centred box method: each node's control volume is its
dual cell. Materials such as $\varepsilon$ and doping are defined **per
primal cell**.

The scheme is an **edge-based core shared with 1D**. Each edge $e=(a,b)$
has a length $h_e$ and a coupling weight $w_e$, which is the dual-face
length. For a cell-wise coefficient $k$ (for example $\varepsilon$), the
edge conductance is $G_e=\sum_{c\ni e}k_c\,w_{e,c}/h_e$, summed over the
two primal cells adjacent to the edge, with $w_{e,c}$ the half-width of
cell $c$ normal to the edge. This makes the flux exact for piecewise-constant
coefficients whose interfaces lie on grid lines. Node volumes are
dual-cell areas. Nodal source and doping quantities are dual-cell
(control-volume) averages of the cell-wise data.

Boundary conditions:

- **Dirichlet:** row replaced by $X-g$.
- **Homogeneous Neumann (zero flux):** natural; boundary dual faces are
  simply omitted.

The physics is that of stages 2, 3 and 7.1, applied edge by edge. That
covers the linear Poisson equation, nonlinear Boltzmann equilibrium, and
quasi-Fermi drift–diffusion with Scharfetter–Gummel fluxes (expm1 form),
using the per-row round-off tolerance of stage 7.1.

## 2. Acceptance criteria

| ID | Case | Criterion | Gating |
|----|------|-----------|--------|
| P1 | Linear MMS on the unit square, $u=\sin\pi x\sin\pi y+xy$, $\varepsilon(x,y)=1+x^2+y^2$ (cell-centre value), Dirichlet, uniform $N\in\{16,32,64,128,256\}$ | max-norm observed orders, last three, in [1.9, 2.1] | yes |
| P2 | Dielectric interface: $\varepsilon=1$ for $x<0.5$ and $\kappa=2.5$ for $x>0.5$, $f=0$, Dirichlet $u=0$ / $u=1$ on the left / right, zero flux on top and bottom; repeat rotated by 90° | exact piecewise-linear solution reproduced to ≤ 1e-12 on $33\times17$ and $32\times16$ meshes. The interface is always on a grid line, as §1 requires | yes |
| P3 | P1 on graded meshes: geometric spacing with ratio 1.05 in $x$ and 0.95 in $y$, each level halving every spacing (same grading pattern, subdivided) | orders, last three, in [1.8, 2.2] | yes |
| P4 | P1 with zero flux on the top edge, $N=128$ | discrete balance: $\sum$ boundary fluxes $+\sum_iV_if_i=0$ to ≤ 1e-12 relative to $\sum_iV_i|f_i|$ | yes |
| P5 | Nonlinear 2D MMS for $-\nabla^2\psi=\delta(e^{-\psi}-e^{\psi})+s$, $\delta=10^{-3}$, $\psi=3\tanh(4(x-\tfrac12))\cos(\pi y/2)$ | orders, last three, in [1.9, 2.1] | yes |
| P6 | **Reduction to 1D, equilibrium:** stage-2 E2-type junction with $x$-dependent doping, zero flux at top and bottom, $N_y=8$, same $x$-mesh as a 1D stage-2 solve | every $y$-row equals the 1D solution to ≤ 1e-11 | yes |
| P7 | **Reduction to 1D, quasi-Fermi DD:** stage-7.1 Q3 diode ($\delta=10^{-24}$, $a=10$, $b=1$, $R\equiv0$), $N_y=4$, zero flux at top and bottom, $V\in\{0,50,100,105\}\,V_T$ | per-row $\psi,\phi_n,\phi_p$ equal the stage-7.1 1D solution to ≤ 1e-10 (abs, in $V_T$); terminal current per unit width equals the 1D current to ≤ 1e-9 relative wherever $|j|\ge10^{-6}$ | yes |
| P8 | 2D quasi-Fermi DD analytic Jacobian vs central FD (per-row, as in stage 7.1), random states, including $\delta=10^{-24}$ and SRH | ≤ 1e-6 | yes |
| P9 | The 1D code path through the shared core reproduces the committed stage-7.1 Q3 results | ≤ 1e-12 relative on currents at every bias | yes |
| R1 | Wall time and peak RSS of linear and quasi-Fermi solves vs number of nodes (up to $257^2$) | reported | no |
| R2 | A genuinely 2D illustration: a diode with a partial-width anode contact (current crowding), field and current maps | reported, not verified beyond P1–P9 | no |

## 3. Gate

Stage 8 passes if P1–P9 hold.

## 4. Amendments

### 2026-10-21: P9 compares only resolvable currents (before the official run, after a development run)

**Disclosure.** A development run, which is not an official result, showed
that P9 as written can never pass. For $V\lesssim25\,V_T$ the committed
stage-7.1 total current `j_mean` is round-off noise, about 1e-28, with an
element-to-element spread of 10–25. Its own pre-registration (stage 7.1
§3) says total current is resolvable only for $|j|\ge10^{-6}$. Reproducing
noise to 1e-12 is meaningless. P9 should have applied that rule.

**P9 (restated).** The 1D path through the shared core reproduces the
committed stage-7.1 Q3 results in two ways:

1. Minority edge currents (`jn` at the first element, `jp` at the last
   element) must agree to ≤ 1e-12 relative at every bias $V\ge1$.
2. The total current must agree to ≤ 1e-12 relative wherever
   $|j|\ge10^{-6}$.

The interpretations of P2 (meshes of 32×16 and 16×32 cells, so that the
interface lies on a grid line), P4 (P1 data with zero flux on top) and P5
(N = 16…256) used by the implementation are recorded as part of the
specification.
