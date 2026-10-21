# Stage 2: nonlinear equilibrium Poisson, abrupt pn junction (1D)

Status: **pre-registered** (2026-10-20). No stage-2 code exists at the time
of this commit. Amendments follow the rules in `docs/PROJECT_SPEC.md` §3.

## 1. Model

Thermal equilibrium, Boltzmann statistics, complete ionization, constant
permittivity, no traps. These are deliberate simplifications. Incomplete
ionization (essential for Al in 4H-SiC) belongs to stage 7.

With $\psi$ in units of the thermal voltage $V_T=k_BT/q$, length in units of
the extrinsic Debye length $L_D=\sqrt{\varepsilon V_T/(qN_0)}$ for a reference
doping $N_0$, and densities in units of $N_0$:

$$
-\psi'' = \delta\,(e^{-\psi}-e^{\psi}) + C(x),\qquad \delta = n_i/N_0,\qquad C=(N_D-N_A)/N_0,
$$

using $n=n_ie^{\psi}$, $p=n_ie^{-\psi}$ (quasi-Fermi level zero).

Ohmic contacts at both ends: $\psi$ is set to its charge-neutral value
$\psi_{eq}=\operatorname{asinh}(C/(2\delta))$. The physical inputs ($N_A$,
$N_D$, $n_i$, $T$, $\varepsilon$) are SI at the API boundary. The scaling
is documented in `docs/equations.md`, which is created by this stage.

## 2. Discretization (reference solver)

Box-method finite volume, reusing the stage-1 operator. The nonlinear term
is lumped at nodes. Newton's method uses the analytic Jacobian, a sparse
tridiagonal solve, and a documented step limiter: the maximum $|\Delta\psi|$
per iteration is capped at $\Delta_{max}=2$ (thermal units). Every limited
step is counted and reported. The initial guess is the local
charge-neutral potential.

## 3. Cases

- **E1 (MMS):** an added source $s(x)$ makes
  $\psi(x)=3\tanh(4(x-\tfrac12))$ exact on $[0,1]$, with $\delta=10^{-3}$
  and $C=0$ apart from $s$.
- **E2 (symmetric abrupt junction):** $N_A=N_D$, $C=\mp1$, $\delta=10^{-8}$,
  domain $[-L,L]$ with $L=40$ (Debye lengths), junction at 0.
- **E3 (asymmetric abrupt junction):** $N_A/N_D=100$, $\delta=10^{-8}$.
  Reported only.

## 4. Acceptance criteria

| ID | Case | Criterion | Gating |
|----|------|-----------|--------|
| S1 | E1 | observed max-norm order over N = 64…1024: last three in [1.9, 2.1] | yes |
| S2 | E1, E2 | Newton converges to $\lVert F\rVert_\infty\le10^{-12}$ in ≤ 40 iterations. The final 3 unlimited iterations show quadratic convergence, $\lVert F_{k+1}\rVert/\lVert F_k\rVert^2\le 10^{3}$ | yes |
| S3 | E2 | bulk neutrality: $|p-n+C|/|C|\le10^{-6}$ at nodes with $|x|\ge 20$ | yes |
| S4 | E2 | total charge $\int(p-n+C)\,dx$ is zero to $\le10^{-8}$ relative to $\int|C|dx$ | yes |
| S5 | E2 | peak field compared with the depletion approximation including the $2V_T$ correction (S. M. Sze, *Physics of Semiconductor Devices*, 2nd ed., Wiley, 1981, ch. 2): relative difference ≤ 10% | yes |
| S6 | E3 | same quantities as S3–S5 | reported |
| S7 | E2 | PINN (stage-1 architecture plus log-density or potential scaling, configuration fixed in an amendment before its first run): rel. L² of $\psi$ vs the FV solution at $N=1024$, all 5 seeds | reported (the PINN does not gate stage 2) |

Wall time, Newton iteration count, limited-step count and peak RSS are
recorded for every run.

## 5. Gate

Stage 2 passes if S1–S5 hold. Otherwise stage 3 does not start.

## 6. Amendments

### 2026-10-20: S2 evaluation mesh and residual definition (before any stage-2 code)

S2 did not fix the mesh it is evaluated on or the scaling of $F$. Both are
fixed here, before implementation:

- $F_i$ is the finite-volume balance at node $i$ (flux difference minus
  control volume times the source), in the scaled units of §1. It is the
  same quantity Newton drives to zero.
- S2 is evaluated on E1 at $N=256$ and on E2 at $N=1024$ (uniform meshes).
- "The final 3 unlimited iterations" means the last three Newton iterations
  in which the step limiter was not active.
- Peak field (S5) is $\max_e |\psi_{e+1}-\psi_e|/h_e$ over elements. The
  depletion-approximation value for one-sided doping magnitudes $a=N_A/N_0$,
  $b=N_D/N_0$ is $E_{max}=\sqrt{2ab(V_{bi}-2)/(a+b)}$ in units of
  $V_T/L_D$, with $V_{bi}=\psi_{eq}(n)-\psi_{eq}(p)$. The derivation is in
  `docs/equations.md`.
