# Stage 4 results

## Run 1, pre-registered criteria v1 (2026-10-20): gate FAILED (G2, G3, G4, G5)

Raw data: `results/stage4/`.

| Criterion | Result | Value |
|---|---|---|
| G1 BV extraction \|B(BV)−1\| ≤ 1e-9 | pass | ≤ 2.4e-14 |
| G2 mesh change 0.025→0.0125 ≤ 1% | **fail** | B1 1.6%, B2 1.3%; observed order 0.79 (B1), 0.77 (B2) |
| G3 BV within 10% of depletion approximation | **fail** (B1) | B1 86.30 vs 97.20 (ratio 0.888); B2 62.80 vs 64.97 (0.967) |
| G4 B1 slope of $E_{max}^2$ vs $V_R$ within 3% of 2ab/(a+b) = 1.980 | **fail** | 1.628 |
| G5 B2 slope of $E_{max}$ vs $V_R$ within 10% of 1/W = 0.125 | **fail** | 0.1057 |
| G6 field slope in the n depletion region = −1 (±1%) | pass | −0.997 |
| G7 analytic Jacobian vs FD | pass | ≤ 1.8e-8 |
| G8 Newton convergence | pass | ≤ 13 iterations |
| G9 no punch-through in B1 | pass | cathode-side field / E_max = 6.6e-8 |

### Diagnosis: the references were wrong, not the solver

The analysis below was done on development runs before the official run,
which then reproduced the predicted failures.

1. **The depletion approximation is inadequate at this bias scale.** The
   heavy p⁺ side ($a=100$) is not depleted. It is compressed within its
   Boltzmann (Debye) tail, and its majority holes spill across the junction
   ($p(0)\approx19$ at $V_R=50$). The exact first integrals on the two sides
   give $E(0)^2/2=\Delta\psi_n-1+p(0)=a\Delta\psi_p-a+p(0)$. From these,
   $d\Delta\psi_n/dV=a/(a+1)$ and $dp(0)/dV=-p(0)\,d\Delta\psi_p/dV$, so
   $d(E^2)/dV=2\big(\tfrac{a}{a+1}-\tfrac{p(0)}{a+1}\big)\approx2(0.990-0.189)=1.60$.
   The simulation converges to **1.603**. The depletion value 1.98 ignores
   $p(0)$. That term is negligible at device-scale voltages ($10^4\,V_T$) but
   not at the $10^2\,V_T$ used here. The punch-through slope (G5) has the
   same cause at both heavy layers.
2. **The element-midpoint maximum field converges at first order.** The
   peak element straddles a region where $|E'|=|\rho|\approx p(0)$, so
   $\max_e|E_e|=|E(0)|-O(\rho(0)h)$. A 1% change between $h=0.025$ and
   $0.0125$ is not attainable at first order with $\rho(0)\approx20$.

Under the project rules the gate stays failed. The criteria are
re-registered in `docs/benchmarks/stage4b_first_integral.md`.

## Run 2, re-registered criteria 4b (2026-10-20): gate PASSED

Raw data: `results/stage4b/` (commit `5af89d4`, clean tree).

| Criterion | Result | Value |
|---|---|---|
| G1 extraction | pass | \|B−1\| ≤ 2.4e-14 on all 5 meshes |
| G2′ BV mesh convergence | pass | B1 orders 0.79, 0.90, 0.95; finest within 0.47% of Richardson BV 83.38. B2 orders 0.77, 0.89, 0.94; within 0.36% of 61.16 |
| G3′ junction field vs exact first integral | pass | B1 3.0e-5 to 5.9e-5; B2 1.4e-5 to 4.9e-5 (threshold 2e-3) |
| G6′ depletion-region field slope | pass | −0.99990 |
| G7, G8, G9 | pass | as in run 1 |
| depletion approximation (reported) | n/a | extrapolated BV / BV_dep = 0.858 (B1), 0.941 (B2) |

**Verified:** the field-based breakdown index $B$ and its BV extraction, as
*numerical quantities*, for 1D abrupt multilayer diodes in the zero-current
approximation, with a constant $E_{crit}=15\,V_T/L_D$ chosen for the study.

**Not claimed:** avalanche breakdown, any physical BV, or any SiC value.
The max-element-field definition converges only at first order. Adjoint
sensitivities of $B$ (stage 5) will have to account for this
non-smoothness: $\max$ is not differentiable, and ZA19 widens it to nodes
with $B>0.9\max B$.
