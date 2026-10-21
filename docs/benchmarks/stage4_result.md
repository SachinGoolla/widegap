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
