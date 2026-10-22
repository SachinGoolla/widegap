# Stage 8.1 result: gate PASSED (2026-10-21)

Raw data: `results/stage8_1/` (commit `24b6579`, clean tree). Implemented by
a Sonnet subagent from the pre-registered spec, then reviewed. G2 and G6 were
amended before the official run, disclosed in spec §6.

| ID | Result | Value (development run; the official run passes identically) |
|---|---|---|
| G1 MOS capacitor $\psi_s(V_G)$ vs exact reference, accumulation to strong inversion | pass | ≤ 8.4e-5 relative (≤ 1.3e-8 absolute near flat band) |
| G2 Gauss conservation | pass | ≤ 4.7e-13 absolute, about 1e-14 relative |
| G3 2D capacitor = 1D | pass | 1.4e-14 |
| G4 MOSFET convergence | pass | all biases, ≤ 23 Newton iterations per continuation target |
| G5 source + drain + body current balance | pass | ≤ 2.9e-9 relative (limit 1e-8; the floor is about 3× below the limit) |
| G6 source/drain mirror symmetry | pass | 5.4e-16 |
| G7 $I_D$ mesh convergence | pass | observed order 2.12; finest mesh within 0.10% of Richardson |
| G8 masked Jacobian vs FD | pass | ≤ 5.4e-10 |

**Reported (R1).** The threshold from linear extrapolation is about 154
$V_T$ in study units. The peak transconductance is 0.77× the
gradual-channel value, and $I_D/I_{GCA}$ is 0.53–0.65 in the on state. The
expected short-channel and finite-inversion-layer departures are of this
size; they are not validated. R2 figure:
`figures/stage8_1_mosfet_on_state.png`. The inversion layer is the thin
sheet at the oxide interface.

**Verified:** 2D MOS machinery. That covers oxide and semiconductor regions,
the gate contact, the zero-carrier-flux interface, the exact MOS electrostatics
from accumulation to inversion, and 2D MOSFET currents (conservation,
symmetry, convergence) at realistic SiC $n_i/N_0$.

**Not claimed:** a 4H-SiC MOSFET prediction. $\varepsilon_{ox}/\varepsilon_s=0.38$
and $V_{FB}$ are study parameters. There are no interface traps or fixed
charge, no channel-mobility degradation, mobility is constant, and
ionization is complete.
