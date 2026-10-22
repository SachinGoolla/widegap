# Stage 8 result: gate PASSED (2026-10-21)

Raw data: `results/stage8/`. Implemented by a Sonnet subagent from the
pre-registered spec, then reviewed. Criterion P9 was amended before the
official run, disclosed in spec §4: a development run showed the original
wording compared round-off noise.

| ID | Result | Value (development run; the official run is identical in criteria) |
|---|---|---|
| P1 MMS, uniform | pass | orders 1.997, 2.000, 2.000 |
| P2 dielectric interface (x and rotated) | pass | 2.3e-15 / 2.1e-15 |
| P3 MMS, graded meshes | pass | orders 1.9994–1.99999 |
| P4 discrete conservation | pass | 1.5e-14 relative |
| P5 nonlinear MMS | pass | orders 2.003, 2.000, 2.000 |
| P6 2D equilibrium = 1D stage 2, per row | pass | 1.5e-13 |
| P7 2D quasi-Fermi DD = 1D stage 7.1, δ = 1e-24 | pass | potentials ≤ 8e-15; resolvable current identical |
| P8 Jacobian vs FD | pass | ≤ 7.8e-7 at ε = 1e-7 (3.9e-9 at ε = 1e-5) |
| P9 1D path through the shared core = committed stage 7.1 (amended) | pass | edge currents ≤ 2.8e-14; resolvable totals exact |

**Reported:**
- R1 performance: 257² nodes for linear Poisson takes 0.38 s and about
  196 MiB. Quasi-Fermi DD at 97² nodes takes 2.7 s (8 Newton iterations).
  All CPU-only.
- R2: a 2D diode whose anode covers only part of one edge converges to
  14 $V_T$. Anode and cathode currents balance to better than 1e-12. The
  current crowding is mild (top/bottom $j_x$ ratio 1.12). Figure:
  `figures/stage8_current_crowding.png`, an illustration only.

**Verified:** 2D tensor-mesh box method (linear, nonlinear equilibrium,
quasi-Fermi DD), sharing one edge-based core with 1D.

**Not yet:**
- oxide or insulator regions, gate contacts, and interface charge (needed
  for a MOSFET);
- incomplete ionization and field-dependent mobility in DD;
- unstructured meshes.
