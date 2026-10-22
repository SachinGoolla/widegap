# Stage 7 result: BLOCKED on provenance (2026-10-21)

Raw data: `results/stage7_run1_blocked/` (commit `8f1c70f`, clean tree). Registry:
`src/widegap/materials/data/sic4h.json`. Search logs:
`docs/literature/stage7/`.

| Criterion | Result | Detail |
|---|---|---|
| M1 provenance, required set R | **fail** | Missing: $E_g$(300 K), $N_c$, $N_v$, N-donor ionization energy. No readable primary source was found for these (§ below) |
| M2, M3, M5 | not evaluable | All need the missing parameters |
| M4 MMS with incomplete ionization | pass | order 2.00 (synthetic levels) |
| M6 Caughey–Thomas identities | pass | registered electron $\parallel c$ fit |

**What was established** (31 registry entries). Each one was read in full
by one agent and independently re-verified against the text by a second:

- **Permittivity:** $\varepsilon_\perp=9.77\pm0.01$ and
  $\varepsilon_\parallel=10.20\pm0.05$ (Li et al., IEEE J. Microwaves 2024;
  millimetre-wave, used as a quasi-static proxy).
- **Al acceptor:** $\Delta E_A=265\,\text{meV}-3.6\times10^{-5}\,\text{meV cm}\,N_A^{1/3}$
  (Koizumi, Suda & Kimoto, JAP 2009, Eq. 7).
- **Degeneracy factors:** $g_A=4$ and $g_D=2$ (conventions as used in read
  sources).
- **Mobility:** electron and hole Hall-mobility fits, $\parallel c$ and
  $\perp c$ (Ishikawa thesis, Kyoto 2025). The hole $\perp c$ value $\mu_{min}$
  is 40 in the thesis, and the unreadable journal version reportedly gives
  20. This conflict is unresolved.
- **Effective masses:** first-principles GW electron band-edge masses
  (Deng, Yang & Pi, PRB 2023; arXiv v1).

**Why blocked.** The readable sources for the missing parameters are
secondary, either reviews or compilations such as the ioffe.ru archive. Or
they are 0 K first-principles values that are not 300 K quantities, or
textbook values quoted in theses. The pre-registered rule forbids
substituting them, and nothing was substituted.

**What would unblock it** (DOIs checked against Crossref). Each needs a
legitimately obtained full text that is read and independently
re-verified:

- $N_c$, $N_v$ / density-of-states masses: C. Persson & U. Lindefelt, JAP
  82, 5496 (1997), [doi:10.1063/1.365578](https://doi.org/10.1063/1.365578);
  JAP 83, 266 (1998),
  [doi:10.1063/1.366681](https://doi.org/10.1063/1.366681). The
  temperature dependence of the hole mass at 300 K needs care, since the
  reviewed literature reports $m_{dh}$ more than doubling by 200 K.
- N donor (h and k sites): W. Götz et al., JAP 73, 3332 (1993),
  [doi:10.1063/1.352983](https://doi.org/10.1063/1.352983); I. G. Ivanov et
  al., PRB 67, 165211 (2003),
  [doi:10.1103/PhysRevB.67.165211](https://doi.org/10.1103/PhysRevB.67.165211).
- $E_g$(300 K): a primary temperature-dependent absorption or
  photoluminescence measurement. Candidates are listed in search log A.

Alternatively, the project owner may decide to accept a documented class of
secondary sources (for example a named review). That would be a new,
disclosed rule change with weaker provenance, and it is not made here.

**Not claimed:** any 4H-SiC device simulation, any SiC-realistic
$n_i/N_0$ result (M5 was not run), or that the registered mobility fits
(Hall mobilities) are drift mobilities.

## Run 2 (2026-10-21): all required parameters sourced; gate FAILED on M5(b) only

The user supplied four full texts: Persson & Lindefelt 1998, Götz et al.
1993, Ivanov et al. 2003 and Grivickas et al. 2007. A fifth file was a
course handout and was rejected. Values were read by the lead and
independently re-verified by a second agent. The registry now has 51
entries, all usable. Raw data: `results/stage7_run2/`.

| Criterion | Result | Value |
|---|---|---|
| M1 provenance | pass | required set complete; 51 entries pass the loader rules |
| M2 closed-form neutrality, 200/300/500 K | pass | worst relative error 1.7e-14 |
| M3 ionization limits | pass | |
| M4 MMS with incomplete ionization | pass | order 2.00 |
| M5(a) convergence at realistic $n_i/N_0$ | pass | $\delta=1.0\times10^{-24}$, 52 Newton iterations |
| M5(c) junction field vs exact first integral | pass | 9.8e-7 relative (threshold 2e-3) |
| M5(b) bulk potentials vs neutrality roots | **fail** | 1.7e-4 $V_T$ (threshold 1e-10) |
| M6 mobility identities | pass | |

**Diagnosis of M5(b).** The test geometry is too short, and the
size of the error was predictable. At $N_A=10^{17}$ cm⁻³ only 7.7% of the Al
is ionized in the neutral p bulk. That agrees with the hand estimate
$p^2/(N_A-p)=(N_v/g_A)e^{-\Delta E_A/kT}$, which gives
$p\approx8\times10^{15}$ cm⁻³. The p-side screening length is
$(p+N_A t(1-t))^{-1/2}\approx0.82\,L_D$. With the p contact only 10 $L_D$
from the junction, the space-charge tail at the check point $x=-9$ is about
$e^{-7.5/0.82}\approx10^{-4}$, which matches the observed 1.7e-4. The
solver and the first-integral reference are consistent to 1e-6.

**Reported (R7a):**
- $E_g(300\,\text{K})=3.2626$ eV.
- $N_c=1.72\times10^{19}$ cm⁻³ (1.03× the value Grivickas et al. quote from Persson & Lindefelt).
- $N_v=3.78\times10^{19}$ cm⁻³, parabolic, which is an underestimate.
- $n_i(300\,\text{K})=1.0\times10^{-8}$ cm⁻³.
