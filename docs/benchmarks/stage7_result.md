# Stage 7 result: BLOCKED on provenance (2026-10-21)

Raw data: `results/stage7/` (commit `0b20f2b`, clean tree). Registry:
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
