# Literature and requirements report

Prepared 2026-10-20. This report sets out what the published sources
establish, what they leave unspecified, and what this repository must
implement and verify before it can make any device-level claim.

## 0. Source verification

Every bibliographic entry below was checked against the Crossref API
(title, authors, venue, year) on 2026-10-20. Open-source solver licences were
read from the GitHub API on the same date. Anything not checked this way is
marked **unverified**.

| Key | Reference | Verified |
|-----|-----------|----------|
| ZA19 | C. Zhu, P. Andrei, "Adjoint method for the optimization of insulated gate bipolar transistors," *AIP Advances* **9**, 095301 (2019). [doi:10.1063/1.5113764](https://doi.org/10.1063/1.5113764). CC BY 4.0. | Crossref ✓, full text read |
| BA26 | P. Brungi, P. Andrei, "Large-scale optimization of doping profiles in 4H–SiC metal–oxide–semiconductor field-effect transistors," *AIP Advances* **16**, 035014 (2026). [doi:10.1063/5.0320204](https://doi.org/10.1063/5.0320204). CC BY 4.0. | Crossref ✓, full text read from the publisher HTML (equations reconstructed from flattened HTML; see §2.3) |
| ZA17 | C. Zhu, P. Andrei, "Adjoint method for increasing the breakdown voltage and reducing the on-state resistance in wide band gap power transistors," *ECS Trans.* **77**, 133–152 (2017). [doi:10.1149/07706.0133ecst](https://doi.org/10.1149/07706.0133ecst) | Crossref ✓, **not read** |
| ZA21 | C. Zhu, P. Andrei, "Analytical estimation of breakdown voltage in insulated-gate bipolar transistors," *J. Comput. Electron.* **20**, 1202–1208 (2021). [doi:10.1007/s10825-021-01691-x](https://doi.org/10.1007/s10825-021-01691-x) | Crossref ✓, **not read** |
| CA98 | J. B. Casady et al., "4H-SiC power devices for use in power electronic motor control," *Solid-State Electron.* **42**, 2165–2176 (1998). [doi:10.1016/S0038-1101(98)00212-3](https://doi.org/10.1016/S0038-1101(98)00212-3) | Crossref ✓, **not read** |
| SG69 | D. L. Scharfetter, H. K. Gummel, "Large-signal analysis of a silicon Read diode oscillator," *IEEE Trans. Electron Devices* **16**, 64–77 (1969). [doi:10.1109/T-ED.1969.16566](https://doi.org/10.1109/T-ED.1969.16566) | Crossref ✓ |
| RPK19 | M. Raissi, P. Perdikaris, G. E. Karniadakis, "Physics-informed neural networks…," *J. Comput. Phys.* **378**, 686–707 (2019). [doi:10.1016/j.jcp.2018.10.045](https://doi.org/10.1016/j.jcp.2018.10.045) | Crossref ✓ |
| K21 | G. E. Karniadakis et al., "Physics-informed machine learning," *Nat. Rev. Phys.* **3**, 422–440 (2021). [doi:10.1038/s42254-021-00314-5](https://doi.org/10.1038/s42254-021-00314-5) | Crossref ✓ |
| WTP21 | S. Wang, Y. Teng, P. Perdikaris, "Understanding and mitigating gradient flow pathologies in physics-informed neural networks," *SIAM J. Sci. Comput.* **43**, A3055–A3081 (2021). [doi:10.1137/20M1318043](https://doi.org/10.1137/20M1318043) | Crossref ✓ |
| R02 | P. J. Roache, "Code verification by the method of manufactured solutions," *J. Fluids Eng.* **124**, 4–10 (2002). [doi:10.1115/1.1436090](https://doi.org/10.1115/1.1436090) | Crossref ✓ |
| S22 | J. Sanchez, "DEVSIM: A TCAD semiconductor device simulator," *J. Open Source Softw.* **7**, 3898 (2022). [doi:10.21105/joss.03898](https://doi.org/10.21105/joss.03898) | Crossref ✓ |

## 1. Which paper is the "SiC optimization paper"?

The PDF supplied for this project is **ZA19**, which optimizes a **silicon**
punch-through IGBT. The 4H–SiC work by the same group is **BA26**, which
applies the same adjoint framework to a 4H–SiC vertical MOSFET. In this
repository:

- ZA19 is the source of the adjoint sensitivity method as first applied to
  power devices, and of the breakdown-voltage index.
- BA26 is the source of the SiC device context and the comparison targets.
- Neither is yet a reproduction target (§4).

## 2. What the papers establish

### 2.1 Transport model (ZA19 §II.A; BA26 §II)

Both papers solve the nonlinear Poisson equation coupled with electron and
hole continuity equations (drift plus diffusion), discretized with finite
elements and solved with Newton's method in the authors' in-house simulator
RandFlux. RandFlux is not open source.

| Model component | ZA19 (Si IGBT) | BA26 (4H–SiC MOSFET) |
|---|---|---|
| Statistics | not stated explicitly | Boltzmann |
| Ionization | "complete ionization model" (wording as printed) | incomplete ionization of donors and acceptors |
| Mobility | Philips model | Caughey–Thomas with field saturation, plus Lombardi-type surface degradation |
| Recombination | Shockley–Read–Hall | SRH plus Auger |
| Impact ionization | Van Overstraeten | Van Overstraeten–De Man, "calibrated for 4H–SiC" |
| Band-gap narrowing | Slotboom | stated, model not named |
| Parameters | in the supplementary material, **not available here** | "taken from the Atlas device simulator" (commercial defaults); numerical values **not printed** |

### 2.2 Optimization method (ZA19 eqs. 1, 6–15)

- **Breakdown index.** $B=\max_\Omega |\nabla\varphi|/E_{crit}$ (eq. 1).
  BV is the collector–emitter voltage at which $B=1$. This field-based
  definition replaces impact-ionization breakdown inside the optimizer. For
  Si, ZA19 uses
  $E_{crit}\approx 4\times10^5/[1-\tfrac13\log_{10}((N_A+N_D)/10^{16})]$ V/cm.
  BA26 uses a constant $E_{crit}\approx3\times10^6$ V/cm for 4H–SiC.
- **Sensitivity functions.** $\gamma_{BV}(\mathbf r)$ and $\gamma_{ON}(\mathbf r)$
  are the Gâteaux derivatives of BV and $I_{ON}$ with respect to the net doping
  $D=N_A-N_D$ (eqs. 6–7).
- **Adjoint.** For a discrete residual $F(X,N_A,N_D)=0$, one linear solve
  with $F_X^{T}$ (eq. 10 for $I_{ON}$, eq. 13 for $B$) gives the sensitivity
  to every doping node. Its cost equals one Newton iteration.
- **BV sensitivity** via the implicit relation
  $\gamma_{BV}=-(\partial V_{CE}/\partial B)\,\gamma_B$ (eq. 14). Only the
  shape of $\gamma_B$ is needed.
- **Update.** $\Delta D^{(k)}=\varepsilon[a\gamma_{BV}+b\gamma_{ON}]$ with
  $a^2+b^2=1$ (eq. 15). The per-node change is capped at 50%. The BV index is
  widened to nodes with $B>0.9\max B$. Reported convergence is 3–10
  iterations.
- **Multi-objective handling.** ZA19 explicitly does *not* compute a Pareto
  front. It solves single-objective problems and inequality-constrained
  variants (eqs. 2–5). BA26 uses fixed weights $f=\alpha BV+\beta I_{ON}$.

### 2.3 Notes on BA26 that need checking against the typeset PDF

These come from the publisher HTML, where equations are flattened into
plain text:

- Eq. (4) appears without the minus sign that follows from eq. (3)
  ($\delta X=-\hat J^{-1}F_{N}\delta N$). This may be a transcription
  artefact.
- The backtracking rule is worded for minimization, but the objective is
  maximized.
- Both runs report $I_{ON}$ rising from 14.7 to 15.1 μA/μm. That is about
  2.7%, not the "about 3%" and "about 14%" stated elsewhere. The
  BV percentages differ between §IV and the conclusion.

None of these change the method. They are recorded so that any future
comparison uses the numbers in the results section, not the summaries.

### 2.4 Reported device results (comparison targets only)

- **ZA19.** Si PT-IGBT, $L=10$ μm, $W=70$ μm, drift doping
  $5\times10^{13}$ cm⁻³, 75 nm oxide. BV rises from 1060 V to about 1140 V
  in 5 iterations. Constrained runs reach about 1150 V. $I_{ON}$ at
  $V_{CE}=1.5$ V rises by more than 50 A/cm².
- **BA26.** 4H–SiC 1.2 kV vertical MOSFET half-cell, 9 μm drift at
  $10^{16}$ cm⁻³, 40 nm oxide. BV 1194 V, $V_{th}=2.7$ V and
  $I_{ON}=14.7$ μA/μm at baseline. The full-device run reaches BV 1347 V,
  $I_{ON}$ 15.1 μA/μm and $V_{th}$ 0.7 V. The drift-only run reaches BV
  1322 V with $V_{th}$ unchanged.

## 3. Open-source solvers surveyed

| Project | Scope | Licence (GitHub API) | Use here |
|---|---|---|---|
| DEVSIM (`devsim/devsim`, S22) | General TCAD, finite volume, Scharfetter–Gummel, Python scripting | Apache-2.0 | Candidate **independent cross-check** for drift–diffusion stages. A licence-compatible dependency if needed |
| ChargeTransport.jl (`WIAS-PDELib/ChargeTransport.jl`) | Drift–diffusion, Voronoi finite volume, Julia | AGPL-3.0 | Methodological reference only. **No code is copied** (copyleft) |
| SIMsalabim (`kostergroup/SIMsalabim`) | 1D drift–diffusion for LEDs and solar cells | GPL-3.0 | Methodological reference only. **No code is copied** |
| Charon (Sandia) | Production TCAD | **unverified**: no public GitHub repository found | Not used |

Commercial simulators (Sentaurus, Atlas) were used by ZA19 and BA26 for
validation and parameters. They are not available here, and their default
parameter sets cannot be cited as primary sources.

## 4. Reproduction feasibility

| Target | Status | Missing |
|---|---|---|
| ZA19 IGBT BV/$I_{ON}$ optimization | **Cannot reproduce** | Supplementary material (models, parameters, geometry script), mesh, impact-ionization coefficients |
| BA26 SiC MOSFET optimization | **Cannot reproduce** | All numerical model parameters (Atlas defaults not printed), full cell width, contact BCs, temperature, mesh, impact-ionization coefficients |
| Adjoint sensitivity method (both) | **Can be independently implemented** | Nothing. The method is fully specified at the discrete level |

The project therefore labels its work as follows:

- *Reproduction*: none claimed.
- *Independent implementation inspired by published work*: the adjoint
  doping-sensitivity method (ZA19, BA26) and the field-based BV index
  (ZA19 eq. 1).
- *New contribution*: none yet. PINN and hybrid solvers for these equations
  remain **unverified hypotheses** until benchmarked.

## 5. Requirements derived

Each stage has to pass its own verification gate before the next one starts.

1. **Linear Poisson with heterogeneous permittivity.** MMS, FV reference and
   PINN. This is benchmark P (`docs/benchmarks/poisson_1d.md`).
2. **Nonlinear equilibrium Poisson** with Boltzmann carriers and a 1D abrupt
   pn junction. Verify with Newton quadratic convergence, comparison with the
   depletion approximation in its regime of validity, and charge neutrality
   in the bulk.
3. **Drift–diffusion** with Scharfetter–Gummel fluxes (SG69). Verify with
   MMS for the coupled system, an ideal-diode current check, and a
   cross-check against DEVSIM on the same 1D structure.
4. **Field-based BV index** (ZA19 eq. 1) in 1D. Verify a punch-through diode
   against the analytical field profile.
5. **Adjoint sensitivities.** Verify that the adjoint gradient matches
   central finite differences of the discrete objective at randomly sampled
   doping nodes (relative agreement to be pre-registered), and that the cost
   is about one linear solve.
6. **4H–SiC material models.** Every parameter must be traced to a primary
   source cited where it is used: permittivity, band gap, intrinsic density,
   incomplete-ionization energies (N, Al), anisotropic mobility,
   impact-ionization coefficients. **No numerical SiC parameter is
   introduced until its primary source is read.** BA26's Atlas defaults do
   not qualify.
7. **PINN and SciML stages** (surrogates, inverse doping problems). Each is
   compared against the verified classical solver on the same discrete
   problem, with all seeds reported. Known failure modes in steep layers
   (WTP21) are tested explicitly, as in benchmark case P4.

## 6. Unresolved physics and open questions

- 4H–SiC is anisotropic (hexagonal). Whether 2D stages need tensor mobility
  and permittivity depends on device orientation. BA26 does not say.
- Incomplete ionization of Al acceptors in 4H–SiC is strongly
  temperature-dependent, and BA26 does not state the temperature.
- A field-based BV index with a constant $E_{crit}$ is a proxy for avalanche
  breakdown. Its agreement with impact-ionization integrals for SiC has to be
  checked by this project, not assumed (ZA21 has not been read).
- The SiC/SiO₂ interface (traps, fixed charge, surface mobility) dominates
  $V_{th}$ and channel resistance. It is out of scope until stage 6.
