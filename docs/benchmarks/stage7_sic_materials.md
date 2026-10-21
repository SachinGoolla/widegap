# Stage 7: 4H-SiC material models with enforced provenance

Status: **pre-registered** (2026-10-21). No stage-7 code exists and no SiC
parameter has been collected at the time of this commit.

## 1. Provenance rule (enforced in code)

Every numerical 4H-SiC parameter lives in one registry file. Each entry has
these fields:

- value and unit;
- temperature and validity range, *as stated in the source*;
- crystal orientation, where relevant;
- source DOI, verified against Crossref (title, authors, year);
- a locator (table, equation, figure or page);
- the access status of the full text;
- the person or agent who extracted it.

Entries record locators and paraphrase. No copyrighted text is stored
beyond a few words needed to identify the value.

The loader **refuses** any parameter whose full text was not read
(`access != "full-text-read"`). It also refuses any parameter not
independently re-verified by a second reader who did not extract it
(`verified_by` present and different from `extracted_by`). Parameter
classes:

- `measured`: a primary experimental or first-principles source;
- `convention`: a theory constant, such as a degeneracy factor. It must be
  cited to a read source that states and uses it, and it is flagged as such.

A parameter that cannot be sourced this way is **excluded**, and any model
that needs it is disabled. It is never substituted.

**Required set R** (the stage-7 gate needs every member):

- relative permittivity, at least one orientation;
- $E_g$ at 300 K;
- $N_c$ and $N_v$ at 300 K, or the density-of-states effective masses;
- ionization energy of the nitrogen donor (at least one lattice site);
- ionization energy of the aluminium acceptor;
- donor and acceptor degeneracy factors.

**Optional set O** (reported): both permittivity orientations,
temperature dependences, and Caughey–Thomas mobility parameters.

## 2. Models

- $n_i=\sqrt{N_cN_v}\,e^{-E_g/2k_BT}$.
- Incomplete ionization:
  $N_D^+=N_D/(1+g_D\,(n/N_c)\,e^{\Delta E_D/k_BT})$ and
  $N_A^-=N_A/(1+g_A\,(p/N_v)\,e^{\Delta E_A/k_BT})$, with Boltzmann
  $n$ and $p$. Optionally there are two donor sites (hexagonal and cubic)
  with stated fractions.
- Equilibrium Poisson (stage-2 operator, scaled):
  $-\psi''=\delta(e^{-\psi}-e^{\psi})+\tilde N_D^+(\psi)-\tilde N_A^-(\psi)$.
  Newton uses the analytic Jacobian and the stage-3b round-off-relative
  tolerance.
- Caughey–Thomas mobility $\mu(N)$, if the set O parameters are available.

## 3. Acceptance criteria

| ID | Criterion | Gating |
|----|-----------|--------|
| M1 | Provenance: every registry entry passes the loader rules of §1 (unit-tested). Set R is complete. Each DOI matches Crossref metadata. Each value is confirmed by the independent re-verifier against the read text | yes |
| M2 | Neutral-bulk solve with incomplete ionization (one donor level) matches the closed-form quadratic $n^2+n\,n_1-N_Dn_1=0$, $n_1=(N_c/g_D)e^{-\Delta E_D/k_BT}$, to ≤ 1e-10 relative, over $N_D\in[10^{14},10^{19}]$ cm⁻³ and $T\in\{200,300,500\}$ K | yes |
| M3 | Limits: the ionized fraction is above $1-10^{-6}$ when $N_D\le10^{-6}n_1$, and below $10^{-3}$ in deep freeze-out ($n_1\le10^{-6}N_D$) | yes |
| M4 | Equilibrium solver with incomplete ionization, MMS on $[0,1]$: observed max-norm order in [1.9, 2.1] (last three of N = 64…1024) | yes |
| M5 | 4H-SiC abrupt pn junction at 300 K ($N_A=10^{17}$, $N_D=10^{16}$ cm⁻³, realistic $n_i/N_0$ about $10^{-25}$). (a) Newton converges, with the stage-3b tolerance. (b) The bulk potentials equal independent neutrality root-solves to ≤ 1e-10 absolute in $V_T$. (c) The junction field, Richardson-extrapolated over two meshes, matches the exact first integral with incomplete ionization (closed form: $\int N_D^+d\psi=N_D[\psi-\log(e^{\psi}+\gamma)]$-type antiderivatives) to ≤ 2e-3 relative | yes |
| M6 | Caughey–Thomas implementation identities: $\mu(N\to0)=\mu_{max}$, $\mu(N_{ref})=(\mu_{max}+\mu_{min})/2$ to 1e-14 | yes, if set O mobility exists; otherwise reported as unavailable |
| R7a | $n_i(300\,\text{K})$ from the registry compared with any $n_i$ value the sources report | reported |
| R7b | Stage-3 drift-diffusion with realistic $n_i/N_0\sim10^{-25}$: expected **unsupported** in $(\psi,n,p)$ variables, because minority currents of size $\delta^2\sim10^{-49}$ are below the round-off of the majority terms. Attempted and the outcome documented | reported |

## 4. Gate and claims

Stage 7 passes if M1–M5 hold, plus M6 where applicable. A pass means the
**implemented models are correct and every number is traceable**. It does
not mean the models are adequate for any particular SiC device. Interface
traps, anisotropic mobility tensors, impact ionization and high-field
transport remain out of scope.

## 5. Amendments

### 2026-10-21: DOI verification is not limited to Crossref (before the official run)

§1 required DOIs to be "verified against Crossref". One read source, R.
Ishikawa's PhD thesis (Kyoto University 2025, DOI 10.14989/doctor.k26017),
has a JaLC-registered DOI that Crossref does not index. Restricting to one
registration agency was arbitrary. The requirement becomes: the DOI's
metadata from its registration agency must match the title, first author
and year. It may be retrieved from Crossref, or through doi.org content
negotiation (`Accept: application/vnd.citationstyles.csl+json`). The
registry field is renamed `doi_metadata_match` and records the agency used.

**Disclosure:** at the time of this amendment the literature search had
already shown that set R is incomplete (no readable primary source for
$E_g$(300 K), $N_c$, $N_v$, or site-resolved N-donor energies). The
amendment therefore cannot change the stage-7 gate outcome. It only decides
whether the thesis-sourced entries ($g_D$ and the Hall-mobility fits) are
usable.
