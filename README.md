# wbg-sciml-tcad

Verified finite-volume and physics-informed neural network (PINN) solvers for semiconductor device equations, built in stages toward adjoint doping optimization of wide-bandgap (4H-SiC) power devices. **The current release solves only the linear 1D Poisson verification benchmark; there is no device simulation yet.**

## Status

Each stage must pass its own verification gate before the next starts (`docs/literature/requirements_report.md`, section 5).

| Stage | Content | Status |
|---|---|---|
| 1 | Linear Poisson with heterogeneous permittivity (MMS, FV reference, PINN) | verified (gate passed) |
| 2 | Nonlinear equilibrium Poisson, Boltzmann carriers, 1D abrupt pn junction | not started |
| 3 | Drift-diffusion with Scharfetter-Gummel fluxes, cross-check against DEVSIM | not started |
| 4 | Field-based breakdown-voltage index in 1D | not started |
| 5 | Adjoint sensitivities checked against finite differences | not started |
| 6 | 4H-SiC material models (each parameter traced to a primary source) | not started |
| 7 | PINN and SciML stages (surrogates, inverse doping problems) | not started |

## Installation

```bash
python -m venv .venv
source .venv/bin/activate
pip install torch --index-url https://download.pytorch.org/whl/cpu   # optional, CPU-only
pip install -e ".[test]"
```

Python >= 3.10. The package is named `wbgtcad`. torch is optional: the finite-volume reference solver depends only on NumPy and SciPy, and the PINN requires torch (CI runs the test suite and example both with and without torch).

## Minimal example

```bash
python examples/minimal_poisson.py
```

It solves benchmark case P1 with the finite-volume solver and, if torch is installed, a shortened PINN run. It prints the max error of the FV solution and the relative L2 error of the short PINN run (with its training time). Without torch it prints the FV line and a message that the PINN was skipped. These numbers are not the benchmark results below.

## Benchmark results (stage 1)

Problem definitions, error measures and acceptance criteria were fixed before any solver result existed: [docs/benchmarks/poisson_1d.md](docs/benchmarks/poisson_1d.md). One dated amendment (R5 flux balance, restated as R5b/R5c) was also made before any solver output. All quantities are dimensionless.

### Finite-volume reference solver (uniform meshes)

Max-norm nodal error against the exact solution (`results/poisson1d/fv.json`):

| Case | N=32 | N=64 | N=128 | N=256 | N=512 | Observed orders (N=32→64 … 256→512) |
|---|---|---|---|---|---|---|
| P1 | 8.0e-04 | 2.0e-04 | 5.0e-05 | 1.3e-05 | 3.1e-06 | 2.00, 2.00, 2.00, 2.00 |
| P2 | 6.5e-04 | 1.6e-04 | 4.0e-05 | 1.0e-05 | 2.5e-06 | 2.00, 2.00, 2.00, 2.00 |
| P4 | 1.7e-01 | 3.9e-02 | 9.7e-03 | 2.4e-03 | 6.1e-04 | 2.13, 2.01, 1.99, 2.00 |

- P1 and P2 observed orders are within the pre-registered window [1.9, 2.1]. P4 is reported only (not gating).
- P3 (piecewise-constant permittivity, harmonic-mean element permittivity): max error 5.0e-15 for N=64 (interface on a node) and 5.3e-14 for N=63 (interface inside an element).
- Linear exact solution: max error 5.7e-15.
- R5b (discrete flux balance, P1, N=512): relative imbalance 2.1e-13.
- R5c (reported, not gating): error of the end-element boundary flux for N=32 to 512 is 2.5e-03, 6.3e-04, 1.6e-04, 3.9e-05, 9.9e-06; observed orders 2.00, 2.00, 2.00, 2.00.

### PINN (3x32 tanh, float64, hard Dirichlet conditions, Adam then L-BFGS)

All five seeds are reported; none was discarded (`results/poisson1d/pinn_<case>_seed<k>.json`). Relative L2 errors use a fixed 2001-point grid.

| Case | Seed | rel. L2 vs exact | rel. L2 vs FV (N=512) | Train time (s) |
|---|---|---|---|---|
| P1 | 0 | 3.2e-07 | 1.6e-06 | 26.7 |
| P1 | 1 | 2.0e-06 | 3.6e-06 | 23.5 |
| P1 | 2 | 5.5e-06 | 7.3e-06 | 27.0 |
| P1 | 3 | 2.5e-06 | 4.1e-06 | 25.9 |
| P1 | 4 | 4.0e-07 | 2.1e-06 | 25.5 |
| P2 | 0 | 1.2e-06 | 2.2e-06 | 25.5 |
| P2 | 1 | 3.1e-07 | 1.7e-06 | 26.4 |
| P2 | 2 | 1.1e-06 | 2.4e-06 | 25.4 |
| P2 | 3 | 1.3e-06 | 2.5e-06 | 27.0 |
| P2 | 4 | 1.0e-06 | 2.1e-06 | 26.8 |
| P4 | 0 | 3.0e-04 | 3.2e-04 | 25.0 |
| P4 | 1 | 9.7e-06 | 1.2e-04 | 26.4 |
| P4 | 2 | 4.2e-04 | 4.3e-04 | 24.3 |
| P4 | 3 | 2.8e-05 | 1.2e-04 | 28.3 |
| P4 | 4 | 1.1e-04 | 1.6e-04 | 24.1 |

- N2: the "vs FV" column is dominated by the FV solution's own error (max 3.1e-06 for P1 at N=512), so it is larger than the "vs exact" column; the two comparisons agree to within that level, as pre-registered.
- P4 (not gating): the rel. L2 error vs exact varies by about 40x across seeds (smallest 9.7e-06, largest 4.2e-04).

### Gate

From `results/poisson1d/summary.json`: R1 (P1 and P2), R2, R3, R4, R5 (linear), R5b and N1 all hold; `gate_passed` is true and `pre_registered_configuration` is true.

### Provenance

- Hardware (`manifest.json`): AMD Ryzen 7 5700U with Radeon Graphics, 16 logical CPUs, 8 torch threads, CPU-only torch (2.14.1+cpu), no GPU. Python 3.12.3, NumPy 2.5.4, SciPy 1.18.1.
- Peak resident memory per PINN run: about 370 MiB (`peak_rss_mib` in the `pinn_*.json` files).
- Results were produced by git commit `91e0ee91c2636c698e6b82d48c92664b4e61f437` (`git_dirty` false in the manifest); an earlier run of the identical training code (commit `3d1990f`, which crashed only when writing `summary.json`; fixed in `91e0ee9`) produced bit-identical errors for all 15 PINN runs.
- Reproduce:

```bash
python benchmarks/poisson1d/run.py --out results/poisson1d
```

PINN timings are wall-clock on this CPU and are not an acceptance criterion.

## Supported features

- 1D vertex-centred (box-method) finite-volume solver for $-(\varepsilon u')' = f$ with Dirichlet conditions, length-weighted harmonic-mean element permittivity, banded direct solve.
- 1D PINN for the same equation: strong-form residual loss via automatic differentiation, Dirichlet conditions imposed exactly by construction.
- Manufactured-solution problems P1 to P4 and error measures (max error, relative L2, observed order).
- Reproducible benchmark runner that writes per-case and per-seed JSON and a manifest.

## Limitations

- 1D only; linear Poisson only. No carriers, doping, drift-diffusion, or breakdown/adjoint machinery.
- No 4H-SiC parameters and no physical units; permittivity contrast and all values are dimensionless.
- PINN timings are CPU-only and were measured on one machine.
- The PINN error on the steep layer P4 varies by about 40x across seeds.
- CI is configured (`.github/workflows/ci.yml`) but its first GitHub Actions run has not yet been observed.
- No licence has been chosen.

## Scientific attribution

The method context comes from the following works. **No published device result is reproduced.**

- C. Zhu, P. Andrei, AIP Advances 9, 095301 (2019), [doi:10.1063/1.5113764](https://doi.org/10.1063/1.5113764): adjoint doping-sensitivity method and field-based breakdown index (silicon IGBT).
- P. Brungi, P. Andrei, AIP Advances 16, 035014 (2026), [doi:10.1063/5.0320204](https://doi.org/10.1063/5.0320204): 4H-SiC MOSFET context and comparison targets.
- M. Raissi, P. Perdikaris, G. E. Karniadakis, J. Comput. Phys. 378, 686-707 (2019), [doi:10.1016/j.jcp.2018.10.045](https://doi.org/10.1016/j.jcp.2018.10.045): PINNs.
- P. J. Roache, J. Fluids Eng. 124, 4-10 (2002), [doi:10.1115/1.1436090](https://doi.org/10.1115/1.1436090): method of manufactured solutions.

Classification of this project's work (`docs/literature/requirements_report.md`, section 4):

- *Reproduction*: none claimed. Both device-level optimization studies cannot be reproduced from the published information (missing supplementary material, model parameters, meshes and impact-ionization coefficients).
- *Independent implementation inspired by published work*: the adjoint doping-sensitivity method and the field-based BV index are planned (stages 4 and 5); neither is implemented yet.
- *New contribution*: none yet. PINN and hybrid solvers for these equations remain unverified hypotheses until benchmarked.

This project is not affiliated with and has not been endorsed by the cited authors.

## Tests

```bash
pytest
```

- R1: P1 and P2 observed max-norm orders in [1.9, 2.1].
- R2, R3: max-error thresholds for P1 and P2 at N=512.
- R4: P3 reproduced to round-off with the interface on a node and inside an element.
- R5: a linear solution is reproduced to round-off.
- R5b: discrete flux balance holds to round-off (P1, N=512).
- N1: PINN on P1 reaches rel. L2 error of at most 1e-3 against the exact solution for all five seeds.

The R1 to R5b checks run on the reference solver; the full N1 check is part of the benchmark runner.

## Documentation

- [Benchmark P1-P4 pre-registration and amendment](docs/benchmarks/poisson_1d.md)
- [Literature and requirements report](docs/literature/requirements_report.md)
- [Raw results and manifest](results/poisson1d/)

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md).

## Citation

See [CITATION.cff](CITATION.cff). Please also cite the original works for the methods used.

## License

Not yet selected (see [docs/repository_metadata.md](docs/repository_metadata.md)). Until a licence file is added, no licence is granted.
