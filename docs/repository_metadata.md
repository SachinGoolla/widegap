# GitHub repository metadata (to be applied by the owner)

Applied with `gh repo create` / `gh repo edit` when the repository is
created; kept here as the source of truth.

**Description**

> Verified FV and physics-informed (PINN) solvers for semiconductor device equations, built toward adjoint doping optimization of 4H-SiC power devices. Benchmarks pre-registered; every seed reported.

**Topics**

`semiconductor-simulation` `tcad` `wide-bandgap` `silicon-carbide` `scientific-machine-learning` `physics-informed-neural-networks` `pinn` `numerical-pde` `finite-volume` `method-of-manufactured-solutions`

**Suggested issues to open** (known limitations, see README)

1. Only linear 1D Poisson is implemented; nonlinear Poisson with Boltzmann carriers is the next stage.
2. PINN timings are CPU-only (no GPU on the development machine).
3. 4H-SiC material parameters: primary sources not yet collected.
4. CI workflow configured but its first run on GitHub Actions has not yet been observed.
