# GitHub repository metadata (to be applied by the owner)

These settings cannot be applied from the development machine (no `gh`
CLI). Apply them in the GitHub UI or with `gh repo edit` after the
repository is created.

**Description**

> Verified FV and physics-informed (PINN) solvers for semiconductor device equations, built toward adjoint doping optimization of 4H-SiC power devices. Benchmarks pre-registered; every seed reported.

**Topics**

`semiconductor-simulation` `tcad` `wide-bandgap` `silicon-carbide` `scientific-machine-learning` `physics-informed-neural-networks` `pinn` `numerical-pde` `finite-volume` `method-of-manufactured-solutions`

**Suggested issues to open** (known limitations, see README)

1. Licence not yet selected (blocks first release).
2. Only linear 1D Poisson is implemented; nonlinear Poisson with Boltzmann carriers is the next stage.
3. PINN timings are CPU-only (no GPU on the development machine).
4. 4H-SiC material parameters: primary sources not yet collected.
5. CI workflow configured but its first run on GitHub Actions has not yet been observed.
