"""Minimal example: solve benchmark case P1 with the reference solver and a short PINN.

Run:  python examples/minimal_poisson.py
The PINN here uses a shortened schedule (about 5 s on a laptop CPU) for
illustration. Benchmark numbers come from benchmarks/poisson1d/run.py.
"""

import numpy as np

from widegap.fv.poisson1d import solve, uniform_mesh
from widegap.problems.poisson1d import evaluate, p1_smooth
from widegap.verification import max_error, rel_l2_error

problem = p1_smooth()

fv = solve(problem, uniform_mesh(256))
print(f"FV  (N=256):      max error = {max_error(fv.u, evaluate(problem.u_exact, fv.x)):.3e}")

try:
    from widegap.pinn.poisson1d import PINNConfig, train
except ImportError:
    print(
        "PINN skipped: install torch (pip install torch --index-url https://download.pytorch.org/whl/cpu)"
    )
else:
    res = train(problem, PINNConfig(adam_steps=500, lbfgs_max_iter=100, seed=0))
    x = np.linspace(0.0, 1.0, 2001)
    err = rel_l2_error(x, res.predict(problem, x), evaluate(problem.u_exact, x))
    print(f"PINN (short run): rel L2 error = {err:.3e}  ({res.train_seconds:.1f} s)")
