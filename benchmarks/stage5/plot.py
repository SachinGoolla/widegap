"""Figure: adjoint doping-sensitivity maps from results/stage5/results.json.

Run:  python benchmarks/stage5/plot.py
Plots gamma_D = -gamma_C (Zhu & Andrei sign convention, D = N_A - N_D) per unit
doping for the diode current (A2) and the field-based BV (A4 is recomputed here).
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from widegap import adjoint as A  # noqa: E402
from widegap.fv import reverse_bias as rb  # noqa: E402
from widegap.fv.equilibrium import control_volumes, cv_average_piecewise  # noqa: E402

REPO = Path(__file__).resolve().parents[2]


def main():
    res = json.loads((REPO / "results/stage5/results.json").read_text())
    man = json.loads((REPO / "results/stage5/manifest.json").read_text())
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))

    x2 = np.linspace(-40.0, 40.0, 1025)
    g2 = np.asarray(res["A2"]["gamma_C"]) / control_volumes(x2)  # per unit length
    ax[0].plot(x2, -g2 / res["A2"]["J"], color="#0072B2")
    ax[0].set_xlabel("x  [L_D]")
    ax[0].set_ylabel(r"$\gamma_D / \bar j$  per unit length")
    ax[0].set_title("Diode current, V = 6 V_T (A2)")
    ax[0].axvline(0, color="0.6", lw=0.8)

    x4 = np.linspace(-5.0, 60.0, 2601)
    c4 = cv_average_piecewise(x4, [0.0], [-100.0, 1.0])
    d = rb.Diode1D(x4, c4, 1e-3)
    sw = rb.sweep(d, range(0, 121))
    bv, psi, _ = rb.breakdown_voltage(d, sw, 15.0)
    g4, _ = A.breakdown_voltage_sensitivity(d, A.polish_reverse(d, psi, bv), bv, 15.0)
    ax[1].plot(x4, -g4 / control_volumes(x4), color="#D55E00")
    ax[1].set_xlim(-2, 20)
    ax[1].set_xlabel("x  [L_D]")
    ax[1].set_ylabel(r"$\gamma_D^{BV}$  [V_T per unit doping-length]")
    ax[1].set_title(f"Field-based BV = {bv:.2f} V_T (A4, proxy)")
    ax[1].axvline(0, color="0.6", lw=0.8)

    fig.text(
        0.5,
        0.005,
        f"widegap | commit {man['git_commit'][:7]} | {man['cpu_model']} | CPU only",
        ha="center",
        fontsize=7,
        color="0.4",
    )
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    out = REPO / "figures"
    out.mkdir(exist_ok=True)
    for ext in ("png", "svg"):
        fig.savefig(out / f"stage5_sensitivity_maps.{ext}", dpi=200)


if __name__ == "__main__":
    main()
