"""Generate figures for the 1-D Poisson benchmark from the raw result JSON."""

import argparse
import json
from pathlib import Path

import matplotlib
import matplotlib.ticker

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

# Okabe-Ito palette
BLUE, ORANGE, GREEN, VERMILION = "#0072B2", "#E69F00", "#009E73", "#D55E00"
CASE_COLOR = {"P1": BLUE, "P2": ORANGE, "P4": GREEN}
CASE_MARKER = {"P1": "o", "P2": "s", "P4": "^"}
CASES = ("P1", "P2", "P4")
SEEDS = range(5)

plt.rcParams.update(
    {
        "font.size": 10,
        "axes.labelsize": 11,
        "axes.grid": True,
        "grid.alpha": 0.3,
        "legend.frameon": False,
    }
)


def load(path):
    return json.loads(Path(path).read_text())


def footer(fig, manifest):
    commit = str(manifest.get("git_commit", "unknown"))[:7]
    cpu = manifest.get("cpu_model", "unknown CPU")
    fig.text(
        0.5,
        0.01,
        f"widegap | commit {commit} | {cpu} | CPU only",
        ha="center",
        va="bottom",
        fontsize=7,
        color="0.35",
    )


def save(fig, out, name):
    for ext in ("png", "svg"):
        fig.savefig(out / f"{name}.{ext}", dpi=200)
    plt.close(fig)


def fig_fv(fv, manifest, out):
    fig, ax = plt.subplots(figsize=(5.5, 4.3))
    h = 1.0 / np.asarray(fv["levels"], dtype=float)
    for case in CASES:
        ax.loglog(
            h,
            fv["cases"][case]["e_inf"],
            marker=CASE_MARKER[case],
            color=CASE_COLOR[case],
            label=case,
        )
    e0 = fv["cases"]["P1"]["e_inf"][0]
    ax.loglog(h, e0 * (h / h[0]) ** 2, "k--", lw=1, label=r"$O(h^2)$ reference")
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xticks(h, [f"1/{n}" for n in fv["levels"]])
    ax.set_xlabel("mesh spacing h")
    ax.set_ylabel("max nodal error")
    ax.set_title("Finite-volume convergence, 1-D Poisson")
    ax.legend()
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    footer(fig, manifest)
    save(fig, out, "poisson1d_fv_convergence")


def fig_errors(results, manifest, out):
    fig, ax = plt.subplots(figsize=(5.5, 4.3))
    offsets = np.linspace(-0.15, 0.15, len(SEEDS))
    for i, case in enumerate(CASES):
        for seed in SEEDS:
            p = load(results / f"pinn_{case}_seed{seed}.json")
            x = i + offsets[seed]
            y = p["e_l2_rel_exact"]
            ax.plot(x, y, CASE_MARKER[case], color=CASE_COLOR[case], ms=6)
            ax.annotate(
                str(seed),
                (x, y),
                textcoords="offset points",
                xytext=(0, 7),
                ha="center",
                fontsize=7,
                color="0.3",
            )
    ax.axhline(1e-3, ls="--", color=VERMILION, lw=1)
    ax.text(
        2.45,
        1.15e-3,
        "N1 threshold (P1 only)",
        color=VERMILION,
        ha="right",
        va="bottom",
        fontsize=8,
    )
    ax.set_yscale("log")
    ax.set_xticks(range(len(CASES)), CASES)
    ax.set_xlim(-0.5, 2.5)
    ax.margins(y=0.25)
    ax.set_xlabel("case (numbers are seeds)")
    ax.set_ylabel("relative $L^2$ error vs exact")
    ax.set_title("PINN accuracy by seed")
    ax.grid(axis="x", visible=False)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    footer(fig, manifest)
    save(fig, out, "poisson1d_pinn_errors")


def fig_loss(results, manifest, out):
    fig, ax = plt.subplots(figsize=(5.5, 4.3))
    colors = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7"]
    for seed in SEEDS:
        p = load(results / f"pinn_P1_seed{seed}.json")
        hist = np.asarray(p["loss_history_every_100"], dtype=float)
        ax.plot(np.arange(len(hist)), hist, color=colors[seed], lw=1.3, label=f"seed {seed}")
    ax.set_yscale("log")
    ax.set_xlabel("loss evaluation index (every 100th)")
    ax.set_ylabel("training loss")
    ax.set_title("P1 PINN training loss (Adam then L-BFGS)")
    ax.legend()
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    footer(fig, manifest)
    save(fig, out, "poisson1d_pinn_loss")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", default="results/poisson1d")
    ap.add_argument("--out", default="figures/")
    args = ap.parse_args()
    results, out = Path(args.results), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    fv = load(results / "fv.json")
    manifest = load(results / "manifest.json")
    fig_fv(fv, manifest, out)
    fig_errors(results, manifest, out)
    fig_loss(results, manifest, out)


if __name__ == "__main__":
    main()
