"""Make the initial geometry and pilot load audit figures."""
import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from ovito.io import import_file

from analyze import load_curve


def plot_initial(cases, output):
    fig, panels = plt.subplots(2, len(cases), figsize=(3.6 * len(cases), 5.5),
                               squeeze=False, layout="constrained")
    for col, case in enumerate(cases):
        for row, stage in enumerate(("initial", "equilibrated")):
            frame = import_file(str(case / f"{stage}.dump")).compute()
            x, y, z = np.asarray(frame.particles.positions).T
            # Keep a thin center slice so projected rows remain legible.
            keep = np.abs(y - frame.cell[1, 1] / 2) < 3.0
            ax = panels[row, col]
            ax.scatter(x[keep] / 10, z[keep] / 10, s=3, alpha=0.48,
                       color="#035e76", linewidths=0)
            ax.set(xlim=(-.05, 5.2), ylim=(-.25, 6.8), aspect="equal",
                   xlabel="Loading axis [110] (nm)",
                   title=f"{case.name} · {stage}")
            if col == 0:
                ax.set_ylabel("Surface normal [001] (nm)")
    fig.savefig(output / "pilot-model-cross-sections.png", dpi=200)
    plt.close(fig)


def plot_loading(cases, audit, output):
    fig, axes = plt.subplots(2, 1, figsize=(9, 6.4), sharex=True, layout="constrained")
    for case in cases:
        curve = load_curve(case / "loading.log")
        rows = list(curve.values())
        strain = np.array([row["v_strain"] for row in rows]) * 100
        stress = [row["v_axial"] for row in rows]
        color = "#b13f25" if case.name.startswith("step") else "#176ba0"
        linestyle = "--" if case.name.endswith("02") else "-"
        axes[0].plot(strain, stress, label=case.name, color=color, ls=linestyle)
        with (audit / case.name / "loading_audit.csv").open() as stream:
            frames = list(csv.DictReader(stream))
        axes[1].plot([100 * float(row["v_strain"]) for row in frames],
                     [float(row["total_line_angstrom"]) / 10 for row in frames],
                     color=color, ls=linestyle)
    axes[0].set_ylabel("Nominal axial stress (GPa)")
    axes[0].legend(frameon=False, ncols=4)
    axes[1].set(xlabel="Engineering strain (%)", ylabel="DXA line length (nm)")
    for ax in axes:
        ax.grid(alpha=.2)
    fig.savefig(output / "pilot-load-dxa.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--audit", type=Path)
    parser.add_argument("cases", nargs="+", type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    plot_initial(args.cases, args.out)
    if args.audit:
        plot_loading(args.cases, args.audit, args.out)
