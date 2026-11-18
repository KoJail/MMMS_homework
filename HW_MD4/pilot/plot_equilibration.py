"""Compare thermal equilibration and residual in-plane stresses across cases."""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def equilibrium_series(path):
    rows = []
    active = False
    for line in path.read_text().splitlines():
        fields = line.split()
        if fields[:4] == ["Step", "Atoms", "Temp", "PotEng"]:
            active = True
            continue
        if active:
            if len(fields) != 11:
                active = False
                continue
            try:
                rows.append([float(x) for x in fields])
            except ValueError:
                active = False
    series = np.array(rows)
    if not len(series):
        raise ValueError(f"No equilibrium data in {path}")
    # The minimizer precedes the NVT run, which resets step to zero.
    resets = np.flatnonzero(np.diff(series[:, 0]) < 0)
    return series[resets[-1] + 1:] if len(resets) else series


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cases", nargs="+", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(3, 1, sharex=True, figsize=(9, 7), layout="constrained")
    for case in args.cases:
        a = equilibrium_series(case / "equil.log")
        time = a[:, 0] * .001
        stress = -a[:, 8:10] * a[:, 7, None] / (18 * 3.615) * 1e-4
        color = "#b13f25" if case.name.startswith("step") else "#176ba0"
        style = "--" if "02" in case.name else "-"
        axes[0].plot(time, a[:, 2], color=color, ls=style, label=case.name, lw=.85)
        axes[1].plot(time, a[:, 3] / a[:, 1], color=color, ls=style, lw=.85)
        axes[2].plot(time, stress[:, 0], color=color, ls=style, lw=.85)
        print(case.name, "final_ps", time[-1])
        for i in (0, 1):
            mask = (time >= time[-1] - (2-i)*5) & (time <= time[-1] - (1-i)*5)
            if mask.any():
                print("  block", i + 1, "T_K", round(a[mask, 2].mean(), 2),
                      "mean_in_plane_GPa", np.round(stress[mask].mean(axis=0), 3))
    axes[0].axhline(300, color="black", lw=.7, alpha=.5)
    axes[0].set_ylabel("Temperature (K)")
    axes[0].legend(frameon=False)
    axes[1].set_ylabel("Potential energy (eV/atom)")
    axes[2].set_ylabel("Nominal axial stress (GPa)")
    axes[2].set_xlabel("NVT equilibration time (ps)")
    for ax in axes:
        ax.grid(alpha=.2)
    fig.savefig(args.out, dpi=180)


if __name__ == "__main__":
    main()
