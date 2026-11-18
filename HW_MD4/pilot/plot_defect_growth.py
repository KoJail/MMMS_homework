"""Compare defect growth on absolute and onset-aligned strain axes."""
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

cases = ("flat-01", "flat-02", "step-01", "step-02")
colors = {"flat": "#245474", "step": "#b05337"}
fig, axes = plt.subplots(2, 2, figsize=(9.2, 6.1), layout="constrained")
for case in cases:
    group, rep = case.split("-")
    with (Path("pilot/audit") / case / "loading_audit.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    strain = np.array([float(r["v_strain"]) * 100 for r in rows])
    length = np.array([float(r["total_line_angstrom"]) / 10 for r in rows])
    hcp = np.array([int(r["hcp"]) / int(r["atom_count"]) * 100 for r in rows])
    onset = strain[np.flatnonzero(length)[0]]
    style = dict(color=colors[group], ls="-" if rep == "01" else "--",
                 lw=1.65, label=case)
    axes[0, 0].plot(strain, length, **style)
    axes[0, 1].plot(strain - onset, length, **style)
    axes[1, 0].plot(strain, hcp, **style)
    axes[1, 1].plot(strain - onset, hcp, **style)
for row, ylabel in zip(axes, ("Total DXA line length (nm)", "PTM HCP atoms (%)")):
    row[0].set_ylabel(ylabel)
    for ax in row:
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", color="#dce4e9", linewidth=.65)
        ax.set_axisbelow(True)
for ax in axes[:, 0]:
    ax.set_xlim(4.5, 12)
for ax in axes[:, 1]:
    ax.set_xlim(-.15, 1.2)
    ax.axvline(0, color="#85919a", lw=.8, ls=":")
axes[0, 0].set(title="(a) Defect population", xticklabels=[])
axes[0, 1].set(title="(b) Early line development", xticklabels=[])
axes[1, 0].set(title="(c) Local HCP fraction", xlabel="Engineering strain (%)")
axes[1, 1].set(title="(d) Early HCP development", xlabel="Strain since first DXA (percentage points)")
axes[0, 0].legend(ncol=2, frameon=False, loc="upper left", fontsize=8)
fig.savefig("结课报告/图/05-缺陷形成与扩展.png", dpi=270)
