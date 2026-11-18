"""Plot an auditable sequence of candidate and persistent DXA signals."""
import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

from analyze import load_curve, pipeline_for


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("case", type=Path)
    p.add_argument("audit", type=Path)
    p.add_argument("--frames", required=True, nargs="+", type=int)
    p.add_argument("--step-detail", action="store_true")
    p.add_argument("--out", required=True, type=Path)
    a = p.parse_args()

    with a.audit.open() as f:
        rows = list(csv.DictReader(f))
    curve = list(load_curve(a.case / "loading.log").values())
    pipe = pipeline_for(a.case / "loading.dump", keep_ptm_structure=True)
    fig = plt.figure(figsize=(max(9, 2.1 * len(a.frames)), 6.9), layout="constrained")
    grid = fig.add_gridspec(2, len(a.frames), height_ratios=[1, 1.4])
    ax = fig.add_subplot(grid[0, :])
    mark = [rows[i] for i in a.frames]
    ex = np.array([100*float(r["v_strain"]) for r in mark])
    ey = [float(r["v_axial"]) for r in mark]
    limits = (ex.min() - .12, ex.max() + .12)
    visible = [r for r in curve if limits[0] <= 100 * r["v_strain"] <= limits[1]]
    ax.plot([100 * r["v_strain"] for r in visible], [r["v_axial"] for r in visible],
            color="#2166ac", lw=1.4)
    ax.scatter(ex, ey, c="#a33222", zorder=3)
    ax.set(xlim=limits, xlabel="Engineering strain (%)",
           ylabel="Nominal axial stress (GPa)", title=f"{a.case.name} · first-event audit")
    ax.grid(alpha=.2)

    for k, i in enumerate(a.frames):
        data = pipe.compute(i)
        panel = fig.add_subplot(grid[1, k])
        xyz = np.asarray(data.particles.positions) / 10
        hcp = np.asarray(data.particles["PTM Structure"]) == 2
        panel.scatter(xyz[::12, 0], xyz[::12, 2], s=2, color="#b1c0c8", alpha=.30, linewidth=0)
        panel.scatter(xyz[hcp, 0], xyz[hcp, 2], s=7, color="#f0a224", alpha=.9, linewidth=0)
        left = data.cell[0, 3] / 10
        lx = data.cell[0, 0] / 10
        film_bounds = Rectangle((left, -0.12), lx, 6.87, transform=panel.transData)
        panel.axvline(left + lx, color="#b1c0c8", linestyle=":", linewidth=0.8)
        for line in data.dislocations.lines:
            points = np.asarray(line.points) / 10
            for shift in (-lx, 0, lx):
                trace, = panel.plot(points[:, 0] + shift, points[:, 2], color="#a33222", lw=1.7)
                trace.set_clip_path(film_bounds)
        x_limits = (-.1, 6.1)
        z_limits = (-.12, 6.75)
        if a.step_detail:
            cut = left + .75 * lx
            panel.axvline(cut, color="#697983", linestyle="--", linewidth=0.9)
            x_limits = (cut - .9, cut + 1.15)
            z_limits = (4.3, 6.55)
        panel.set(xlim=x_limits, ylim=z_limits, aspect="equal",
                  title=f"{ex[k]:.2f}% · frame {i}\nDXA lines: {len(data.dislocations.lines)}",
                  xlabel="[110] (nm)")
        if k == 0:
            panel.set_ylabel("[001] (nm)")
        else:
            panel.tick_params(labelleft=False)
    panel.scatter([], [], s=14, color="#f0a224", label="PTM HCP")
    panel.plot([], [], color="#a33222", label="DXA line")
    if a.step_detail:
        panel.plot([], [], color="#697983", linestyle="--", label="Initial cut (affine)")
    panel.legend(loc="lower left", frameon=False, fontsize=8)
    fig.savefig(a.out, dpi=185)


if __name__ == "__main__":
    main()
