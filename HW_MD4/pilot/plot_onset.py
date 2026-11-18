"""Publication panels for pre-yield stress and independently sampled event metrics."""
import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from analyze import load_curve

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.labelcolor": "#29333a", "xtick.color": "#29333a",
                     "ytick.color": "#29333a"})
root = Path("pilot")
fig = plt.figure(figsize=(9.5, 6.7), layout="constrained")
gs = fig.add_gridspec(2, 3, width_ratios=(1.15, 1.15, 1))
ax = fig.add_subplot(gs[:, :2])
strain_ax = fig.add_subplot(gs[0, 2])
stress_ax = fig.add_subplot(gs[1, 2])
colors = {"flat": "#245474", "step": "#b05337"}
for idx, case in enumerate(("flat-01", "flat-02", "step-01", "step-02")):
    group, number = case.split("-")
    color = colors[group]
    curve = list(load_curve(root / "data" / case / "loading.log").values())
    with (root / "audit" / case / "loading_audit.csv").open() as stream:
        frames = list(csv.DictReader(stream))
    first = next(row for row in frames if int(row["line_count"]))
    peak = max(curve, key=lambda row: row["v_axial"])
    ax.plot([row["v_strain"]*100 for row in curve], [row["v_axial"] for row in curve],
            color=color, ls="-" if number == "01" else "--", lw=1.6)
    px, py = 100 * peak["v_strain"], peak["v_axial"]
    dx, dy = 100 * float(first["v_strain"]), float(first["v_axial"])
    ax.scatter(px, py, color=color, marker="D", s=29, zorder=4)
    ax.scatter(dx, dy, color=color, marker="o", s=51, facecolors="white", zorder=4)
    xpos = (0 if group == "flat" else 1) + (-.10 if number == "01" else .10)
    strain_ax.plot([xpos, xpos], [px, dx], color=color, alpha=.48, lw=1.3)
    strain_ax.scatter(xpos, px, color=color, marker="D", s=29, zorder=3)
    strain_ax.scatter(xpos, dx, color=color, marker="o", facecolors="white", s=51, zorder=3)
    stress_ax.scatter(xpos, py, marker="D", color=color, s=32, zorder=3)
    if group == "step":
        corrected = py * 18/17
        stress_ax.plot([xpos, xpos], [py, corrected], color=color, alpha=.48, lw=1.3)
        stress_ax.scatter(xpos, corrected, marker="s", color=color,
                          facecolors="white", s=38, zorder=3)

ax.set(xlim=(0, 9), ylim=(-1.1, 6.5), xlabel="Engineering strain (%)",
       ylabel="Nominal axial stress (GPa)", title="(a) Stress trajectory and first detection")
ax.legend(handles=[
    Line2D([], [], color=colors["flat"], label="Flat"),
    Line2D([], [], color=colors["step"], label="Step"),
    Line2D([], [], color="#43525c", marker="D", ls="", label="Peak stress"),
    Line2D([], [], color="#43525c", marker="o", markerfacecolor="white", ls="", label="First DXA"),
], loc="upper left", frameon=False, ncol=2)
strain_ax.set(ylabel="Engineering strain (%)", ylim=(4.55, 8.05),
              title="(b) Onset strain")
stress_ax.set(ylabel="Peak stress (GPa)", ylim=(5.12, 6.11),
              title="(c) Section sensitivity")
stress_ax.legend(handles=[Line2D([], [], ls="", marker="s", markerfacecolor="white",
                                color=colors["step"], label="Minimum section")],
                 loc="upper right", fontsize=8, frameon=False)
for panel in (strain_ax, stress_ax):
    panel.set_xticks([0, 1], ["Flat", "Step"])
    panel.set_xlim(-.48, 1.48)
for panel in (ax, strain_ax, stress_ax):
    panel.yaxis.grid(True, lw=.5, color="#dce4e9")
    panel.set_axisbelow(True)

fig.savefig("结课报告/图/03-应力与首事件.png", dpi=300)
fig.savefig("结课报告/图/03-应力与首事件.svg")
