"""Render the spatial distribution of OVITO-classified planar faults."""
from pathlib import Path
from tempfile import TemporaryDirectory

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from ovito.io import import_file
from ovito.modifiers import (DeleteSelectedModifier, IdentifyFCCPlanarFaultsModifier,
                             PolyhedralTemplateMatchingModifier)

from render_ovito import panel, viewport

COLORS = np.array([(0.82, 0.87, 0.89), (0.60, 0.60, 0.60), (0.89, 0.38, 0.12),
                   (0.10, 0.44, 0.66), (0.53, 0.32, 0.58)])


def color_planar_faults(frame, data):
    types = data.particles["Planar Fault Type"]
    x, y, z = data.particles.positions.T
    backdrop = ((y > 48) & (z < 4)) | ((y > 48) & (z > 57)) | ((x < 3) & (y > 48))
    data.particles_.create_property("Color", data=COLORS[types])
    data.particles_.create_property("Selection", data=~(backdrop | (types >= 2)))


def main():
    fig, axes = plt.subplots(2, 3, figsize=(9.6, 6.5), layout="constrained")
    view = viewport((32, 13, 33), (-.65, 2.0, -.48), 43)
    with TemporaryDirectory() as temp:
        for row, case in zip(axes, ("step-01", "step-02")):
            pipe = import_file(str(Path("pilot/data") / case / "loading.dump"))
            pipe.modifiers.append(PolyhedralTemplateMatchingModifier(
                output_orientation=True, output_interatomic_distance=True))
            pipe.modifiers.append(IdentifyFCCPlanarFaultsModifier())
            pipe.modifiers.append(color_planar_faults)
            pipe.modifiers.append(DeleteSelectedModifier())
            data = pipe.compute()
            data.cell.vis.enabled = False
            data.particles.vis.radius = .90
            for ax, frame in zip(row, (130, 170, 240)):
                panel(pipe, view, frame, ax, temp, size=(740, 660))
                ax.set_title(f"{case}   |   {frame * .05:.2f}%", fontsize=10)
    fig.legend(handles=[
        Line2D([], [], ls="", marker="o", color=COLORS[2], label="Intrinsic stacking fault"),
        Line2D([], [], ls="", marker="o", color=COLORS[3], label="Coherent twin boundary"),
        Line2D([], [], ls="", marker="o", color=COLORS[4], label="Multilayer fault"),
    ], loc="outside lower center", ncol=3, frameon=False, fontsize=9)
    fig.savefig("结课报告/图/08-OVITO层错与孪晶界.png", dpi=260)


if __name__ == "__main__":
    main()
