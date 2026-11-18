"""Distinguish intrinsic stacking faults and twin boundaries with OVITO."""
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from ovito.io import import_file
from ovito.modifiers import IdentifyFCCPlanarFaultsModifier, PolyhedralTemplateMatchingModifier

cases = ("flat-01", "flat-02", "step-01", "step-02")
labels = ("Isolated / indeterminate HCP", "Intrinsic stacking fault",
          "Coherent twin boundary", "Multilayer faults")
colors = ("#a3a7a8", "#a44b35", "#255475", "#8664a2")
fig, axes = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(9.4, 6.0),
                         layout="constrained")
with Path("结课报告/数据/平面缺陷分类.csv").open("w", newline="") as stream:
    writer = csv.writer(stream)
    writer.writerow(("case", "strain_percent", "atom_count", "indeterminate_hcp_atoms",
                     "intrinsic_fault_atoms", "coherent_twin_boundary_atoms", "multilayer_fault_atoms"))
    for case, ax in zip(cases, axes.flat):
        path = Path("pilot/data") / case / "loading.dump"
        pipe = import_file(str(path))
        pipe.modifiers.append(PolyhedralTemplateMatchingModifier(output_interatomic_distance=True,
                                                                  output_orientation=True))
        pipe.modifiers.append(IdentifyFCCPlanarFaultsModifier())
        values = []
        strain = []
        for i in range(90, pipe.source.num_frames, 2):
            data = pipe.compute(i)
            count = data.particles.count
            defects = np.bincount(data.particles["Planar Fault Type"], minlength=5)[1:]
            assert len(defects) == 4
            x = data.attributes["Timestep"] * .0005
            writer.writerow((case, x, count, *defects))
            strain.append(x)
            values.append(defects * 100 / count)
        for y, name, color in zip(np.array(values).T, labels, colors):
            ax.plot(strain, y, label=name, color=color, lw=1.6)
        ax.set_title(case, loc="left", fontweight="bold", fontsize=11)
        ax.set_xlim(4.5, 12)
        ax.set_axisbelow(True)
        ax.grid(axis="y", color="#dce4e9", lw=.6)
        ax.spines[["top", "right"]].set_visible(False)
fig.supxlabel("Engineering strain (%)")
fig.supylabel("Atoms assigned to planar-defect classes (%)")
fig.legend(*axes[0, 0].get_legend_handles_labels(), loc="outside upper center",
           frameon=False, ncol=2, fontsize=9)
fig.savefig("结课报告/图/07-层错与孪晶界.png", dpi=270)
