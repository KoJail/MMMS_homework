"""Compare loading curves and defect measurements from completed pilot cases."""
import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt

from analyze import load_curve


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("cases", type=Path, nargs="+")
    args = parser.parse_args()
    pilot = Path(__file__).resolve().parent

    fig, axes = plt.subplots(3, 1, figsize=(8, 8), sharex=True, layout="constrained")
    stress, lines, hcp = axes
    for case in args.cases:
        log = load_curve(case / "loading.log")
        with (pilot / "audit" / case.name / "loading_audit.csv").open() as stream:
            frames = list(csv.DictReader(stream))
        expected_steps = list(range(0, 24001, 100))
        if [int(row["step"]) for row in frames] != expected_steps or max(log) != 24000:
            raise ValueError(f"{case.name}: incomplete R2 trajectory")

        morphology, replicate = case.name.split("-")
        color = {"flat": "#2166ac", "step": "#b24927"}[morphology]
        style = "-" if replicate == "01" else "--"
        strain = [100 * row["v_strain"] for row in log.values()]
        stress.plot(strain, [row["v_axial"] for row in log.values()],
                    label=case.name, color=color, linestyle=style, linewidth=1.4)
        sampled_strains = [100 * float(row["v_strain"]) for row in frames]
        lines.plot(sampled_strains, [float(row["total_line_angstrom"]) / 10 for row in frames],
                   color=color, linestyle=style, linewidth=1.2)
        hcp.plot(sampled_strains, [int(row["hcp"]) for row in frames],
                 color=color, linestyle=style, linewidth=1.2)
        candidate = next((row for row in frames if int(row["line_count"])), None)
        if candidate:
            stress.scatter(100 * float(candidate["v_strain"]), float(candidate["v_axial"]),
                           s=38, facecolor="white", edgecolor=color, zorder=3)
    stress.set(ylabel="Nominal axial stress (GPa)",
               title="Pilot loading curves · open circles mark first DXA detections")
    stress.legend(frameon=False, ncol=2)
    lines.set_ylabel("Total DXA line length (nm)")
    hcp.set(xlabel="Engineering strain (%)", ylabel="PTM HCP atoms", xlim=(0, 12))
    for ax in axes:
        ax.grid(alpha=.2)
    fig.savefig(args.out, dpi=185)


if __name__ == "__main__":
    main()
