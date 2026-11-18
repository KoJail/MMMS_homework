"""Diagnose the coupling of film thickness and axial stress before DXA detections."""
import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from ovito.io import import_file

from analyze import load_curve


def detrend(times, values):
    return values - np.polyval(np.polyfit(times, values, 1), times)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("cases", nargs=4, type=Path)
    args = parser.parse_args()
    fig, axes = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(9, 5.7),
                             layout="constrained")
    with args.csv.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("case", "time_ps", "p99_minus_p01_z_angstrom", "axial_gpa",
                         "detrended_thickness_angstrom", "detrended_axial_gpa"))
        for case, ax in zip(args.cases, axes.flat):
            pipe = import_file(str(case / "loading.dump"))
            curve = load_curve(case / "loading.log")
            times, thickness, axial = [], [], []
            for i in range(99):
                frame = pipe.compute(i)
                row = curve[frame.attributes["Timestep"]]
                low, high = np.quantile(frame.particles.positions[:, 2], (.01, .99))
                times.append(row["Time"])
                thickness.append(high - low)
                axial.append(row["v_axial"])
            z = detrend(times, thickness)
            s = detrend(times, axial)
            for row in zip(times, thickness, axial, z, s):
                writer.writerow((case.name, *row))
            corr = np.corrcoef(z, s)[0, 1]
            ax.plot(times, s / s.std(), color="#22639a", label="Axial stress")
            ax.plot(times, z / z.std(), color="#b24927", ls="--", label="Film thickness")
            ax.set_title(f"{case.name} · r = {corr:.3f}")
            ax.grid(alpha=.2)
    axes[0, 0].legend(frameon=True, loc="lower center", ncol=2)
    fig.supxlabel("Time after loading starts (ps)")
    fig.supylabel("Detrended signal (standard deviation units)")
    fig.savefig(args.out, dpi=185)


if __name__ == "__main__":
    main()
