"""Summarize completed high-frequency repeats alongside their pilot references."""
import argparse
import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt

pilot = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(pilot))
from analyze import load_curve


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--csv", type=Path, required=True)
    args = parser.parse_args()
    pairs = [line.split()[:2] for line in (pilot / "phase/active-jobs.txt").read_text().splitlines()]
    fig, axes = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(9, 6),
                             layout="constrained")
    with args.csv.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("case", "delay_ps", "init_axial_gpa", "peak_strain_pct", "peak_gpa",
                         "first_dxa_frame", "first_dxa_strain_pct", "first_dxa_gpa",
                         "prior_frame_strain_pct", "cutoff_pct"))
        for ax, (name, delay_steps) in zip(axes.flat, pairs):
            repeat = f"phase-{name}-{delay_steps}"
            for case, delay in ((name, 0.0), (repeat, int(delay_steps) * .001)):
                audit = pilot / "audit" / case
                if not (audit / "loading_audit.csv").is_file():
                    raise FileNotFoundError(f"{case}: missing audited trajectory")
                curve = list(load_curve(pilot / "data" / case / "loading.log").values())
                with (audit / "loading_audit.csv").open(newline="") as data:
                    frames = list(csv.DictReader(data))
                expected_steps = list(range(0, 18001, 20) if delay else range(0, 24001, 100))
                if [int(row["step"]) for row in frames] != expected_steps:
                    raise ValueError(f"{case}: incomplete trajectory")
                with (audit / "initial_audit.csv").open(newline="") as data:
                    initial = list(csv.DictReader(data))[-1]
                if int(initial["line_count"]) or int(initial["hcp"]):
                    raise ValueError(f"{case}: initial defects require manual review")
                peak = max(curve, key=lambda row: row["v_axial"])
                first = next((i for i, row in enumerate(frames) if int(row["line_count"])), None)
                x = y = before = None
                color = "#b24927" if delay else "#22639a"
                if first is not None:
                    x = 100 * float(frames[first]["v_strain"])
                    y = float(frames[first]["v_axial"])
                    if first > 0:
                        before = 100 * float(frames[first - 1]["v_strain"])
                    ax.scatter(x, y, facecolor="white", edgecolor=color, zorder=3)
                ax.scatter(100 * peak["v_strain"], peak["v_axial"], color=color, marker="D", s=25, zorder=3)
                ax.plot([row["v_strain"] * 100 for row in curve], [row["v_axial"] for row in curve],
                        label=f"Delay {delay:g} ps", color=color, ls="--" if delay else "-")
                writer.writerow((name, delay, curve[0]["v_axial"], 100 * peak["v_strain"], peak["v_axial"],
                                 first, x, y, before, 100 * curve[-1]["v_strain"]))
            ax.set(title=name, xlim=(0, 9))
            ax.spines[["top", "right"]].set_visible(False)
            ax.grid(axis="y", color="#dce4e9", linewidth=.6)
            ax.legend(frameon=False)
    fig.supxlabel("Engineering strain (%)")
    fig.supylabel("Nominal axial stress (GPa)")
    fig.savefig(args.out, dpi=185)


if __name__ == "__main__":
    main()
