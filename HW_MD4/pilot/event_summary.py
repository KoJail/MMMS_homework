"""Audit first detections and consecutive spatial matches in the eight trajectories."""
import csv
from itertools import permutations, product
from pathlib import Path

import numpy as np

from analyze import load_curve, pipeline_for, write_rows

ROOT = Path("pilot")
REPORT = Path("结课报告/数据")
A = 3.615
CRYSTAL_TO_BOX = np.array([[1, 1, 0], [-1, 1, 0], [0, 0, 1]], dtype=float)
CRYSTAL_TO_BOX /= np.linalg.norm(CRYSTAL_TO_BOX, axis=1)[:, None]


def slip_systems():
    directions = sorted(set(v for m in permutations((1, 1, 2))
                            for signs in product((-1, 1), repeat=3)
                            if (v := tuple(a*b for a, b in zip(m, signs)))[0] > 0))
    planes = [(1, y, z) for y, z in product((-1, 1), repeat=2)]
    return [(n, s) for n in planes for s in directions if np.dot(n, s) == 0]


def nominal_diagonal(row):
    # Use the same material thickness to normalize all three box-pressure components.
    return -row["Lz"] / (18 * A) * 1e-4 * np.array([row[key] for key in ("Pxx", "Pyy", "Pzz")])


def classify(points, cell, low, high, stepped):
    x_fraction = (points[:, 0] - cell[0, 3]) / cell[0, 0]
    upper = points[:, 2] >= high - 6
    if stepped:
        for edge in (.25, .75):
            distance = abs((x_fraction - edge + .5) % 1 - .5) * cell[0, 0]
            if np.any(upper & (distance <= 5)):
                return f"upper step at x/Lx={edge:.2f}"
    if upper.any():
        return "upper surface"
    if points[:, 2].min() <= low + 6:
        return "lower surface"
    return "film interior"


def candidates(data, stepped):
    low, high = np.quantile(data.particles.positions[:, 2], (.005, .995))
    output = []
    for line in data.dislocations.lines:
        b = np.abs(np.array(line.true_burgers_vector))
        if not np.allclose(np.sort(b), [1/6, 1/6, 1/3], atol=.002):
            continue
        points = np.asarray(line.points)
        output.append((points.mean(axis=0), classify(points, data.cell, low, high, stepped), line.length))
    return output


def adjacent(left, right, lengths):
    delta = left - right
    delta[:2] -= lengths * np.round(delta[:2] / lengths)
    return np.linalg.norm(delta) <= 8


def track(pipe, first, stepped):
    frames = []
    for i in range(first, min(first + 30, pipe.source.num_frames)):
        frame = pipe.compute(i)
        lines = candidates(frame, stepped)
        frames.append(lines)
        if len(frames) < 3:
            continue
        lengths = np.array([frame.cell[0, 0], frame.cell[1, 1]])
        for a, b, c in product(*frames[-3:]):
            if adjacent(a[0], b[0], lengths) and adjacent(b[0], c[0], lengths):
                return frames[0], i - 2, a[1]
    raise ValueError("No three-frame spatial match in candidate window")


def main():
    summaries, shear = [], []
    systems = slip_systems()
    assert len(systems) == 12
    with (REPORT / "加载起点对照.csv").open() as f:
        for row in csv.DictReader(f):
            delay = float(row["delay_ps"])
            name = row["case"]
            case = name if delay == 0 else f"phase-{name}-{round(delay*1000)}"
            logs = load_curve(ROOT / "data" / case / "loading.log")
            pipe = pipeline_for(ROOT / "data" / case / "loading.dump")
            first = int(row["first_dxa_frame"])
            lines, matched, location = track(pipe, first, name.startswith("step"))
            stride = 20 if delay else 100
            diagonal = nominal_diagonal(logs[first * stride])
            max_shear = 0.
            for n, s in systems:
                nv = CRYSTAL_TO_BOX @ (np.array(n) / np.sqrt(3))
                sv = CRYSTAL_TO_BOX @ (np.array(s) / np.sqrt(6))
                assert abs(np.dot(nv, sv)) < 1e-12
                tau = float(np.dot(sv * nv, diagonal))
                shear.append(dict(case=case, plane=str(n), direction=str(s), tau_diagonal_gpa=tau))
                max_shear = max(max_shear, abs(tau))
            summaries.append(dict(case=case, delay_ps=delay,
                last_clear_strain_pct=row["prior_frame_strain_pct"], first_strain_pct=row["first_dxa_strain_pct"],
                first_locations="; ".join(dict.fromkeys(x[1] for x in lines)), first_axial_gpa=diagonal[0],
                first_transverse_gpa=diagonal[1], first_normal_gpa=diagonal[2],
                max_abs_diagonal_rss_gpa=max_shear, peak_gpa=row["peak_gpa"],
                matched_start_strain_pct=logs[matched * stride]["v_strain"] * 100, matched_location=location))
    write_rows(REPORT / "首事件与分解剪应力.csv", summaries)
    write_rows(REPORT / "候选滑移系投影.csv", shear)
    for row in summaries:
        print(row)


if __name__ == "__main__":
    main()
