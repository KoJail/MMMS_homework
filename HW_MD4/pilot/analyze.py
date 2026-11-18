"""Read LAMMPS pilot cases and write per-frame audit tables."""
import argparse
import csv
from pathlib import Path

import numpy as np
from ovito.io import import_file
from ovito.modifiers import DislocationAnalysisModifier, PolyhedralTemplateMatchingModifier


def load_curve(log_path):
    rows = {}
    columns = ()
    for line in log_path.read_text().splitlines():
        fields = line.split()
        if fields[:4] == ["Step", "Time", "v_strain", "v_axial"]:
            columns = fields
        elif columns:
            try:
                vals = [float(item) for item in fields]
            except ValueError:
                columns = ()
                continue
            if len(vals) == len(columns):
                rows[int(vals[0])] = dict(zip(columns, vals))
    return rows


def preserve_ptm_structure(frame, data):
    # DXA replaces the standard Structure Type particle property.
    data.particles_.create_property("PTM Structure", data=np.asarray(data.particles["Structure Type"]).copy())


def pipeline_for(path, *, keep_ptm_structure=False):
    pipeline = import_file(str(path))
    pipeline.modifiers.append(PolyhedralTemplateMatchingModifier())
    if keep_ptm_structure:
        pipeline.modifiers.append(preserve_ptm_structure)
    pipeline.modifiers.append(DislocationAnalysisModifier(
        input_crystal_structure=DislocationAnalysisModifier.Lattice.FCC
    ))
    return pipeline


def flatness(data, bins=20):
    x, _, z = np.asarray(data.particles.positions).T
    length = data.cell[0, 0]
    indices = np.floor(np.mod(x - data.cell[0, 3], length) / length * bins).astype(int)
    heights = np.array([np.quantile(z[indices == idx], .975) for idx in range(bins)])
    return float(np.median(heights[:5]) - np.median(heights[5:15])), heights


def read_frame(data, frame, curve_row=None):
    lines = data.dislocations.lines
    largest = max(lines, key=lambda line: line.length, default=None)
    if largest is not None:
        points = np.asarray(largest.points)
        frac_x = np.mod((points[:, 0] - data.cell[0, 3]) / data.cell[0, 0], 1)
        # Distance to the geometrical cut planes. Verify actual step positions visually.
        cut_distance = min(
            np.min(np.abs(np.mod(frac_x - face + .5, 1) - .5))
            for face in (.25, .75)
        ) * data.cell[0, 0]
    else:
        points = None
        cut_distance = None
    row = dict(
        frame=frame, step=data.attributes["Timestep"],
        atom_count=data.particles.count,
        total_line_angstrom=sum(line.length for line in lines),
        line_count=len(lines),
        longest_line_angstrom=largest.length if largest is not None else None,
        longest_line_z_min_angstrom=float(points[:, 2].min()) if points is not None else None,
        longest_line_z_max_angstrom=float(points[:, 2].max()) if points is not None else None,
        distance_to_initial_cut_plane_angstrom=cut_distance,
        fcc=data.attributes.get("PolyhedralTemplateMatching.counts.FCC"),
        hcp=data.attributes.get("PolyhedralTemplateMatching.counts.HCP"),
        burgers="; ".join(str(np.round(line.true_burgers_vector, 4)) for line in lines),
    )
    if curve_row is not None:
        row.update((key, curve_row[key]) for key in ("Time", "v_strain", "v_axial", "v_lateral", "Temp"))
    return row


def write_rows(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def analyze_case(case, out):
    case_out = out / case.name
    case_out.mkdir(parents=True, exist_ok=True)
    rows = []
    for stage in ("initial", "minimized", "equilibrated"):
        dump = case / f"{stage}.dump"
        if not dump.exists():
            continue
        data = pipeline_for(dump).compute()
        row = read_frame(data, 0)
        contrast, heights = flatness(data)
        row.update(stage=stage, terrace_contrast_angstrom=contrast)
        rows.append(row)
        if stage == "equilibrated":
            np.savetxt(case_out / "terrace_profile.csv", heights, delimiter=",", header="height_angstrom")
            print(case.name, f"atoms={data.particles.count}", f"line_count={row['line_count']}", f"contrast={contrast:.3f} A", flush=True)
    if rows:
        write_rows(case_out / "initial_audit.csv", rows)

    dump = case / "loading.dump"
    if dump.exists():
        pipeline = pipeline_for(dump)
        curve = load_curve(case / "loading.log")
        frames = []
        for index in range(pipeline.source.num_frames):
            data = pipeline.compute(index)
            step = data.attributes["Timestep"]
            if step not in curve:
                raise ValueError(f"{case.name}: no thermo record for dump step {step}")
            frames.append(read_frame(data, index, curve[step]))
        write_rows(case_out / "loading_audit.csv", frames)
        first = next((row for row in frames if row["line_count"]), None)
        print(case.name, f"aligned {len(frames)} frames", f"first DXA candidate: {first}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cases", nargs="+", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    for case in args.cases:
        analyze_case(case, args.out)


if __name__ == "__main__":
    main()
