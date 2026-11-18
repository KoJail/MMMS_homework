"""Trace the positions and Burgers vectors of individual lines around a candidate."""
import argparse
from pathlib import Path

import numpy as np

from analyze import load_curve, pipeline_for


def inspect(dump, start, stop):
    pipe = pipeline_for(dump)
    log = load_curve(dump.with_name("loading.log"))
    for index in range(start, min(stop + 1, pipe.source.num_frames)):
        frame = pipe.compute(index)
        step = frame.attributes["Timestep"]
        xyz = np.asarray(frame.particles.positions)
        low, high = np.quantile(xyz[:, 2], [.005, .995])
        thermo = log[step]
        print(f"frame {index}   strain {100 * thermo['v_strain']:.2f}%"
              f"   stress {thermo['v_axial']:.3f} GPa"
              f"   HCP {frame.attributes['PolyhedralTemplateMatching.counts.HCP']}")
        for line in frame.dislocations.lines:
            points = np.asarray(line.points)
            center = points.mean(axis=0)
            for axis in (0, 1):
                origin = frame.cell[axis, 3]
                center[axis] = origin + (center[axis] - origin) % frame.cell[axis, axis]
            print(" ", "b =", np.round(line.true_burgers_vector, 4),
                  "length (A) =", round(line.length, 2),
                  "center (A) =", np.round(center, 2),
                  "z-range (A) =", np.round([points[:, 2].min(), points[:, 2].max()], 2),
                  "surface z (A) =", np.round([low, high], 2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dump", type=Path)
    parser.add_argument("start", type=int)
    parser.add_argument("stop", type=int)
    args = parser.parse_args()
    inspect(args.dump, args.start, args.stop)
