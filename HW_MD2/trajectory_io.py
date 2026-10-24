"""Self-contained trajectory archives and human-readable exports; no plotting."""
import csv
from dataclasses import dataclass
import json
from pathlib import Path

import numpy as np

THERMAL_COLUMNS = ('step', 'time', 'potential_energy', 'kinetic_energy',
                   'total_energy', 'temperature_reduced')


@dataclass
class Trajectory:
    positions: np.ndarray
    velocities: np.ndarray
    forces: np.ndarray
    thermal: np.ndarray
    positions_unwrapped: np.ndarray
    metadata: dict

    @property
    def boundary(self):
        return self.metadata['boundary']

    @property
    def time(self):
        return self.thermal[:, 1]

    def validate(self):
        if self.metadata.get('schema_version') != 1:
            raise ValueError('Unsupported trajectory schema version')
        if self.boundary not in ('free', 'periodic'):
            raise ValueError('Invalid boundary in trajectory metadata')
        if self.metadata.get('thermal_columns') != list(THERMAL_COLUMNS):
            raise ValueError('Unexpected thermal column definitions')
        if self.positions.ndim != 2 or self.positions.shape[0] < 2:
            raise ValueError('Positions must be a 2D history including initial and final states')
        shape = self.positions.shape
        if shape[1] != self.metadata['config']['n']:
            raise ValueError('Particle count does not match metadata')
        for array in (self.velocities, self.forces, self.positions_unwrapped):
            if array.shape != shape:
                raise ValueError('All particle histories must have the same 2D shape')
        if self.thermal.shape != (shape[0], len(THERMAL_COLUMNS)):
            raise ValueError('Invalid thermal history shape')
        for array in (self.positions, self.velocities, self.forces,
                      self.positions_unwrapped, self.thermal):
            if not np.isfinite(array).all():
                raise ValueError('Trajectory contains non-finite values')
        if np.any(np.diff(self.time) <= 0):
            raise ValueError('Saved times must be strictly increasing')


def save_trajectory(trajectory, directory):
    """NPZ is authoritative; JSON, CSV and XYZ are convenient exports."""
    trajectory.validate()
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    metadata_json = json.dumps(trajectory.metadata, indent=2, allow_nan=False)
    np.savez_compressed(
        directory / 'trajectory.npz', positions=trajectory.positions,
        velocities=trajectory.velocities, forces=trajectory.forces,
        thermal=trajectory.thermal, positions_unwrapped=trajectory.positions_unwrapped,
        box_length=trajectory.metadata['box_length'], metadata_json=metadata_json,
    )
    (directory / 'metadata.json').write_text(metadata_json + '\n', encoding='utf-8')
    with (directory / 'thermal.csv').open('w', newline='', encoding='utf-8') as fp:
        writer = csv.writer(fp)
        writer.writerow(THERMAL_COLUMNS)
        writer.writerows(trajectory.thermal)
    # Retain the original plain XYZ export. Periodic-cell metadata is in the NPZ.
    with (directory / 'trajectory.xyz').open('w', encoding='utf-8') as fp:
        for positions, row in zip(trajectory.positions, trajectory.thermal):
            fp.write(f'{len(positions)}\nstep={int(row[0])} time={row[1]:.6f} '
                     f'boundary={trajectory.boundary}\n')
            for x in positions:
                fp.write(f'Ar {x:.9f} 0.0 0.0\n')


def load_trajectory(path):
    """Load without pickle or access to the simulation module or its defaults."""
    path = Path(path)
    if path.is_dir():
        path = path / 'trajectory.npz'
    with np.load(path, allow_pickle=False) as data:
        if 'metadata_json' not in data or 'positions_unwrapped' not in data:
            raise ValueError(f'{path}: legacy archive lacks metadata; rerun md1d.py')
        trajectory = Trajectory(
            positions=data['positions'], velocities=data['velocities'],
            forces=data['forces'], thermal=data['thermal'],
            positions_unwrapped=data['positions_unwrapped'],
            metadata=json.loads(str(data['metadata_json'].item())),
        )
    trajectory.validate()
    return trajectory
