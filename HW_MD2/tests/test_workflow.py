"""Regression checks for the calculation/storage/visualization boundaries."""
from dataclasses import replace
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

from analysis import com_removed_displacement, com_removed_velocity
from md1d import SimulationConfig, simulate
from trajectory_io import load_trajectory, save_trajectory

ROOT = Path(__file__).resolve().parents[1]


class WorkflowTests(unittest.TestCase):
    def test_compute_without_matplotlib_or_files(self):
        code = '''
import builtins
from pathlib import Path
original_import = builtins.__import__
def isolated_import(name, *args, **kwargs):
    if name.split('.')[0] == 'matplotlib':
        raise ImportError('matplotlib intentionally unavailable')
    return original_import(name, *args, **kwargs)
builtins.__import__ = isolated_import
from md1d import simulate, SimulationConfig
simulate('free', SimulationConfig(nsteps=3))
assert not list(Path('.').iterdir()), 'simulate wrote files'
'''
        with tempfile.TemporaryDirectory() as directory:
            env = dict(os.environ, PYTHONPATH=str(ROOT), PYTHONDONTWRITEBYTECODE='1')
            subprocess.run([sys.executable, '-c', code], cwd=directory, env=env, check=True)

    def test_self_contained_roundtrip_and_uneven_sampling(self):
        config = SimulationConfig(n=8, mass=2.0, epsilon=0.8, sigma=1.2,
                                  nsteps=11, output_every=4)
        result = simulate('periodic', config)
        np.testing.assert_array_equal(result.thermal[:, 0], [0, 4, 8, 11])
        with tempfile.TemporaryDirectory() as directory:
            save_trajectory(result, directory)
            Path(directory, 'metadata.json').unlink()  # NPZ must suffice on its own.
            loaded = load_trajectory(directory)
        self.assertEqual(loaded.metadata, result.metadata)
        for name in ('positions', 'positions_unwrapped', 'velocities', 'forces', 'thermal'):
            np.testing.assert_array_equal(getattr(loaded, name), getattr(result, name))
        np.testing.assert_allclose(com_removed_displacement(loaded).mean(axis=1), 0, atol=1e-15)
        np.testing.assert_allclose(com_removed_velocity(loaded).mean(axis=1), 0, atol=1e-15)

    def test_unwrapped_positions_count_multiple_crossings(self):
        # A synthetic one-step case tests bookkeeping, not physical stability.
        for sign in (-1, 1):
            config = SimulationConfig(nsteps=1, output_every=10)
            config = replace(config, kick_velocity=sign * (2 * config.box_length + .05) / config.dt)
            result = simulate('periodic', config)
            expected = np.zeros(config.n)
            expected[config.n // 2] = config.kick_velocity * config.dt
            np.testing.assert_allclose(
                result.positions_unwrapped[1] - result.positions_unwrapped[0], expected, atol=1e-12)

    def test_plotting_does_not_import_solver_or_modify_input(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            archive_dir = directory / 'input'
            save_trajectory(simulate('free', SimulationConfig(nsteps=11)), archive_dir)
            archive = archive_dir / 'trajectory.npz'
            before = archive.read_bytes()
            code = '''
import builtins
import sys
original_import = builtins.__import__
def isolated_import(name, *args, **kwargs):
    if name == 'md1d':
        raise ImportError('solver intentionally unavailable')
    return original_import(name, *args, **kwargs)
builtins.__import__ = isolated_import
from trajectory_io import load_trajectory
from visualize import plot_spacetime
plot_spacetime([load_trajectory(sys.argv[1])], sys.argv[2])
'''
            env = dict(os.environ, PYTHONPATH=str(ROOT), PYTHONDONTWRITEBYTECODE='1',
                       MPLCONFIGDIR=str(directory / 'mpl-cache'))
            subprocess.run([sys.executable, '-c', code, str(archive), str(directory / 'plots')],
                           cwd=directory, env=env, check=True)
            self.assertEqual(archive.read_bytes(), before)
            self.assertEqual(len(list((directory / 'plots').glob('*.png'))), 2)


if __name__ == '__main__':
    unittest.main()
