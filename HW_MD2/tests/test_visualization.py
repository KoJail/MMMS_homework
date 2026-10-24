"""Physical meaning and offline-export checks for linked visualization."""
from dataclasses import replace
from html.parser import HTMLParser
from pathlib import Path
import tempfile
import unittest

import numpy as np

from analysis import particle_kinetic_energy, time_edges, view_fields
from interactive import dashboard_data, write_dashboard
from md1d import SimulationConfig, simulate


class ScriptInspector(HTMLParser):
    def __init__(self):
        super().__init__()
        self.sources = []

    def handle_starttag(self, tag, attrs):
        if tag == 'script':
            self.sources.append(dict(attrs).get('src'))


class VisualizationTests(unittest.TestCase):
    def test_internal_kinetic_energy_matches_thermal_definition(self):
        result = simulate('free', SimulationConfig(n=8, mass=2.0, nsteps=11, output_every=4))
        internal = particle_kinetic_energy(result).sum(axis=1)
        com_energy = 0.5 * 2.0 * 8 * result.velocities.mean(axis=1) ** 2
        np.testing.assert_allclose(internal, result.thermal[:, 3] - com_energy, atol=1e-15)
        np.testing.assert_allclose(internal, 0.5 * 7 * result.thermal[:, 5], atol=1e-15)

    def test_fields_are_invariant_under_uniform_translation_and_boost(self):
        result = simulate('periodic', SimulationConfig(nsteps=20, output_every=4))
        boosted = replace(result, velocities=result.velocities + 2.3,
                          positions_unwrapped=result.positions_unwrapped + 7 + 2.3 * result.time[:, None])
        for name, field in view_fields(result).items():
            np.testing.assert_allclose(field, view_fields(boosted)[name], atol=1e-14)

    def test_actual_time_edges_and_different_sampling_grids(self):
        np.testing.assert_allclose(time_edges(np.array([0, .008, .016, .022])),
                                   [0, .004, .012, .019, .022])
        a = simulate('free', SimulationConfig(nsteps=11, output_every=4))
        b = simulate('periodic', SimulationConfig(nsteps=11, output_every=3))
        data = dashboard_data([a, b])
        np.testing.assert_array_equal(data['time'], np.unique(np.r_[a.time, b.time]))
        self.assertEqual(data['cases'][1]['time_edges'][-1], b.time[-1])
        for name in ('displacement', 'velocity', 'kinetic'):
            self.assertGreaterEqual(data['limits'][name],
                                    max(np.abs(view_fields(item)[name]).max() for item in [a, b]))

    def test_standalone_html_contains_library_and_data_without_external_scripts(self):
        result = simulate('free', SimulationConfig(nsteps=11, output_every=4))
        with tempfile.TemporaryDirectory() as directory:
            path = write_dashboard([result], Path(directory) / 'dashboard.html')
            html = path.read_text(encoding='utf-8')
        parser = ScriptInspector()
        parser.feed(html)
        self.assertEqual(parser.sources, [None, None, None])
        self.assertIn('plotly.js', html)
        self.assertIn('"positions_unwrapped":', html)
        self.assertNotIn('__TRAJECTORY_DATA__', html)


if __name__ == '__main__':
    unittest.main()
