"""1D Lennard-Jones dynamics; computation and command-line simulation only.

Velocity Verlet follows example/MD1HO.m; the pair loop follows ForceCal.m.
All quantities use Lennard-Jones reduced units.
"""
import argparse
from dataclasses import asdict, dataclass
from pathlib import Path
import platform

import numpy as np

from analysis import energy_drift
from trajectory_io import THERMAL_COLUMNS, Trajectory, save_trajectory


@dataclass(frozen=True)
class SimulationConfig:
    n: int = 10
    epsilon: float = 1.0
    sigma: float = 1.0
    mass: float = 1.0
    dt: float = 0.002
    nsteps: int = 10000
    output_every: int = 10
    rcut_sigma: float = 2.5
    spacing_sigma: float = 2.0 ** (1.0 / 6.0)
    kick_velocity: float = 0.2

    @property
    def rcut(self):
        return self.rcut_sigma * self.sigma

    @property
    def spacing(self):
        return self.spacing_sigma * self.sigma

    @property
    def box_length(self):
        return self.n * self.spacing

    def validate(self, boundary):
        if boundary not in ('free', 'periodic'):
            raise ValueError('boundary must be free or periodic')
        for name in ('n', 'nsteps', 'output_every'):
            value = getattr(self, name)
            if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
                raise ValueError(f'{name} must be a positive integer')
        if self.n <= 5:
            raise ValueError('The assignment requires more than five particles')
        for name in ('epsilon', 'sigma', 'mass', 'dt', 'rcut_sigma', 'spacing_sigma'):
            value = getattr(self, name)
            if not np.isfinite(value) or value <= 0:
                raise ValueError(f'{name} must be finite and positive')
        if not np.isfinite(self.kick_velocity):
            raise ValueError('kick_velocity must be finite')
        if boundary == 'periodic' and self.box_length <= 2 * self.rcut:
            raise ValueError('Minimum image requires box_length > 2*rcut')


def raw_lj(r, config):
    """Return potential and outward radial force, -dU/dr."""
    q6 = (config.sigma / r) ** 6
    q12 = q6 * q6
    return (4.0 * config.epsilon * (q12 - q6),
            24.0 * config.epsilon / r * (2.0 * q12 - q6))


def shifted_force_lj(r, config):
    """Return the shifted-force pair potential and radial force (zero at rc)."""
    if r >= config.rcut:
        return 0.0, 0.0
    u, fr = raw_lj(r, config)
    uc, fc = raw_lj(config.rcut, config)
    return u - uc + (r - config.rcut) * fc, fr - fc


def evaluate_forces(x, boundary, config):
    """Compute pair forces and count each pair's potential once."""
    f = np.zeros(config.n)
    ep = 0.0
    for i in range(config.n - 1):
        for j in range(i + 1, config.n):
            dx = x[j] - x[i]
            if boundary == 'periodic':
                dx -= config.box_length * np.floor(dx / config.box_length + 0.5)
            r = abs(dx)
            if r < 1.0e-10:
                raise ValueError('Overlapping particles; cannot evaluate LJ force')
            if r >= config.rcut:
                continue
            u, fr = shifted_force_lj(r, config)
            fij = -fr * np.sign(dx)
            f[i] += fij
            f[j] -= fij
            ep += u
    return f, ep


def energies(v, ep, config):
    ek = 0.5 * config.mass * np.dot(v, v)
    v_thermal = v - np.mean(v)
    ek_thermal = 0.5 * config.mass * np.dot(v_thermal, v_thermal)
    return ek, ep + ek, 2.0 * ek_thermal / (config.n - 1)


def simulate(boundary, config=None):
    """Return sampled data and metadata without writing files or drawing plots."""
    if config is None:
        config = SimulationConfig()
    config.validate(boundary)
    x = np.arange(config.n, dtype=float) * config.spacing
    v = np.zeros(config.n)
    v[config.n // 2] = config.kick_velocity
    f, ep = evaluate_forces(x, boundary, config)
    a = f / config.mass
    # Count every boundary crossing during integration, not between saved frames.
    images = np.zeros(config.n, dtype=np.int64)

    nsamples = config.nsteps // config.output_every + 1
    if config.nsteps % config.output_every:
        nsamples += 1  # Retain the final state, even with uneven sampling.
    X = np.zeros((nsamples, config.n))
    XU = np.zeros_like(X)
    V = np.zeros_like(X)
    F = np.zeros_like(X)
    thermal = np.zeros((nsamples, len(THERMAL_COLUMNS)))

    def record(k, step, ep_current):
        ek, et, temp = energies(v, ep_current, config)
        X[k] = x
        XU[k] = x + images * config.box_length if boundary == 'periodic' else x
        V[k] = v
        F[k] = f
        thermal[k] = (step, step * config.dt, ep_current, ek, et, temp)

    record(0, 0, ep)
    idx = 1
    for step in range(1, config.nsteps + 1):
        x = x + v * config.dt + 0.5 * a * config.dt * config.dt
        if boundary == 'periodic':
            images += np.floor(x / config.box_length).astype(np.int64)
            x = np.mod(x, config.box_length)
        fn, ep = evaluate_forces(x, boundary, config)
        an = fn / config.mass
        v = v + 0.5 * config.dt * (a + an)
        f, a = fn, an
        if step % config.output_every == 0 or step == config.nsteps:
            record(idx, step, ep)
            idx += 1

    metadata = {
        'schema_version': 1,
        'boundary': boundary,
        'config': asdict(config),
        'rcut': config.rcut,
        'spacing': config.spacing,
        'box_length': config.box_length,
        'dimensions': 1,
        'integrator': 'velocity_verlet',
        'potential': 'lennard_jones_shifted_force',
        'ensemble': 'NVE',
        'initialization': {
            'positions': 'x[i] = i * spacing',
            'velocities': 'single_particle_kick; all other velocities zero',
            'kick_particle_index_zero_based': config.n // 2,
            'random_seed': None,
            'center_of_mass_velocity_removed': False,
        },
        'units': {
            'system': 'Lennard-Jones reduced units',
            'length': 'sigma', 'energy': 'epsilon', 'mass': 'm',
            'time': 'sigma * sqrt(m / epsilon)',
            'velocity': 'sqrt(epsilon / m)', 'force': 'epsilon / sigma',
            'temperature': 'epsilon / k_B', 'k_B': 1.0,
        },
        'array_axes': ['sample', 'particle'],
        'particle_order': 'fixed; column i is particle i (zero based)',
        'positions_convention': 'wrapped_[0,L)' if boundary == 'periodic' else 'unbounded',
        'positions_unwrapped_convention': 'continuous laboratory coordinates',
        'thermal_columns': list(THERMAL_COLUMNS),
        'kinetic_energy_convention': 'includes center-of-mass translation',
        'temperature_convention': 'COM removed; N-1 degrees of freedom; k_B=1',
        'sampling': {'initial_state': True, 'final_state': True,
                     'regular_interval': config.output_every * config.dt},
        'software': {'python': platform.python_version(), 'numpy': np.__version__},
    }
    return Trajectory(X, V, F, thermal, XU, metadata)


def main():
    parser = argparse.ArgumentParser(description='Compute 1D LJ trajectories (no plotting)')
    parser.add_argument('--boundary', choices=['free', 'periodic', 'both'], default='both')
    parser.add_argument('--out', type=Path, default=Path(__file__).resolve().parent,
                        help='Output root; data are written to ROOT/results/BOUNDARY')
    args = parser.parse_args()
    targets = ['free', 'periodic'] if args.boundary == 'both' else [args.boundary]
    for boundary in targets:
        trajectory = simulate(boundary)
        directory = args.out / 'results' / boundary
        save_trajectory(trajectory, directory)
        drift = energy_drift(trajectory)
        print(f"{boundary:8s}: E(0)={drift['initial_energy']:.8f}, "
              f"max |E(t)-E(0)|={drift['max_absolute_error']:.3e}, "
              f"relative={drift['relative_error']:.3e}; saved {directory}")


if __name__ == '__main__':
    main()
