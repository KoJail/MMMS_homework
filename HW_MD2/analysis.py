"""Read-only numerical analysis of saved trajectories; no integration or plotting."""
import numpy as np


def com_removed_displacement(trajectory):
    """Displacement relative to the initial state, with COM translation removed."""
    displacement = trajectory.positions_unwrapped - trajectory.positions_unwrapped[0]
    return displacement - displacement.mean(axis=1, keepdims=True)


def com_removed_velocity(trajectory):
    return trajectory.velocities - trajectory.velocities.mean(axis=1, keepdims=True)


def particle_kinetic_energy(trajectory):
    """Per-particle internal kinetic energy; not local total energy or heat flux."""
    mass = trajectory.metadata['config']['mass']
    return 0.5 * mass * com_removed_velocity(trajectory) ** 2


def view_fields(trajectory):
    """Physical fields shared by static and interactive views, in archive units."""
    return {'displacement': com_removed_displacement(trajectory),
            'velocity': com_removed_velocity(trajectory),
            'kinetic': particle_kinetic_energy(trajectory)}


def time_edges(time):
    """Cell boundaries for real sampling times, including an uneven final step."""
    midpoints = 0.5 * (time[:-1] + time[1:])
    return np.r_[time[0], midpoints, time[-1]]


def energy_drift(trajectory):
    energy = trajectory.thermal[:, 4]
    absolute = float(np.max(np.abs(energy - energy[0])))
    return {'initial_energy': float(energy[0]), 'max_absolute_error': absolute,
            'relative_error': absolute / max(abs(float(energy[0])), 1e-12)}


def report_metrics(trajectory):
    """Measured finite-run diagnostics; no assumption of thermal equilibrium."""
    config = trajectory.metadata['config']
    mass, n = config['mass'], config['n']
    time = trajectory.time
    momentum = mass * trajectory.velocities.sum(axis=1)
    center = trajectory.positions_unwrapped.mean(axis=1)
    predicted = momentum[0] / (mass * n) * (time - time[0])
    temperature = trajectory.thermal[:, 5]
    return {
        'boundary': trajectory.boundary,
        'initial_total_energy': float(trajectory.thermal[0, 4]),
        'max_energy_error': energy_drift(trajectory)['max_absolute_error'],
        'initial_momentum': float(momentum[0]),
        'max_momentum_error': float(np.max(np.abs(momentum - momentum[0]))),
        'final_com_displacement': float(center[-1] - center[0]),
        'max_com_ballistic_error': float(np.max(np.abs(center - center[0] - predicted))),
        'initial_max_force': float(np.max(np.abs(trajectory.forces[0]))),
        'initial_internal_kinetic_energy': float(particle_kinetic_energy(trajectory)[0].sum()),
        'max_internal_displacement': float(np.max(np.abs(com_removed_displacement(trajectory)))),
        'initial_temperature': float(temperature[0]),
        'min_temperature': float(temperature.min()),
        'max_temperature': float(temperature.max()),
        'time_average_temperature': float(np.trapezoid(temperature, time) / (time[-1] - time[0])),
    }
