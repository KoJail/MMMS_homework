"""Static report evidence from existing archives; no solver or new experiments."""
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from analysis import com_removed_displacement, com_removed_velocity, report_metrics

COLORS = {'free': '#2868a5', 'periodic': '#b96026'}
PROBE_COLORS = ['#c45b27', '#2868a5', '#37836b', '#795796']


def _save(fig, directory, name):
    fig.savefig(directory / f'{name}.png', dpi=300)
    plt.close(fig)


def _finish(ax, xlabel, ylabel):
    ax.set(xlabel=xlabel, ylabel=ylabel)
    ax.grid(alpha=.18)
    ax.spines[['top', 'right']].set_visible(False)


def _initial_conditions(trajectories, directory):
    fig = plt.figure(figsize=(9.2, max(4.5, 2.3 * len(trajectories))), layout='constrained')
    grid = fig.add_gridspec(len(trajectories), 2, width_ratios=[1.5, 1])
    force_ax = fig.add_subplot(grid[:, 1])
    for row, trajectory in enumerate(trajectories):
        ax = fig.add_subplot(grid[row, 0])
        meta = trajectory.metadata
        x = trajectory.positions_unwrapped[0]
        kick = meta['initialization']['kick_particle_index_zero_based']
        spacing, length = meta['spacing'], meta['box_length']
        color = COLORS[trajectory.boundary]
        ax.axvspan(x[kick] - meta['rcut'], x[kick] + meta['rcut'], color=color, alpha=.08)
        ax.text(.02, .94, 'Shading: cutoff range of kicked particle',
                transform=ax.transAxes, va='top', fontsize=8)
        ax.plot(x, np.zeros_like(x), '-', color='#c1c8d1', lw=1, zorder=1)
        ax.scatter(x, np.zeros_like(x), s=65, color=color, zorder=3)
        for j, position in enumerate(x):
            ax.text(position, -.12, str(j + 1), ha='center', va='top', fontsize=9)
        ax.annotate('', xy=(x[kick] + .65 * spacing, .17), xytext=(x[kick], .17),
                    arrowprops={'arrowstyle': '->', 'color': '#273d59', 'lw': 1.5})
        ax.text(x[kick], .28, f'Kick: v = {trajectory.velocities[0, kick]:g}', fontsize=9)
        if trajectory.boundary == 'periodic':
            ghosts = [x[-1] - length, x[0] + length]
            ax.scatter(ghosts, [0, 0], s=65, facecolors='none', edgecolors=color)
            ax.plot([ghosts[0], x[0]], [0, 0], '--', color=color, lw=1)
            ax.plot([x[-1], ghosts[1]], [0, 0], '--', color=color, lw=1)
            for position, label in zip(ghosts, [f'{len(x)} image', '1 image']):
                ax.text(position, -.12, label, ha='center', va='top', fontsize=8)
        ax.set(xlim=(-1.7 * spacing, length + .9 * spacing), ylim=(-.45, .65),
               yticks=[], xlabel='Initial x (reduced)', title=f'{trajectory.boundary.capitalize()}: IDs and initial kick')
        ax.spines[['left', 'top', 'right']].set_visible(False)
        force_ax.plot(np.arange(1, len(x) + 1), trajectory.forces[0],
                      marker='o' if trajectory.boundary == 'free' else 's',
                      color=color, label=trajectory.boundary.capitalize(), lw=1.3, ms=5)
    _finish(force_ax, 'Particle ID', 'Initial force (reduced)')
    force_ax.set_title('Initial force balance')
    force_ax.axhline(0, color='0.5', lw=.7)
    force_ax.legend()
    fig.suptitle('Initial conditions: a common kick, different boundary force balance', fontsize=14)
    _save(fig, directory, 'initial_conditions')


def _conservation(trajectories, directory):
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 6.5), layout='constrained')
    for trajectory in trajectories:
        color, label = COLORS[trajectory.boundary], trajectory.boundary.capitalize()
        t, thermal = trajectory.time, trajectory.thermal
        config = trajectory.metadata['config']
        p = config['mass'] * trajectory.velocities.sum(axis=1)
        com = trajectory.positions_unwrapped.mean(axis=1)
        drift = com - com[0]
        error = (thermal[:, 4] - thermal[0, 4]) * 1e6
        axes[0, 0].plot(t, error, color=color, lw=.6, alpha=.22)
        axes[0, 0].plot(t, np.maximum.accumulate(np.abs(error)), color=color,
                        label=f'{label}: max {np.max(np.abs(error)):.3f}', lw=1.3)
        axes[0, 1].plot(t, (p - p[0]) * 1e15, color=color, label=label, lw=1)
        axes[1, 0].plot(t, drift, color=color, label=label, linestyle='none',
                       marker='o' if trajectory.boundary == 'free' else '^',
                       markevery=max(1, len(t) // 8), ms=4)
        axes[1, 0].plot(t, p[0] / (config['mass'] * config['n']) * (t - t[0]),
                       '--', color=color, lw=.9, alpha=.65)
        u = com_removed_displacement(trajectory)
        axes[1, 1].plot(t, np.sqrt(np.mean(u ** 2, axis=1)), color=color, label=label, lw=1)
    descriptions = [
        ('A. Energy error: faint = signed; bold = running max', r'Energy error $/10^{-6}$'),
        ('B. Total-momentum error', r'$[P(t)-P(0)]/10^{-15}$'),
        ('C. COM: markers = data; dashed = prediction', 'COM displacement (reduced)'),
        ('D. Internal motion after removing COM translation', 'RMS displacement (reduced)'),
    ]
    for ax, (title, ylabel) in zip(axes.flat, descriptions):
        _finish(ax, 'Time (reduced)', ylabel)
        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=9)
    fig.suptitle('Numerical conservation and separation of collective / internal motion', fontsize=14)
    _save(fig, directory, 'conservation')


def _thermodynamics(trajectories, directory):
    fig, axes = plt.subplots(len(trajectories), 2, squeeze=False,
                             figsize=(9.2, 3.3 * len(trajectories)), layout='constrained')
    energy_limit = max(1e-15, max(float(np.max(np.abs(item.thermal[:, 2:4] - item.thermal[0, 2:4])))
                                for item in trajectories))
    temperature_limit = max(1e-15, max(float(item.thermal[:, 5].max()) for item in trajectories))
    for row, trajectory in enumerate(trajectories):
        t, thermal = trajectory.time, trajectory.thermal
        ax = axes[row, 0]
        ax.plot(t, thermal[:, 2] - thermal[0, 2], label=r'$U(t)-U(0)$', color='#2868a5', lw=.9)
        ax.plot(t, thermal[:, 3] - thermal[0, 3], label=r'$K(t)-K(0)$', color='#b96026', lw=.9)
        ax.set(ylim=(-1.1 * energy_limit, 1.1 * energy_limit), title=f'{trajectory.boundary.capitalize()}: energy exchange')
        _finish(ax, 'Time (reduced)', 'Energy change (reduced)')
        ax.legend(fontsize=9)
        ax = axes[row, 1]
        mean = report_metrics(trajectory)['time_average_temperature']
        ax.plot(t, thermal[:, 5], color=COLORS[trajectory.boundary], lw=.9, label='Instantaneous T*')
        ax.axhline(mean, color='#4e5867', ls='--', lw=1, label=f'Finite-run mean = {mean:.5f}')
        ax.set(ylim=(0, 1.1 * temperature_limit), title=f'{trajectory.boundary.capitalize()}: kinetic temperature')
        _finish(ax, 'Time (reduced)', 'Temperature T*')
        ax.legend(fontsize=9)
    fig.suptitle('NVE energy exchange and fluctuating kinetic temperature', fontsize=14)
    _save(fig, directory, 'thermodynamics')


def _snapshots(trajectories, directory):
    horizon = min(2.0, min(item.time[-1] for item in trajectories))
    targets = np.array([.1, .3, .5, .7]) * horizon
    saved = []
    for trajectory in trajectories:
        velocity = com_removed_velocity(trajectory)
        samples = [int(np.argmin(np.abs(trajectory.time - target))) for target in targets]
        saved.append((trajectory, velocity, samples))
    vmax = max(1e-15, max(float(np.abs(velocity[samples]).max()) for _, velocity, samples in saved))
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 6.2), layout='constrained')
    for column, (ax, target) in enumerate(zip(axes.flat, targets)):
        for trajectory, velocity, samples in saved:
            k = samples[column]
            ids = np.arange(1, velocity.shape[1] + 1)
            label = f'{trajectory.boundary.capitalize()}, t={trajectory.time[k]:.3f}'
            ax.plot(ids, velocity[k], 'o-' if trajectory.boundary == 'free' else 's--',
                    color=COLORS[trajectory.boundary], lw=1.2, ms=4, label=label)
        ax.axhline(0, color='0.5', lw=.7)
        ax.set(ylim=(-1.1 * vmax, 1.1 * vmax), title=f'Target time {target:.3f}')
        _finish(ax, 'Particle ID', "COM-removed velocity v'")
        ax.legend(fontsize=9)
    fig.suptitle('Saved velocity profiles: redistribution of the initial disturbance', fontsize=14)
    _save(fig, directory, 'wave_snapshots')


def _responses(trajectories, directory):
    end = min(2.0, min(item.time[-1] for item in trajectories))
    fig, axes = plt.subplots(len(trajectories), 2, squeeze=False,
                             figsize=(9.2, 3.3 * len(trajectories)), layout='constrained')
    umax = max(1e-15, max(float(np.abs(com_removed_displacement(item)).max()) for item in trajectories))
    vmax = max(1e-15, max(float(np.abs(com_removed_velocity(item)).max()) for item in trajectories))
    for row, trajectory in enumerate(trajectories):
        n = trajectory.positions.shape[1]
        kick = trajectory.metadata['initialization']['kick_particle_index_zero_based']
        probes = list(dict.fromkeys([kick, min(kick + 1, n - 1), 0, n - 1]))
        for column, (values, limit, ylabel) in enumerate([
            (com_removed_displacement(trajectory), umax, 'COM-removed displacement u'),
            (com_removed_velocity(trajectory), vmax, "COM-removed velocity v'")]):
            ax = axes[row, column]
            for color, particle in zip(PROBE_COLORS, probes):
                label = f'ID {particle + 1}' + (' (kicked)' if particle == kick else '')
                ax.plot(trajectory.time, values[:, particle], color=color, label=label, lw=1.1)
            ax.set(xlim=(0, end), ylim=(-1.1 * limit, 1.1 * limit), title=trajectory.boundary.capitalize())
            _finish(ax, 'Time (reduced)', ylabel)
            ax.legend(ncol=2, fontsize=8)
    fig.suptitle('Representative particle responses under the two boundaries', fontsize=14)
    _save(fig, directory, 'particle_responses')


def write_report_figures(trajectories, directory):
    """Create complementary evidence figures and a machine-readable measurement table."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    with plt.rc_context({'font.size': 10, 'axes.titlesize': 11, 'axes.labelsize': 10}):
        for plot in (_initial_conditions, _conservation, _thermodynamics, _snapshots, _responses):
            plot(trajectories, directory)
    metrics = [report_metrics(item) for item in trajectories]
    with (directory / 'report_metrics.csv').open('w', newline='', encoding='utf-8') as fp:
        writer = csv.DictWriter(fp, fieldnames=list(metrics[0]))
        writer.writeheader()
        writer.writerows(metrics)
