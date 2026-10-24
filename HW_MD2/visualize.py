"""Plot saved trajectories independently; never run the simulation."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from analysis import time_edges, view_fields
from trajectory_io import load_trajectory

FIELDS = ('displacement', 'velocity', 'kinetic')
LABELS = ('COM-removed displacement', 'COM-removed velocity', 'Internal kinetic energy per particle')
UNITS = ('Displacement (reduced)', 'Velocity (reduced)', 'Energy (reduced)')


def field_limits(trajectories):
    """One fixed color scale per quantity across boundaries and all saved times."""
    limits = {key: 0.0 for key in FIELDS}
    for trajectory in trajectories:
        for key, values in view_fields(trajectory).items():
            limits[key] = max(limits[key], float(np.max(np.abs(values))))
    return {key: max(value, 1e-15) for key, value in limits.items()}


def plot_spacetime(trajectories, directory, early_time=2.0):
    """Comparison figures: full histories and an early-time zoom, no interpolation."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    limits = field_limits(trajectories)
    common_end = min(float(item.time[-1]) for item in trajectories)
    for suffix, end in (('spacetime', common_end), ('spacetime_early', min(early_time, common_end))):
        fig, axes = plt.subplots(len(trajectories), 3, squeeze=False,
                                 figsize=(16, 4 * len(trajectories)), layout='constrained')
        for row, trajectory in enumerate(trajectories):
            fields = view_fields(trajectory)
            ids = np.arange(1, trajectory.positions.shape[1] + 1)
            for col, key in enumerate(FIELDS):
                ax = axes[row, col]
                positive = key == 'kinetic'
                mesh = ax.pcolormesh(time_edges(trajectory.time),
                                     np.arange(len(ids) + 1) + 0.5, fields[key].T,
                                     shading='flat', cmap='magma' if positive else 'RdBu_r',
                                     vmin=0 if positive else -limits[key], vmax=limits[key],
                                     rasterized=True)
                ax.set(xlabel='Time (reduced)', ylabel='Particle ID', xlim=(0, end),
                       title=f'{trajectory.boundary.capitalize()}: {LABELS[col]}')
                if len(ids) <= 16:
                    ax.set_yticks(ids)
                fig.colorbar(mesh, ax=ax, label=UNITS[col], shrink=0.85)
        fig.suptitle('1D LJ chain: propagation and boundary effects', fontsize=17)
        fig.savefig(directory / f'{suffix}.png', dpi=180)
        plt.close(fig)


def main():
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description='Plot existing 1D LJ results without recomputing')
    parser.add_argument('--input', type=Path, default=root / 'results',
                        help='Results root, one boundary directory, or one trajectory.npz')
    parser.add_argument('--out', type=Path, default=root / 'plots', help='Plot directory')
    parser.add_argument('--boundary', choices=['free', 'periodic', 'both'], default='both')
    parser.add_argument('--format', choices=['png', 'html', 'both'], default='both',
                        help='Static figures, offline linked dashboard, or both (default)')
    args = parser.parse_args()
    if args.input.is_file() or (args.input / 'trajectory.npz').is_file():
        paths = [args.input]
    else:
        targets = ['free', 'periodic'] if args.boundary == 'both' else [args.boundary]
        paths = [args.input / boundary / 'trajectory.npz' for boundary in targets]
    # Read all inputs before generating output, so missing/legacy data fail early.
    trajectories = [load_trajectory(path) for path in paths]
    for trajectory in trajectories:
        if args.boundary != 'both' and trajectory.boundary != args.boundary:
            parser.error('--boundary does not match the input trajectory')
    if args.format in ('png', 'both'):
        from report_figures import write_report_figures
        plot_spacetime(trajectories, args.out)
        write_report_figures(trajectories, args.out)
        print(f'Static plots saved to {args.out}')
    if args.format in ('html', 'both'):
        from interactive import write_dashboard
        path = write_dashboard(trajectories, args.out / 'dashboard.html')
        print(f'Offline interactive dashboard saved to {path}')


if __name__ == '__main__':
    main()
