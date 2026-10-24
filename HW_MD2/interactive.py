"""Export a fully offline, linked dashboard without importing the MD solver."""
import json
from pathlib import Path

import numpy as np

from analysis import time_edges, view_fields


def dashboard_data(trajectories):
    if not trajectories:
        raise ValueError('At least one trajectory is required')
    cases = []
    for trajectory in trajectories:
        trajectory.validate()
        fields = view_fields(trajectory)
        cases.append({
            'boundary': trajectory.boundary,
            'time': trajectory.time.tolist(),
            'time_edges': time_edges(trajectory.time).tolist(),
            'positions': trajectory.positions.tolist(),
            'positions_unwrapped': trajectory.positions_unwrapped.tolist(),
            'thermal': trajectory.thermal.tolist(),
            'fields': {name: values.tolist() for name, values in fields.items()},
            'metadata': trajectory.metadata,
        })
    start = max(item.time[0] for item in trajectories)
    end = min(item.time[-1] for item in trajectories)
    if end <= start:
        raise ValueError('Trajectories must have an overlapping time interval')
    # Use physical time, not frame number, to link different saved sampling grids.
    times = np.unique(np.concatenate([item.time for item in trajectories]))
    times = times[(times >= start) & (times <= end)]
    limits = {key: max(1e-15, max(float(np.max(np.abs(view_fields(item)[key])))
                                for item in trajectories))
              for key in ('displacement', 'velocity', 'kinetic')}
    return {'cases': cases, 'time': times.tolist(), 'limits': limits}


def write_dashboard(trajectories, path):
    # Import Plotly only for HTML export, keeping PNG-only workflows independent.
    from plotly.offline import get_plotlyjs

    assets = Path(__file__).resolve().parent / 'dashboard'
    template = (assets / 'index.html').read_text(encoding='utf-8')
    payload = json.dumps(dashboard_data(trajectories), ensure_ascii=False, allow_nan=False,
                         separators=(',', ':')).replace('<', '\\u003c')
    html = template.replace('__DASHBOARD_STYLE__', (assets / 'style.css').read_text(encoding='utf-8'))
    html = html.replace('__PLOTLY_BUNDLE__', get_plotlyjs())
    html = html.replace('__TRAJECTORY_DATA__', payload)
    html = html.replace('__DASHBOARD_SCRIPT__', (assets / 'app.js').read_text(encoding='utf-8'))
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(html, encoding='utf-8')
    return path
