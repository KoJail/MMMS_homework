"""Render actual simulation structures using the official OVITO Python package.

From the workspace root:
    .venv/bin/python HW_MD3/scripts/render_ovito.py
No GUI interactions and no format autodetection are required.
"""
import os
os.environ.setdefault('OVITO_THREAD_COUNT', '4')

from pathlib import Path
import argparse
import json
import numpy as np
import ovito
from ovito.io import import_file
from ovito.modifiers import CommonNeighborAnalysisModifier, ColorCodingModifier
from ovito.vis import Viewport, TachyonRenderer, TextLabelOverlay, ColorLegendOverlay
from ovito.qt_compat import QtCore

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'figures/ovito'
STATES = ROOT / 'ovito/states'
OUT.mkdir(parents=True, exist_ok=True)
STATES.mkdir(parents=True, exist_ok=True)
SIZE = (1600, 1000)

def render(hkl, kind, mode):
    name = f'Cu{hkl}_{kind}_{mode}'
    path = ROOT / f'results/cluster/inputs/{hkl}_n12_v40/{kind}.dump'
    pipeline = import_file(str(path), input_format='lammps/dump')
    initial = pipeline.compute()
    n = initial.particles.count
    if n <= 0 or 'c_pea' not in initial.particles:
        raise ValueError(f'Invalid particle input: {path}')
    pipeline.source.data.particles.vis.radius = .9
    pipeline.source.data.cell.vis.enabled = True
    pipeline.source.data.cell.vis.line_width = .065
    pipeline.source.data.cell.vis.render_cell = True
    cna = CommonNeighborAnalysisModifier(mode=CommonNeighborAnalysisModifier.Mode.AdaptiveCutoff)
    cna.structures[CommonNeighborAnalysisModifier.Type.FCC].color = (.24, .65, .46)
    cna.structures[CommonNeighborAnalysisModifier.Type.OTHER].color = (.72, .74, .77)
    pipeline.modifiers.append(cna)
    data = pipeline.compute()
    counts = {key: int(value) for key,value in data.attributes.items() if key.startswith('CommonNeighborAnalysis.counts.')}
    coloring = None
    if mode == 'energy':
        coloring = ColorCodingModifier(property='c_pea', start_value=-3.55, end_value=-2.8,
                                       gradient=ColorCodingModifier.Viridis())
        pipeline.modifiers.append(coloring)
    pipeline.add_to_scene()
    vp = ovito.scene.viewports.active_vp
    vp.type = Viewport.Type.Ortho
    vp.camera_dir = (.12, -.28, -1) if kind == 'slab' else (-.8, -.65, -1)
    vp.camera_up = (0, 1, 0)
    vp.overlays.clear()
    vp.zoom_all(size=SIZE)
    vp.fov *= 1.17
    title = f'Cu({hkl}) | {"Periodic bulk" if kind == "bulk" else "Relaxed slab + vacuum"}'
    vp.overlays.append(TextLabelOverlay(text=title, font_size=.033,
        alignment=QtCore.Qt.AlignmentFlag.AlignLeft | QtCore.Qt.AlignmentFlag.AlignTop,
        offset_x=.03, offset_y=-.025, text_color=(.1,.15,.2)))
    if mode == 'cna':
        description = f'Adaptive CNA: FCC (green) {counts.get("CommonNeighborAnalysis.counts.FCC",0)} / {n}; Other (gray) {counts.get("CommonNeighborAnalysis.counts.OTHER",0)}'
        vp.overlays.append(TextLabelOverlay(text=description, font_size=.024,
            alignment=QtCore.Qt.AlignmentFlag.AlignLeft | QtCore.Qt.AlignmentFlag.AlignBottom,
            offset_x=.03, offset_y=.03, text_color=(.1,.15,.2)))
    else:
        vp.overlays.append(ColorLegendOverlay(modifier=coloring, caption='Potential energy (eV/atom)',
            alignment=QtCore.Qt.AlignmentFlag.AlignHCenter | QtCore.Qt.AlignmentFlag.AlignBottom,
            orientation=QtCore.Qt.Orientation.Horizontal, offset_y=.04, font_size=.06,
            legend_size=.38, aspect_ratio=18, text_color=(.1,.15,.2), format_string='%.2f'))
    output = OUT / f'{name}.png'
    vp.render_image(filename=str(output), size=SIZE, background=(1,1,1),
                    renderer=TachyonRenderer(ambient_occlusion=True, shadows=True,
                                            antialiasing_samples=4, ambient_occlusion_samples=12))
    # Explicitly saves the chosen importer, modifiers, and fitted camera.
    state = STATES / f'{name}.ovito'
    ovito.scene.save(str(state))
    pipeline.remove_from_scene()
    record = dict(name=name, source=str(path.relative_to(ROOT)), atoms=n, CNA_counts=counts,
                  energy_min_eV=float(np.min(data.particles['c_pea'])),
                  energy_max_eV=float(np.max(data.particles['c_pea'])),
                  image=str(output.relative_to(ROOT)), state=str(state.relative_to(ROOT)))
    print(json.dumps(record), flush=True)
    return record

if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--preview', action='store_true', help='Render the bulk CNA image only')
    parser.add_argument('--appendix', action='store_true', help='Render CNA and energy images for Cu(110) and Cu(111)')
    args=parser.parse_args()
    if args.appendix:
        jobs=[(hkl,'slab',mode) for hkl in ('110','111') for mode in ('cna','energy')]
    else:
        jobs=[('100','bulk','cna')]
    if not args.preview and not args.appendix:
        jobs += [('100','slab','cna'),('100','slab','energy'),('110','slab','energy'),('111','slab','energy')]
    records=[render(*job) for job in jobs]
    manifest='appendix_manifest.json' if args.appendix else ('preview.json' if args.preview else 'render_manifest.json')
    (OUT / manifest).write_text(
        json.dumps(dict(ovito_version=ovito.version_string, renderer='TachyonRenderer', size=SIZE,
                        input_format='lammps/dump', energy_range_eV=[-3.55,-2.8], images=records), indent=2)+'\n')
