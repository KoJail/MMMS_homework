"""Render published atomistic panels and export an editable OVITO scene."""
import argparse
from pathlib import Path
from tempfile import TemporaryDirectory

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from ovito import scene
from ovito.io import export_file, import_file
from ovito.modifiers import DeleteSelectedModifier
from ovito.vis import TachyonRenderer, Viewport

from analyze import pipeline_for

FCC = (0.32, 0.47, 0.57)
SURFACE = (0.58, 0.69, 0.74)
HCP = (0.98, 0.55, 0.10)


def look(frame, data):
    xyz = data.particles.positions
    hcp = np.asarray(data.particles["PTM Structure"]) == 2
    surface = (xyz[:, 2] > 53) | (xyz[:, 2] < 6)
    colors = np.tile(FCC, (len(xyz), 1))
    colors[surface] = SURFACE
    colors[hcp] = HCP
    # Show one crystallographic background slice and all HCP atoms.
    background = (xyz[:, 1] > 23) & (xyz[:, 1] < 26)
    data.particles_.create_property("Color", data=colors)
    data.particles_.create_property("Selection", data=~(background | hcp))


def configure_defects(path):
    pipe = pipeline_for(path, keep_ptm_structure=True)
    pipe.modifiers.append(look)
    pipe.modifiers.append(DeleteSelectedModifier())
    data = pipe.compute()
    data.cell.vis.enabled = False
    data.particles.vis.radius = 0.62
    data.surfaces["dxa-defect-mesh"].vis.enabled = False
    network = data.dislocations.vis
    network.line_width = 1.0
    network.show_line_directions = False
    network.shading = network.Shading.Normal
    return pipe


def viewport(position, direction, fov):
    return Viewport(type=Viewport.Type.Ortho, camera_pos=position,
                    camera_dir=direction, camera_up=(0, 0, 1), fov=fov)


def geometry(out, temp):
    fig = plt.figure(figsize=(9.5, 3.5), layout="constrained")
    grid = fig.add_gridspec(1, 3, width_ratios=(1, 1, 1.05))
    axes = [fig.add_subplot(grid[0, i]) for i in range(3)]
    view = viewport((32, 12, 33), (-1.0, 2.3, -.62), 54)
    for ax, (label, case) in zip(axes[:2], (("(a) Flat (001)", "flat-01"), ("(b) Double-layer step pair", "step-01"))):
        pipe = import_file(str(Path("pilot/data") / case / "equilibrated.dump"))
        data = pipe.compute()
        data.cell.vis.enabled = False
        data.particles.vis.radius = 1.10

        def color_surface(frame, data):
            z = data.particles.positions[:, 2]
            colors = np.tile(FCC, (len(z), 1))
            colors[z > 58.8] = SURFACE
            data.particles_.create_property("Color", data=colors)

        pipe.modifiers.append(color_surface)
        panel(pipe, view, 0, ax, temp)
        ax.set_title(label, loc="left", fontweight="bold", fontsize=11)
    ax = axes[2]
    ax.set_title("(c) Section through step pair", loc="left", fontweight="bold", fontsize=11)
    ax.fill([0, .25, .25, .75, .75, 1, 1, 0],
            [1, 1, .88, .88, 1, 1, 0, 0], color="#dce6eb", ec="#344c5a", lw=1.2)
    for edge in (.25, .75):
        ax.plot([edge, edge], [1.02, 1.12], color="#b05337", lw=1)
        ax.text(edge, 1.16, str(edge), ha="center", fontsize=8, color="#9c4a35")
    ax.annotate("Tension [110]", xy=(1.05, .47), xytext=(.34, .47),
                fontsize=9, va="center", arrowprops=dict(arrowstyle="->", color="#245474"))
    ax.text(.5, -.15, "Lower free surface", ha="center", fontsize=8)
    ax.text(.5, 1.29, "(001) upper surface", ha="center", fontsize=8)
    ax.text(.5, -.32, r"$x/L_x$    ·    schematic", ha="center", fontsize=8)
    ax.set(xlim=(-.08, 1.12), ylim=(-.42, 1.45))
    ax.set_box_aspect(680 / 800)
    ax.set_axis_off()
    fig.savefig(out, dpi=260)
    plt.close(fig)


def panel(pipe, view, frame, ax, temp, *, size=(800, 680)):
    pipe.add_to_scene()
    try:
        image = Path(temp) / "panel.png"
        view.render_image(filename=str(image), frame=frame, size=size,
                          renderer=TachyonRenderer(ambient_occlusion=True, ambient_occlusion_samples=12, shadows=False, antialiasing_samples=8),
                          background=(1, 1, 1))
        ax.imshow(plt.imread(image))
        ax.set_axis_off()
    finally:
        pipe.remove_from_scene()


def defects(out, temp, project):
    rows = [
        ("flat-01", (142, 143, 147, 149), (26, 25, 31), 40),
        ("step-01", (104, 106, 108, 112), (40, 25, 56), 16),
        ("step-02", (98, 99, 101, 113), (40, 25, 56), 16),
        ("phase-step-02-1500", (483, 486, 488, 496), (39, 24, 57), 16),
    ]
    fig, axes = plt.subplots(len(rows), 4, figsize=(11, 9.8), layout="constrained")
    for row_axes, (case, frames, center, fov) in zip(axes, rows):
        path = Path("pilot/data") / case / "loading.dump"
        pipe = configure_defects(path)
        view = viewport(center, (0., 1., 0.), fov)
        sample = .01 if case.startswith("phase") else .05
        for ax, frame in zip(row_axes, frames):
            panel(pipe, view, frame, ax, temp)
            ax.set_title(f"{sample * frame:.2f}%", loc="right", fontsize=11)
        row_axes[0].text(0, 1.04, case.upper(), transform=row_axes[0].transAxes,
                         fontsize=10, fontweight="bold")
        # Orthographic view height = 2*fov. The frame is 800 by 680 px.
        fraction = 10 / (2 * fov * 800 / 680)
        ax = row_axes[-1]
        ax.plot([.08, .08 + fraction], [.08, .08], transform=ax.transAxes,
                color="#1d303a", lw=1.6)
        ax.text(.08, .11, "1 nm", transform=ax.transAxes, fontsize=8, color="#1d303a")
    fig.legend(handles=[
        Line2D([], [], marker="o", color=HCP, linestyle="none", label="PTM HCP"),
        Line2D([], [], color="#26732e", linewidth=2, label="DXA Shockley partials"),
    ], loc="outside lower center", frameon=False, ncol=2)
    fig.savefig(out, dpi=245)
    plt.close(fig)

    pipe = configure_defects(Path("pilot/data/step-01/loading.dump"))
    pipe.add_to_scene()
    frame = 108
    # A native project built from static files uses built-in modifiers only.
    scene_data = project.parent / "ovito-data"
    scene_data.mkdir(exist_ok=True)
    export_file(pipe, str(scene_data / "step01-frame108.xyz"), "xyz", frame=frame,
                columns=["Particle Identifier", "Position.X", "Position.Y", "Position.Z",
                         "Color.R", "Color.G", "Color.B"])
    export_file(pipe, str(scene_data / "step01-frame108.ca"), "ca", frame=frame)
    pipe.remove_from_scene()
    atoms = import_file(str(scene_data / "step01-frame108.xyz"))
    lines = import_file(str(scene_data / "step01-frame108.ca"))
    atoms.compute().particles.vis.radius = .62
    atoms.compute().cell.vis.enabled = False
    dxa = lines.compute()
    dxa.cell.vis.enabled = False
    dxa.dislocations.vis.line_width = 1
    for surface in dxa.surfaces.values():
        surface.vis.enabled = False
    atoms.add_to_scene()
    lines.add_to_scene()
    scene.viewports.active_vp.type = Viewport.Type.Ortho
    scene.viewports.active_vp.camera_pos = (40, 5, 56)
    scene.viewports.active_vp.camera_dir = (0, 1, 0)
    scene.viewports.active_vp.fov = 16
    scene.save(str(project))
    atoms.remove_from_scene()
    lines.remove_from_scene()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geometry", type=Path, required=True)
    parser.add_argument("--defects", type=Path, required=True)
    parser.add_argument("--project", type=Path, required=True)
    args = parser.parse_args()
    with TemporaryDirectory() as temp:
        geometry(args.geometry, temp)
        defects(args.defects, temp, args.project)


if __name__ == "__main__":
    main()
