"""Parse actual cluster results, validate bookkeeping, and construct Wulff facets.

Uses NumPy and Matplotlib only. Plane-triplet enumeration avoids extra packages.
"""
from pathlib import Path
from itertools import combinations, permutations, product
import csv
import json
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / 'results/cluster/inputs'
FIG = ROOT / 'figures'
FIG.mkdir(exist_ok=True)
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                     'axes.spines.top': False, 'axes.spines.right': False})
COLORS = {'100': '#e5a137', '110': '#57aa86', '111': '#4f7fbc'}
CONV = 16.02176634
NIST = {'100': 1.34525, '110': 1.47539, '111': 1.23946}

def read_dump(path):
    lines = path.read_text().splitlines()
    if lines.count('ITEM: TIMESTEP') != 1:
        raise ValueError(f'Expected exactly one frame in {path}')
    n = int(lines[3])
    bounds = np.array([[float(x) for x in l.split()] for l in lines[5:8]])
    fields = lines[8].split()[2:]
    data = np.array([[float(x) for x in l.split()] for l in lines[9:]])
    if len(data) != n or len(set(data[:, fields.index('id')])) != n:
        raise ValueError(f'Atom count/id mismatch: {path}')
    return bounds, {f: data[:, i] for i, f in enumerate(fields)}

def stops(log):
    return re.findall(r'Stopping criterion = (.+)', log)

def result_records(log):
    records = {}
    for line in log.splitlines():
        if line.startswith('RESULT '):
            parts = line.split()
            records[parts[1]] = {k: float(v) for k, v in (x.split('=') for x in parts[2:])}
    if set(records) != {'bulk', 'unrelaxed', 'slab'}:
        raise ValueError('Missing or unexpected RESULT records')
    return records

rows = []
for case in json.loads((ROOT / 'inputs/cases.json').read_text()):
    if case['kind'] != 'production':
        continue
    folder = RAW / case['name']
    log = (folder / 'log.lammps').read_text()
    if 'ERROR' in log or 'Total wall time:' not in log:
        raise ValueError(f'Failed/incomplete case: {folder}')
    records = result_records(log)
    b, u, s = (records[k] for k in ('bulk', 'unrelaxed', 'slab'))
    if not b['N'] == u['N'] == s['N']:
        raise ValueError('A/B atom numbers differ')
    area = b['Ly'] * b['Lz']
    if not np.allclose([b['Ly'], b['Lz']], [s['Ly'], s['Lz']], atol=1e-10, rtol=0):
        raise ValueError('Inconsistent area')
    bounds, atoms = read_dump(folder / 'slab.dump')
    energy_error = abs(atoms['c_pea'].sum() - s['E'])
    if energy_error > 1e-7:
        raise ValueError('Per-atom energies do not sum to total potential energy')
    if s['E'] > u['E'] + 1e-7:
        raise ValueError('Relaxation increased the energy')
    gap = s['Lx'] - np.ptp(atoms['x'])
    if gap < 5.50679:
        raise ValueError('Vacuum shorter than interaction cutoff')
    row = {**case, 'N': int(b['N']), 'thickness_box_A': b['Lx'],
           'Ly_A': b['Ly'], 'Lz_A': b['Lz'], 'area_A2': area,
           'E_bulk_eV': b['E'], 'E_unrelaxed_eV': u['E'], 'E_slab_eV': s['E'],
           'e_bulk_eV': b['E']/b['N'], 'bulk_pressure_bar': b['Pressure'],
           'gamma_J_m2': (s['E']-b['E'])/(2*area)*CONV,
           'gamma_unrelaxed_J_m2': (u['E']-b['E'])/(2*area)*CONV,
           'fmax_eV_A': s['Fmax'], 'stopping': stops(log)[-1],
           'dump_energy_error_eV': energy_error, 'geometric_gap_A': gap}
    rows.append(row)
with (ROOT / 'results/surface_energies.csv').open('w') as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)

BASE = {r['hkl']: r for r in rows if r['nx'] == 12 and r['vacuum'] == 40 and 'tight' not in r['name']}
if np.ptp([r['e_bulk_eV'] for r in rows]) > 1e-9:
    raise ValueError('Bulk reference energy inconsistent across orientations')

def wulff(energies):
    normals, distances, families = [], [], []
    for family, energy in energies.items():
        ints = set(permutations(map(int, family)))
        signed = sorted({tuple(np.array(p)*sign) for p in ints for sign in product((-1, 1), repeat=3)})
        for p in signed:
            p = np.array(p, dtype=float)
            normals.append(p / np.linalg.norm(p))
            distances.append(energy)
            families.append(family)
    normals, distances = np.array(normals), np.array(distances)
    vertices = []
    for ix in combinations(range(len(normals)), 3):
        mat = normals[list(ix)]
        if abs(np.linalg.det(mat)) < 1e-10:
            continue
        v = np.linalg.solve(mat, distances[list(ix)])
        if np.all(normals @ v <= distances + 1e-9) and not any(np.linalg.norm(v-w) < 1e-8 for w in vertices):
            vertices.append(v)
    vertices = np.array(vertices)
    facets = []
    for normal, distance, family in zip(normals, distances, families):
        poly = vertices[abs(vertices @ normal - distance) < 1e-8]
        if len(poly) < 3:
            continue
        center = poly.mean(axis=0)
        u = poly[0]-center
        u /= np.linalg.norm(u)
        v = np.cross(normal, u)
        angles = np.arctan2((poly-center)@v, (poly-center)@u)
        poly = poly[np.argsort(angles)]
        area = sum(np.linalg.norm(np.cross(poly[i]-center, poly[(i+1)%len(poly)]-center))/2 for i in range(len(poly)))
        facets.append(dict(family=family, polygon=poly.tolist(), area=float(area), distance=float(distance)))
    areas = {h: sum(f['area'] for f in facets if f['family'] == h) for h in energies}
    total = sum(areas.values())
    # Euler topology and divergence-theorem volume provide independent bookkeeping checks.
    edges = set()
    for facet in facets:
        ids = [int(np.argmin(np.linalg.norm(vertices-np.array(p), axis=1))) for p in facet['polygon']]
        edges.update(tuple(sorted((ids[i], ids[(i+1)%len(ids)]))) for i in range(len(ids)))
    if len(vertices)-len(edges)+len(facets) != 2:
        raise ValueError('Wulff polyhedron fails Euler topology check')
    volume = sum(f['distance']*f['area']/3 for f in facets)
    return dict(energies=energies, area_fractions={h: a/total for h, a in areas.items()},
                face_counts={h: sum(f['family']==h for f in facets) for h in energies},
                vertices=vertices.tolist(), facets=facets, edges=len(edges), scaled_volume=volume)

shape = wulff({h: r['gamma_J_m2'] for h, r in BASE.items()})
(ROOT / 'results/wulff.json').write_text(json.dumps(shape, indent=2)+'\n')

validation = {}
for h, base in BASE.items():
    byname = {r['name']: r for r in rows}
    thick = byname[f'{h}_n18_v40']
    vacuum = byname[f'{h}_n12_v60']
    tight = byname[f'{h}_tight']
    validation[h] = dict(thickness_relative_change=abs(thick['gamma_J_m2']/base['gamma_J_m2']-1),
                         vacuum_relative_change=abs(vacuum['gamma_J_m2']/base['gamma_J_m2']-1),
                         tolerance_relative_change=abs(tight['gamma_J_m2']/base['gamma_J_m2']-1),
                         NIST_J_m2=NIST[h], NIST_relative_difference=base['gamma_J_m2']/NIST[h]-1)
teacher = {}
for h in ('100', '111'):
    folder = RAW / f'teacher_{h}'
    energies = []
    for kind in ('bulk', 'slab'):
        log = (folder / f'{kind}.log').read_text()
        match = re.search(r'Energy initial, next-to-last, final =\s*\n\s*(\S+)\s+(\S+)\s+(\S+)', log)
        if not match:
            raise ValueError(f'Missing teacher minimization energy: {folder}')
        energies.append(float(match[3]))
    data = (folder / f'data.p{h}.original').read_text()
    lengths = {}
    for axis in 'xyz':
        m = re.search(rf'^(\S+)\s+(\S+)\s+{axis}lo {axis}hi', data, re.M)
        lengths[axis] = float(m[2])-float(m[1])
    teacher[h] = dict(E_bulk_eV=energies[0], E_slab_eV=energies[1],
                      gamma_J_m2=(energies[1]-energies[0])/(2*lengths['y']*lengths['z'])*CONV)
summary = dict(baseline=BASE, validation=validation, teacher_reproduction=teacher,
               area_fractions=shape['area_fractions'], face_counts=shape['face_counts'],
               h_ratio_to_111={h: r['gamma_J_m2']/BASE['111']['gamma_J_m2'] for h, r in BASE.items()})
(ROOT / 'results/summary.json').write_text(json.dumps(summary, indent=2)+'\n')

def save(fig, name):
    fig.savefig(FIG / f'{name}.png', dpi=200, bbox_inches='tight')
    fig.savefig(FIG / f'{name}.svg', bbox_inches='tight')
    plt.close(fig)

fig, ax = plt.subplots(figsize=(6.8, 3.8))
hs = ['111', '100', '110']
x = np.arange(3)
ax.bar(x-.17, [BASE[h]['gamma_unrelaxed_J_m2'] for h in hs], .32, color='#c7cdd3', label='Unrelaxed (this work)')
ax.bar(x+.17, [BASE[h]['gamma_J_m2'] for h in hs], .32, color=[COLORS[h] for h in hs], label='Relaxed (this work)')
ax.scatter(x+.17, [NIST[h] for h in hs], marker='_', s=220, color='black', zorder=3, label='NIST EAM1 reference')
for i,h in enumerate(hs):
    ax.text(i+.17, BASE[h]['gamma_J_m2']+.025, f"{BASE[h]['gamma_J_m2']:.4f}", ha='center', fontsize=9)
ax.set(xticks=x, xticklabels=[f'({h})' for h in hs], ylabel=r'Surface energy (J m$^{-2}$)', ylim=(0,1.95))
ax.legend(frameon=False, fontsize=8, loc='upper left')
save(fig, 'surface_energies')

fig, axes = plt.subplots(1, 3, figsize=(10, 3.1))
for ax,h in zip(axes,hs):
    selected = sorted([r for r in rows if r['hkl']==h and r['vacuum']==40 and 'tight' not in r['name']], key=lambda r:r['nx'])
    vals = [r['gamma_J_m2'] for r in selected]
    ax.plot([r['thickness_box_A'] for r in selected], vals, 'o-', color=COLORS[h])
    ax.set(title=f'Cu({h})', xlabel='Bulk-box thickness (A)', ylabel=r'$\gamma$ (J m$^{-2}$)')
    span = max(np.ptp(vals)*1.6, .00004)
    ax.set_ylim(np.mean(vals)-span, np.mean(vals)+span)
    ax.ticklabel_format(axis='y', style='plain', useOffset=False)
fig.tight_layout()
save(fig, 'thickness_convergence')

fig = plt.figure(figsize=(10, 4.5))
ax = fig.add_subplot(121, projection='3d')
ax.add_collection3d(Poly3DCollection([f['polygon'] for f in shape['facets']],
    facecolors=[COLORS[f['family']] for f in shape['facets']], edgecolors='white', linewidths=.8, alpha=1))
lim = max(np.abs(shape['vertices']).ravel())*1.07
ax.set(xlim=(-lim,lim), ylim=(-lim,lim), zlim=(-lim,lim))
ax.set_box_aspect((1,1,1))
ax.view_init(elev=23, azim=32)
ax.set_axis_off()
ax.set_title('Wulff shape: three facet families', pad=0)
ax2 = fig.add_subplot(122)
fractions = [shape['area_fractions'][h]*100 for h in hs]
ax2.barh([f'{{{h}}}' for h in hs], fractions, color=[COLORS[h] for h in hs])
for i,v in enumerate(fractions):
    ax2.text(v+1, i, f'{v:.2f}%', va='center')
ax2.set(xlim=(0,100), xlabel='Total exposed area (%)')
fig.tight_layout()
save(fig, 'wulff_shape')

# Actual atomic coordinates, not an OVITO screenshot. Plot one thin z slice.
fig, axes = plt.subplots(2, 1, figsize=(10, 5.4))
for ax,kind in zip(axes,('bulk','slab')):
    bounds, atoms = read_dump(RAW / '100_n12_v40' / f'{kind}.dump')
    z = atoms['z']
    select = abs(z-np.median(z)) < 1.0
    scat = ax.scatter(atoms['x'][select], atoms['y'][select], c=atoms['c_pea'][select], cmap='viridis', vmin=-3.55,vmax=-2.95, s=24, edgecolors='none')
    ax.set(xlim=bounds[0], ylim=bounds[1], xlabel='x / [100] (A)', ylabel='y / [010] (A)', title=f'Cu(100) {kind}: real-coordinate slice')
    ax.set_aspect('equal')
    if kind=='slab':
        ax.text(-20,11,'Vacuum',ha='center',color='#666666')
fig.tight_layout()
fig.colorbar(scat, ax=axes, label='Potential energy per atom (eV)', fraction=.02, pad=.025)
save(fig, 'bulk_slab')
print(json.dumps(summary, indent=2))
