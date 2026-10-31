"""Generate isolated LAMMPS cases; run from HW_MD3. No simulation runs here."""
from pathlib import Path
import json
import shutil

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'reference/mmmsfall/0Surface'
ORIENT = {
    '100': 'orient x 1 0 0 orient y 0 1 0 orient z 0 0 1',
    '110': 'orient x 1 1 0 orient y -1 1 0 orient z 0 0 1',
    '111': 'orient x 1 1 1 orient y 1 -1 0 orient z 1 1 -2',
}
CASES = []

def build(hkl, nx, vacuum=40, ftol='1e-8', label=None):
    name = label or f'{hkl}_n{nx}_v{vacuum}'
    folder = ROOT / 'inputs' / name
    folder.mkdir(parents=True, exist_ok=True)
    shutil.copy2(SOURCE / 'CuYM.eam.alloy', folder)
    common = '''pair_style eam/alloy
pair_coeff * * CuYM.eam.alloy Cu
neighbor 0.3 bin
neigh_modify every 1 delay 0 check yes
compute pea all pe/atom
thermo 20
thermo_style custom step atoms pe lx ly lz press fmax fnorm
thermo_modify format float %.16g
'''
    text = f'''# Generated from teacher geometry; same a=3.615 A and EAM file.
units metal
atom_style atomic
boundary p p p
lattice fcc 3.615 {ORIENT[hkl]}
region box block 0 {nx} 0 6 0 6
create_box 1 box
create_atoms 1 box
{common}
min_style cg
min_modify dmax 0.1
minimize 0.0 {ftol} 10000 100000
variable EB equal $(pe:%.16g)
print "RESULT bulk N=$(count(all):%.0f) E=$(pe:%.16g) Lx=$(lx:%.16g) Ly=$(ly:%.16g) Lz=$(lz:%.16g) Fmax=$(fmax:%.16g) Pressure=$(press:%.16g)"
dump final all custom 1 bulk.dump id type x y z c_pea
dump_modify final sort id format float %.16g
run 0
undump final
write_data bulk.data
change_box all x delta -{vacuum} 0 units box
run 0
print "RESULT unrelaxed N=$(count(all):%.0f) E=$(pe:%.16g) Lx=$(lx:%.16g) Ly=$(ly:%.16g) Lz=$(lz:%.16g) Fmax=$(fmax:%.16g)"
minimize 0.0 {ftol} 10000 100000
print "RESULT slab N=$(count(all):%.0f) E=$(pe:%.16g) Lx=$(lx:%.16g) Ly=$(ly:%.16g) Lz=$(lz:%.16g) Fmax=$(fmax:%.16g)"
dump final all custom 1 slab.dump id type x y z c_pea
dump_modify final sort id format float %.16g
run 0
undump final
write_data slab.data
'''
    (folder / 'in.surface').write_text(text)
    CASES.append(dict(name=name, hkl=hkl, nx=nx, vacuum=vacuum, ftol=ftol, kind='production'))

for hkl in ORIENT:
    for nx in (6, 12, 18):
        build(hkl, nx)
    build(hkl, 12, vacuum=60)
    build(hkl, 12, ftol='1e-10', label=f'{hkl}_tight')

for hkl in ('100', '111'):
    name = f'teacher_{hkl}'
    folder = ROOT / 'inputs' / name
    folder.mkdir(parents=True, exist_ok=True)
    for filename in ('CuYM.eam.alloy', f'in.p{hkl}', f'in.s{hkl}'):
        shutil.copy2(SOURCE / filename, folder)
    CASES.append(dict(name=name, hkl=hkl, kind='teacher'))

(ROOT / 'inputs/cases.json').write_text(json.dumps(CASES, indent=2) + '\n')
print(f'Generated {len(CASES)} cases')
