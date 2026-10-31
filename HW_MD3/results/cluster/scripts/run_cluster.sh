#!/bin/bash
#SBATCH --job-name=HW_MD3_Cu
#SBATCH --partition=64c512g
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err
set -eo pipefail
if ! type module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh
fi
module load lammps/20240829.1-intel-2021.4.0
set -u
export OMP_NUM_THREADS=1
date -Is
hostname
module list 2>&1
command -v lmp
# Each case has its own working directory; no teacher file is overwritten.
for hkl in 100 111; do
    cd "$SLURM_SUBMIT_DIR/inputs/teacher_$hkl"
    srun -n 1 lmp -in "in.p$hkl" -log bulk.log > bulk.screen
    cp "data.p$hkl" "data.p$hkl.original"
    # Preserve all atom coordinates and extend only xlo by 40 A.
    awk '/xlo xhi/ {$1=$1-40} {print}' "data.p$hkl.original" > "data.p$hkl"
    srun -n 1 lmp -in "in.s$hkl" -log slab.log > slab.screen
done
for input in "$SLURM_SUBMIT_DIR"/inputs/*/in.surface; do
    cd "$(dirname "$input")"
    srun -n 1 lmp -in in.surface -log log.lammps > screen.txt
done
cd "$SLURM_SUBMIT_DIR"
date -Is
touch COMPLETE
