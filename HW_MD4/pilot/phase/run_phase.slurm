#!/bin/bash
#SBATCH --partition=64c512g
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:25:00

set -eo pipefail
source /etc/profile.d/modules.sh
module use /dssg/share/spack/modules/icelake/linux-centos8-icelake /dssg/share/modulefiles
module load lammps/20240829.1-intel-2021.4.0
set -u
export OMP_NUM_THREADS=1
case_name=${1:?reference case}
delay=${2:?phase delay in timesteps}
root=$PWD
reference="$root/runs/$case_name/equilibrated.restart"
test -s "$reference"
case_dir="runs/phase-${case_name}-${delay}"
mkdir -p "$case_dir"
cd "$case_dir"
test ! -e phase.complete
ln -sfn "$reference" reference.restart
sha256sum "$reference" ../../CuYM.eam.alloy ../../in.phase ../../in.load
srun -n 1 lmp -in ../../in.phase -log hold.log -var delay_steps "$delay"
srun -n 1 lmp -in ../../in.load -log loading.log -var rate_ps 0.005 -var load_steps 18000 -var dump_every 20
touch phase.complete
