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
mode=${1:?mode}
case_name=${2:?case name}
morphology=${3:?morphology}
seed=${4:?seed}
mkdir -p "runs/$case_name"
cd "runs/$case_name"
sha256sum ../../CuYM.eam.alloy ../../in.prepare ../../in.load
if [[ "$mode" == smoke || "$mode" == prepare ]]; then
    equil_steps=30000
    [[ "$mode" != smoke ]] || equil_steps=400
    srun -n 1 lmp -in ../../in.prepare -log equil.log -var morphology "$morphology" -var seed "$seed" -var equil_steps "$equil_steps"
fi
if [[ "$mode" == smoke || "$mode" == load ]]; then
    load_steps=24000
    [[ "$mode" != smoke ]] || load_steps=400
    srun -n 1 lmp -in ../../in.load -log loading.log -var rate_ps 0.005 -var load_steps "$load_steps" -var dump_every 100
fi
touch "${mode}.complete"
