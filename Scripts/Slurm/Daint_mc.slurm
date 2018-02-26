#!/bin/bash -l
#
#SBATCH --job-name="DynShell"
#SBATCH --partition=normal
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --ntasks-per-core=1
#SBATCH --constraint=mc
#SBATCH --output=slurm_out.%j
#SBATCH --error=slurm_err.%j

module load daint-mc
module switch PrgEnv-cray PrgEnv-gnu
module load fftw
module load cray-hdf5-parallel
module load cray-python
module load cray-tpsl
module load CMake

srun -n $SLURM_NTASKS --ntasks-per-node $SLURM_NTASKS_PER_NODE ./BoussinesqDynamoShellStdModel
