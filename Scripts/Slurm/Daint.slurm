#!/bin/bash -l
#
#SBATCH --job-name="DynCouette"
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm_out.%j
#SBATCH --error=slurm_err.%j

module switch PrgEnv-cray PrgEnv-gnu
module load Python/3.4.3-CrayGNU-5.2.40 
#module load Python/2.7.11-CrayGNU-2016.03 
#module load Python/3.5.1-CrayGNU-2015.11 
module load cray-hdf5-parallel
module load fftw
module load cray-tpsl
srun -n $SLURM_NTASKS -c $SLURM_CPUS_PER_TASK --cpu_bind=rank ./BoussinesqDynamoCouetteShellStdModel
