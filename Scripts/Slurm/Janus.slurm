#!/bin/bash
#SBATCH -J FPr1Ra9
#SBATCH --time=0:15:00
#SBATCH -N 1
#SBATCH --tasks-per-node 12
#SBATCH -o FPr1Ra9-%j.out
#SBATCH --qos=janus-debug
#SBATCH --signal=10@60

module load python/anaconda-2.0.0
module load intel/intel-13.0.0
module load openmpi/openmpi-1.8.2_intel-13.0.0
module load hdf5/hdf5-1.8.13_szip-2.1_zlib-1.2.8_jpeglib-9a_openmpi-1.8.2_intel-13.0.0
module load fftw/fftw-3.3.4_openmpi-1.8.2_intel-13.0.0_double

#module load python/anaconda-2.0.0
#module load gcc/gcc-4.9.1
#module load openmpi/openmpi-1.8.2_gcc-4.9.1
#module load hdf5/hdf5-1.8.13_szip-2.1_zlib-1.2.8_jpeglib-9a_openmpi-1.8.2_gcc-4.9.1
#module load fftw/fftw-3.3.4_openmpi-1.8.2_gcc-4.9.1_double

mpirun ./BoussinesqFPlane3DQGModel
