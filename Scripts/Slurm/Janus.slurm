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
module load openmpi/openmpi-1.6.4_intel-13.0.0_ib
module load hdf5/hdf5-1.8.11_szip-2.1_zlib-1.2.78_jpeglib-8d_openmpi-1.6.4_intel-13.0.0_ib
module load fftw/fftw-3.3.3_openmpi-1.6.4_intel-13.0.0_double_ib

#module load python/anaconda-2.0.0
#module load gcc/gcc-4.7.2
#module load openmpi/openmpi-1.6.4_gcc-4.7.2_ib
#module load hdf5/hdf5-1.8.11_szip-2.1_zlib-1.2.78_jpeglib-8d_openmpi-1.6.4_gcc-4.7.2_ib
#module load fftw/fftw-3.3.3_openmpi-1.6.4_gcc-4.7.2_double_ib

mpirun ./BoussinesqFPlane3DQGModel
