#!/bin/bash
#SBATCH -J ShellC0	   # job name
#SBATCH -o ShellC0.o%j       # output and error file name (%j expands to jobID)
#SBATCH -n 16              # total number of mpi tasks requested
#SBATCH -p development     # queue (partition) -- normal, development, etc.
#SBATCH -t 01:00:00        # run time (hh:mm:ss)
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes

#export LD_PRELOAD=/opt/apps/python/epd/7.3.2/lib/libmkl_core.so:/opt/apps/python/epd/7.3.2/lib/libmkl_intel_thread.so
#export LD_PRELOAD=${TACC_MKL_LIB}/libmkl_core.so:${TACC_MKL_LIB}/libmkl_intel_thread.so
#export LD_PRELOAD=${TACC_MKL_LIB}/libmkl_rt.so
export PYTHONPATH=$HOME/local/lib/python2.7/site-packages/:$PYTHONPATH

ibrun ./BoussinesqRTCShellStdModel   # run the MPI executable named a.out
