#!/bin/bash

if [ "$#" -lt 1 ]; then
    echo "Requires at least 1 arguments"
    echo "Usage:"
    echo "benchmark_run.sh benchmark_list.dat"
    exit 1
fi

#set -euf
#set -o pipefail
#set -x

# Store current directory
current_dir=`pwd`
build_dir="Master-Benchmarks"

BENCHMARKS_LIST=$1

while IFS=" " read -u "${file_fd}" -r QUICC_APPROX QUICC_GEOMETRY QUICC_MODEL QUICC_MPIALGO QUICC_TIMESTEPPER QUICC_SPLINALG QUICC_FFT QUICC_THREADS QUICC_MEMORYUSAGE QUICC_PLATFORM QUICC_DIM1D QUICC_DIM2D QUICC_DIM3D QUICC_DT QUICC_TRASH
do
   # Ignore lines starting with #
   if [ "${QUICC_APPROX:0:1}" == "#" ];
   then
      continue
   fi

   if [ ! -z "$QUICC_TRASH" ];
   then
      echo "Extra information is present in benchmark list!"
      exit 1
   fi

   QUICC_EXEC=${QUICC_APPROX}${QUICC_GEOMETRY}${QUICC_MODEL}Model

   runpath=${build_dir}/${QUICC_APPROX}/${QUICC_GEOMETRY}/${QUICC_MODEL}/${QUICC_MPIALGO}/${QUICC_TIMESTEPPER}/${QUICC_SPLINALG}/${QUICC_FFT}/${QUICC_THREADS}/${QUICC_MEMORYUSAGE}/${QUICC_PLATFORM}/${QUICC_DIM1D}/${QUICC_DIM2D}/${QUICC_DIM3D}/${QUICC_DT}/
   cd ${runpath}

   for QUICC_GITHASH in */;
   do
      cd ${QUICC_GITHASH}
      QUICC_NCPU=`grep '<cpus>[1-9]*</cpus>' parameters.cfg | sed -e 's,<cpus>\([1-9]*\)</cpus>,\1,g'`
      echo "Running benchmark in ${runpath}/${QUICC_GITHASH} ..."
      if [ ${QUICC_NCPU} -lt 1 ];
	   then
	      ./${QUICC_EXEC}
	   else
	      mpirun -n ${QUICC_NCPU} ./${QUICC_EXEC}
	   fi
	   cd ..
   done

   cd ${current_dir}
done  {file_fd}<"${BENCHMARKS_LIST}"

cd ${current_dir}
echo "Done"
