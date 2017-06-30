#!/bin/bash

if [ "$#" -lt 4 ]; then
    echo "Requires at least 3 arguments\n Usage benchmark_build.sh /project/dir /benchmark/dir benchmark_list.dat [ncpu] [parallel_make]"
    exit 1
fi

if [ "$#" -eq 3 ]; then
   QUICC_NCPU=1 
else
   QUICC_NCPU=$4 
fi

if [ "$#" -eq 4 ]; then
   CI_QUICC_PARAMAKE=1
else
   CI_QUICC_PARAMAKE=$5
fi

set -euf
set -o pipefail
#set -x

# make sorted list of benchmarks
SORTED_BENCHMARKS=$(mktemp "/tmp/sorted_benchmarks_XXXXX.dat")
sort -b -k 1,9 "${3}" > ${SORTED_BENCHMARKS}

# Store current directory
current_dir=`pwd`
PROJECT_DIR=$1
BENCHMARK_DIR=$2

# Get git hash
cd ${PROJECT_DIR}
QUICC_GITHASH=`git rev-parse --short HEAD`
cd ${current_dir}

echo "Make directory structure ..."
build_dir="Release-Benchmark"
mkdir -p ${build_dir}
cd ${build_dir}

comppath=""
while IFS=" " read -r QUICC_APPROX QUICC_GEOMETRY QUICC_MODEL QUICC_MPIALGO QUICC_TIMESTEPPER QUICC_SPLINALG QUICC_THREADS QUICC_MEMORYUSAGE QUICC_PLATFORM QUICC_DIM1D QUICC_DIM2D QUICC_DIM3D QUICC_DT QUICC_TRASH
do
   # Ignore lines starting with#
   if [ "${QUICC_APPROX:0:1}" == "#" ];
   then
      continue
   fi

   if [ ! -z "$QUICC_TRASH" ];
   then
      echo "Extra information is present in benchmark list!"
      exit 1
   fi

   # Store part of path that forces new compilation
   if [ ! "${comppath}" == "${QUICC_APPROX}/${QUICC_GEOMETRY}/${QUICC_MODEL}/${QUICC_MPIALGO}/${QUICC_TIMESTEPPER}/${QUICC_SPLINALG}/${QUICC_THREADS}/${QUICC_MEMORYUSAGE}/${QUICC_PLATFORM}" ];
   then

      QUICC_EXEC=${QUICC_APPROX}${QUICC_GEOMETRY}${QUICC_MODEL}Model

      echo "CMake configuration ..."
      QUICC_BOUNDARYMETHOD="Tau"
      QUICC_FFT="FFTW"
      QUICC_FFTPLAN="Medium"
      QUICC_PYTHON="python34"
      QUICC_SMARTPTR="TR1"
      if [ "${QUICC_MPIALGO}" = "Tubular" ]
      then
         QUICC_TRANSGROUPER="Transform"
      else
         QUICC_TRANSGROUPER="Equation"
      fi
      cmake -DCMAKE_BUILD_TYPE=Release -DQUICC_PLATFORM=${QUICC_PLATFORM} ${PROJECT_DIR} >& /dev/null || true
      cmake -DQUICC_BOUNDARYMETHOD=${QUICC_BOUNDARYMETHOD} -DQUICC_FFTPLAN=${QUICC_FFT} -DQUICC_FFTPLAN=${QUICC_FFTPLAN} -DQUICC_THREADS=${QUICC_THREADS} -DQUICC_MEMORYUSAGE=${QUICC_MEMORYUSAGE} -DQUICC_MPIALGO=${QUICC_MPIALGO} -DQUICC_PYTHON=${QUICC_PYTHON} -DQUICC_SMARTPTR=${QUICC_SMARTPTR} -DQUICC_SPLINALG=${QUICC_SPLINALG} -DQUICC_TIMESTEPPER=${QUICC_TIMESTEPPER} ${PROJECT_DIR} >& /dev/null || true
      cmake -DQUICC_TRANSGROUPER=${QUICC_TRANSGROUPER} ${PROJECT_DIR}

      # Clean buid
      make clean

      echo "Installing Python files ..."
      make install > /dev/null

      echo "Building ${QUICC_EXEC} ..."
      make -j ${CI_QUICC_PARAMAKE} ${QUICC_EXEC} > /dev/null
   fi

   runpath=${QUICC_APPROX}/${QUICC_GEOMETRY}/${QUICC_MODEL}/${QUICC_MPIALGO}/${QUICC_TIMESTEPPER}/${QUICC_SPLINALG}/${QUICC_THREADS}/${QUICC_MEMORYUSAGE}/${QUICC_PLATFORM}/${QUICC_DIM1D}/${QUICC_DIM2D}/${QUICC_DIM3D}/${QUICC_DT}/${QUICC_GITHASH}/
   if [ -d ${runpath} ];
   then
      echo "Destination directory already exists. Aborting"
      exit 1
   fi

   # Create run directory and install data
   echo "Installing benchmark into ${runpath} ..."
   mkdir -p ${runpath}
   cp Executables/${QUICC_EXEC} ${runpath}
   BENCHMARK_DATA_DIR=${BENCHMARK_DIR}/${QUICC_APPROX}/${QUICC_GEOMETRY}/${QUICC_MODEL}/setup/Master
   cp ${BENCHMARK_DATA_DIR}/state_initial.hdf5 ${runpath}/
   sed -e "s/++DIM1D++/${QUICC_DIM1D}/g" -e "s/++DIM2D++/${QUICC_DIM2D}/g" -e "s/++DIM3D++/${QUICC_DIM3D}/g" -e "s/++NCPU++/${QUICC_NCPU}/g"  -e "s/++DT++/${QUICC_DT}/g" ${BENCHMARK_DATA_DIR}/parameters_template.cfg > ${runpath}/parameters.cfg

   comppath=${QUICC_APPROX}/${QUICC_GEOMETRY}/${QUICC_MODEL}/${QUICC_MPIALGO}/${QUICC_TIMESTEPPER}/${QUICC_SPLINALG}/${QUICC_THREADS}/${QUICC_MEMORYUSAGE}/${QUICC_PLATFORM} 

done < "${SORTED_BENCHMARKS}"

cd ${current_dir}
rm ${SORTED_BENCHMARKS}
echo "Done"
