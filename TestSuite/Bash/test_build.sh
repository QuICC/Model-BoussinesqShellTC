#!/bin/bash

set -euf
set -o pipefail
#set -x

# Store current directory
current_dir=`pwd`

echo "Make directory structure ..."
mkdir Release-${3}-${4}
cd Release-${3}-${4}

echo "CMake configuration ..."
PROJECT_DIR=$1
QUICC_PLATFORM=$2
QUICC_MPIALGO=$3
QUICC_BOUNDARYMETHOD="Tau"
QUICC_FFT="FFTW"
QUICC_PYTHON="python34"
QUICC_SMARTPTR="TR1"
QUICC_SPLINALG="UmfPack"
QUICC_TIMESTEPPER="ImExRK3"
if [ "${QUICC_MPIALGO}" = "Tubular" ]
then
	QUICC_TRANSGROUPER="Transform"
else
	QUICC_TRANSGROUPER="Equation"
fi
cmake -DCMAKE_BUILD_TYPE=Release -DQUICC_PLATFORM=${QUICC_PLATFORM} ${PROJECT_DIR} >& /dev/null || true
cmake -DQUICC_BOUNDARYMETHOD=${QUICC_BOUNDARYMETHOD} -DQUICC_FFTPLAN=${QUICC_FFT} -DQUICC_FFTPLAN=Medium -DQUICC_MEMORYUSAGE=High -DQUICC_MPIALGO=${QUICC_MPIALGO} -DQUICC_PYTHON=${QUICC_PYTHON} -DQUICC_SMARTPTR=${QUICC_SMARTPTR} -DQUICC_SPLINALG=${QUICC_SPLINALG} -DQUICC_TIMESTEPPER=${QUICC_TIMESTEPPER} ${PROJECT_DIR} >& /dev/null || true
cmake -DQUICC_TRANSGROUPER=${QUICC_TRANSGROUPER} ${PROJECT_DIR}

echo "Installing Python files ..."
make install > /dev/null

echo "Building $4 ..."
CI_QUICC_PARAMAKE=$5
#make -j ${CI_QUICC_PARAMAKE} ${4}Config
#make -j ${CI_QUICC_PARAMAKE} ${4}State
make -j ${CI_QUICC_PARAMAKE} ${4}Model > /dev/null
#make -j ${CI_QUICC_PARAMAKE} ${4}Visu

cd ${current_dir}
echo "Done"
