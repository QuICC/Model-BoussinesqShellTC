#!/bin/bash

set -euf
set -o pipefail
#set -x

echo "Installing Python files ..."
make install > /dev/null

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

echo "Building ..."
short_cartesian_models="BoussinesqFPlane3DQG BoussinesqRRBCPlane "
short_cylindrical_models="BoussinesqRBCCylinder"
short_spherical_models="BoussinesqDynamoShellStd BoussinesqDynamoSphere"
long_cartesian_models="BoussinesqFPlane3DQG BoussinesqTiltedFPlane3DQG BoussinesqNoTiltedFPlane3DQG BoussinesqFPlaneNHBGE BoussinesqRRBCPlane"
long_cylindrical_models="BoussinesqRBCCylinder BoussinesqRRBCCylinder"
long_spherical_models="BoussinesqTCShellStd BoussinesqRTCShellStd BoussinesqDynamoShellStd BoussinesqTCShell BoussinesqRTCShell BoussinesqDynamoShell BoussinesqTCSphere BoussinesqRTCSphere BoussinesqDynamoSphere BoussinesqTCSphereStd BoussinesqRTCSphereStd BoussinesqDynamoSphereStd"
builds_list="${4}_models"
builds=${!builds_list}

if [ $# -le 4 ]
then
	CI_QUICC_PARAMAKE=1
else
	CI_QUICC_PARAMAKE=$5
fi

for b in ${builds}
do
	echo "   $b"
	#make -j ${CI_QUICC_PARAMAKE} ${BUILD_NAME}Config
	#make -j ${CI_QUICC_PARAMAKE} ${BUILD_NAME}State
	make -j ${CI_QUICC_PARAMAKE} ${b}Model > /dev/null
	#make -j ${CI_QUICC_PARAMAKE} ${BUILD_NAME}Visu
done
echo "Done"
