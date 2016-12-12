set(MHDModelSources
	${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqDynamo3DQGMeanHeat.cpp
	${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqDynamo3DQGfbx.cpp
	${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqDynamo3DQGfby.cpp
	${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqDynamo3DQGfbz.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqDynamo3DQGBy.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqDynamo3DQGBx.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqDynamo3DQGEmfx.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqDynamo3DQGEmfy.cpp
	${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqDynamo3DQGStreamfunction.cpp
	${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqDynamo3DQGPressure.cpp
	${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqDynamo3DQGVelocityZ.cpp
	${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqDynamo3DQGTransport.cpp
	${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqDynamo3DQGVorticityZ.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Tools/Eigen2DTools.cpp
   ${GEOMHDISCC_SRC_DIR}/IoVariable/NusseltTags.cpp
   ${GEOMHDISCC_SRC_DIR}/IoStats/AvgTags.cpp
   ${GEOMHDISCC_SRC_DIR}/IoStats/Cartesian1DScalarAvgWriter.cpp
   ${GEOMHDISCC_SRC_DIR}/IoVariable/Cartesian1DNusseltZWriter.cpp
   ${GEOMHDISCC_SRC_DIR}/IoVariable/EnergyTags.cpp
   ${GEOMHDISCC_SRC_DIR}/IoVariable/Cartesian1DScalarEnergyWriter.cpp
   ${GEOMHDISCC_SRC_DIR}/IoVariable/Cartesian1DStreamEnergyWriter.cpp
   ${GEOMHDISCC_SRC_DIR}/IoVariable/Cartesian1DMagneticEnergyWriter.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/CartesianExactStateIds.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/CartesianExactScalarState.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/RandomScalarState.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/Visualizers/ScalarFieldVisualizer.cpp
)