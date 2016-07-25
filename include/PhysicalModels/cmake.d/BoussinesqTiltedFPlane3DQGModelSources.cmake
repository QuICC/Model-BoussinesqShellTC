set(MHDModelSources
   ${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqTiltedFPlane3DQGMeanHeat.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqTiltedFPlane3DQGStreamfunction.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqTiltedFPlane3DQGVelocityZ.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqTiltedFPlane3DQGTransport.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqTiltedFPlane3DQGNoStreamfunction.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqTiltedFPlane3DQGNoVelocityZ.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqTiltedFPlane3DQGNoVorticityZ.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Tools/Eigen2DTools.cpp
   ${GEOMHDISCC_SRC_DIR}/IoVariable/EnergyTags.cpp
   ${GEOMHDISCC_SRC_DIR}/IoVariable/Cartesian1DScalarEnergyWriter.cpp
   ${GEOMHDISCC_SRC_DIR}/IoVariable/Cartesian1DStreamEnergyWriter.cpp
   ${GEOMHDISCC_SRC_DIR}/IoVariable/NusseltTags.cpp
   ${GEOMHDISCC_SRC_DIR}/IoVariable/Cartesian1DNusseltZWriter.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/CartesianExactStateIds.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/CartesianExactScalarState.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/RandomScalarState.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/Visualizers/ScalarFieldVisualizer.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/Visualizers/TiltedScalarFieldVisualizer.cpp
)
