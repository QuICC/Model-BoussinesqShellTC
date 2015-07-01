set(MHDModelSources
   ${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqNoTiltedFPlane3DQGMeanHeat.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqNoTiltedFPlane3DQGTransport.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqNoTiltedFPlane3DQGNoStreamfunction.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqNoTiltedFPlane3DQGNoVelocityZ.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqNoTiltedFPlane3DQGNoVorticityZ.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Tools/EquationEigen2DTools.cpp
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
