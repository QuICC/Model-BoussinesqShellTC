set(MHDModelSources
   ${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBeta3DQGPerMeanHeat.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBeta3DQGPerStreamfunction.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBeta3DQGPerVelocityZ.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBeta3DQGPerTransport.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBeta3DQGPerVorticityZ.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Tools/Eigen2DTools.cpp
   ${GEOMHDISCC_SRC_DIR}/IoVariable/NusseltTags.cpp
   ${GEOMHDISCC_SRC_DIR}/IoVariable/Cartesian1DNusseltXWriter.cpp
   ${GEOMHDISCC_SRC_DIR}/IoVariable/EnergyTags.cpp
   ${GEOMHDISCC_SRC_DIR}/IoVariable/Cartesian1DScalarEnergyWriter.cpp
   ${GEOMHDISCC_SRC_DIR}/IoVariable/Cartesian1DStreamEnergyWriter.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/CartesianExactStateIds.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/CartesianExactScalarState.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/RandomScalarState.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/Visualizers/ScalarFieldVisualizer.cpp
)
