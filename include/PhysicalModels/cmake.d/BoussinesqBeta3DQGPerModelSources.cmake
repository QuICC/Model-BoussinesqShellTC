set(MHDModelSources
   ${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBeta3DQGPerMeanHeat.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBeta3DQGPerStreamfunction.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBeta3DQGPerVelocityZ.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBeta3DQGPerTransport.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBeta3DQGPerVorticityZ.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBeta3DQGPerKineticNRG.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBeta3DQGPerZonalKineticNRG.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBeta3DQGPerNonZonalKineticNRG.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBeta3DQGPerMeanVelocityY.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBeta3DQGPerMeanVelocityZ.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Tools/EquationEigen2DTools.cpp
   ${GEOMHDISCC_SRC_DIR}/IoVariable/NusseltTags.cpp
   ${GEOMHDISCC_SRC_DIR}/IoVariable/NusseltWriter.cpp
   ${GEOMHDISCC_SRC_DIR}/IoVariable/EnergyTags.cpp
   ${GEOMHDISCC_SRC_DIR}/IoVariable/Cartesian1DScalarEnergyWriter.cpp
   ${GEOMHDISCC_SRC_DIR}/IoVariable/Cartesian1DStreamEnergyWriter.cpp
   ${GEOMHDISCC_SRC_DIR}/IoVariable/KineticEnergyTags.cpp
   ${GEOMHDISCC_SRC_DIR}/IoVariable/NusseltBeta3DQGPerWriter.cpp
   ${GEOMHDISCC_SRC_DIR}/IoVariable/KineticEnergyBeta3DQGPerWriter.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/CartesianExactStateIds.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/CartesianExactScalarState.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/RandomScalarState.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/Visualizers/ScalarFieldVisualizer.cpp
)
