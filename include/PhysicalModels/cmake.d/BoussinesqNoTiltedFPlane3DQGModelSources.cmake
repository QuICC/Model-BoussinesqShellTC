set(MHDModelSources
   ${QUICC_SRC_DIR}/Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqNoTiltedFPlane3DQGMeanHeat.cpp
   ${QUICC_SRC_DIR}/Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqNoTiltedFPlane3DQGTransport.cpp
   ${QUICC_SRC_DIR}/Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqNoTiltedFPlane3DQGNoStreamfunction.cpp
   ${QUICC_SRC_DIR}/Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqNoTiltedFPlane3DQGNoVelocityZ.cpp
   ${QUICC_SRC_DIR}/Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqNoTiltedFPlane3DQGNoVorticityZ.cpp
   ${QUICC_SRC_DIR}/Equations/Tools/Eigen2DTools.cpp
   ${QUICC_SRC_DIR}/IoVariable/EnergyTags.cpp
   ${QUICC_SRC_DIR}/IoVariable/Cartesian1DScalarEnergyWriter.cpp
   ${QUICC_SRC_DIR}/IoVariable/Cartesian1DStreamEnergyWriter.cpp
   ${QUICC_SRC_DIR}/IoVariable/NusseltTags.cpp
   ${QUICC_SRC_DIR}/IoVariable/Cartesian1DNusseltZWriter.cpp
   ${QUICC_SRC_DIR}/Generator/States/CartesianExactStateIds.cpp
   ${QUICC_SRC_DIR}/Generator/States/CartesianExactScalarState.cpp
   ${QUICC_SRC_DIR}/Generator/States/RandomScalarState.cpp
   ${QUICC_SRC_DIR}/Generator/Visualizers/ScalarFieldVisualizer.cpp
   ${QUICC_SRC_DIR}/Generator/Visualizers/TiltedScalarFieldVisualizer.cpp
)
