set(MHDModelSources
   ${QUICC_SRC_DIR}/Equations/Box/Boussinesq/BoussinesqDynamoPlaneTransport.cpp
   ${QUICC_SRC_DIR}/Equations/Box/Boussinesq/BoussinesqDynamoPlaneMomentum.cpp
   ${QUICC_SRC_DIR}/Equations/Box/Boussinesq/BoussinesqDynamoPlaneInduction.cpp
   ${QUICC_SRC_DIR}/Equations/Tools/Eigen2DTools.cpp
   ${QUICC_SRC_DIR}/Generator/States/RandomScalarState.cpp
   ${QUICC_SRC_DIR}/Generator/States/RandomVectorState.cpp
   ${QUICC_SRC_DIR}/Generator/States/CartesianExactStateIds.cpp
   ${QUICC_SRC_DIR}/Generator/States/CartesianExactScalarState.cpp
   ${QUICC_SRC_DIR}/Generator/States/CartesianExactVectorState.cpp
   ${QUICC_SRC_DIR}/Generator/Visualizers/ScalarFieldVisualizer.cpp
   ${QUICC_SRC_DIR}/Generator/Visualizers/ScalarFieldTrivialVisualizer.cpp
   ${QUICC_SRC_DIR}/Generator/Visualizers/VectorFieldVisualizer.cpp
   ${QUICC_SRC_DIR}/IoVariable/EnergyTags.cpp
   ${QUICC_SRC_DIR}/IoVariable/Cartesian1DScalarEnergyWriter.cpp
   ${QUICC_SRC_DIR}/IoVariable/Cartesian1DTorPolEnergyWriter.cpp
   ${QUICC_SRC_DIR}/IoVariable/NusseltTags.cpp
   ${QUICC_SRC_DIR}/IoVariable/Cartesian1DNusseltDZWriter.cpp
)
