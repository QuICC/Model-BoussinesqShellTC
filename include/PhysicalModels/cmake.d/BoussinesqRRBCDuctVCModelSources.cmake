set(MHDModelSources
   ${QUICC_SRC_DIR}/Equations/Box/Boussinesq/BoussinesqRRBCDuctVCTransport.cpp
   ${QUICC_SRC_DIR}/Equations/Box/Boussinesq/BoussinesqRRBCDuctVCMomentum.cpp
   ${QUICC_SRC_DIR}/Equations/Box/Boussinesq/BoussinesqRRBCDuctVCContinuity.cpp
   ${QUICC_SRC_DIR}/Equations/Tools/Eigen1DTools.cpp
   ${QUICC_SRC_DIR}/Generator/States/RandomScalarState.cpp
   ${QUICC_SRC_DIR}/Generator/States/RandomVectorState.cpp
   ${QUICC_SRC_DIR}/Generator/States/CartesianExactStateIds.cpp
   ${QUICC_SRC_DIR}/Generator/States/CartesianExactScalarState.cpp
   ${QUICC_SRC_DIR}/Generator/States/CartesianExactVectorState.cpp
   ${QUICC_SRC_DIR}/Generator/Visualizers/ScalarFieldVisualizer.cpp
   ${QUICC_SRC_DIR}/Generator/Visualizers/VectorFieldVisualizer.cpp
   ${QUICC_SRC_DIR}/IoVariable/ContinuityTags.cpp
   ${QUICC_SRC_DIR}/IoVariable/ContinuityWriter.cpp
   ${QUICC_SRC_DIR}/IoVariable/EnergyTags.cpp
)
