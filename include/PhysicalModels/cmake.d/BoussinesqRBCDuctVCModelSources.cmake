set(MHDModelSources
   ${GEOMHDISCC_SRC_DIR}/Equations/Box/Boussinesq/BoussinesqRBCDuctVCTransport.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Box/Boussinesq/BoussinesqRBCDuctVCMomentum.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Box/Boussinesq/BoussinesqRBCDuctVCContinuity.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Tools/EquationEigen1DTools.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/RandomScalarState.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/RandomVectorState.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/CartesianExactStateIds.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/CartesianExactScalarState.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/CartesianExactVectorState.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/Visualizers/ScalarFieldVisualizer.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/Visualizers/VectorFieldVisualizer.cpp
   ${GEOMHDISCC_SRC_DIR}/IoVariable/ContinuityTags.cpp
   ${GEOMHDISCC_SRC_DIR}/IoVariable/ContinuityWriter.cpp
   ${GEOMHDISCC_SRC_DIR}/IoVariable/EnergyTags.cpp
)
