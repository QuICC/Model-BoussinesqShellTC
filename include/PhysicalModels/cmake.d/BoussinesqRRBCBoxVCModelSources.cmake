set(MHDModelSources
   ${GEOMHDISCC_SRC_DIR}/Equations/Box/Boussinesq/BoussinesqRRBCBoxVCTransport.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Box/Boussinesq/BoussinesqRRBCBoxVCMomentum.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Box/Boussinesq/BoussinesqRRBCBoxVCContinuity.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Tools/NoEigenTools.cpp
   ${GEOMHDISCC_SRC_DIR}/IoVariable/NusseltTags.cpp
   ${GEOMHDISCC_SRC_DIR}/IoVariable/Cartesian3DNusseltZWriter.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/RandomScalarState.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/RandomVectorState.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/CartesianExactStateIds.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/CartesianExactScalarState.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/CartesianExactVectorState.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/Visualizers/ScalarFieldVisualizer.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/Visualizers/VectorFieldVisualizer.cpp
)
