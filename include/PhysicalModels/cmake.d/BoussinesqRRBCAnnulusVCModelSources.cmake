set(MHDModelSources
   ${GEOMHDISCC_SRC_DIR}/Equations/Annulus/Boussinesq/BoussinesqRRBCAnnulusVCTransport.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Annulus/Boussinesq/BoussinesqRRBCAnnulusVCMomentum.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Annulus/Boussinesq/BoussinesqRRBCAnnulusVCContinuity.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Tools/Eigen1DTools.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/RandomScalarState.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/RandomVectorState.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/AnnulusExactStateIds.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/AnnulusExactScalarState.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/AnnulusExactVectorState.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/Visualizers/ScalarFieldVisualizer.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/Visualizers/VectorFieldVisualizer.cpp
)
