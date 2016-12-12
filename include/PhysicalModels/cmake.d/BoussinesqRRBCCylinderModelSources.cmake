set(MHDModelSources
   ${GEOMHDISCC_SRC_DIR}/Equations/Cylinder/Boussinesq/BoussinesqRRBCCylinderTransport.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Cylinder/Boussinesq/BoussinesqRRBCCylinderMomentum.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Tools/Eigen1DTools.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/RandomScalarState.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/RandomVectorState.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/CylinderExactStateIds.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/CylinderExactScalarState.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/CylinderExactVectorState.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/Visualizers/ScalarFieldVisualizer.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/Visualizers/ScalarFieldTrivialVisualizer.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/Visualizers/VectorFieldVisualizer.cpp
)