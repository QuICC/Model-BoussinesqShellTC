set(MHDModelSources
   ${GEOMHDISCC_SRC_DIR}/Equations/Sphere/Boussinesq/BoussinesqTCSphereTransport.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Sphere/Boussinesq/BoussinesqTCSphereMomentum.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Tools/EquationEigenSHmTools.cpp
   ${GEOMHDISCC_SRC_DIR}/IoVariable/EnergyTags.cpp
   ${GEOMHDISCC_SRC_DIR}/IoVariable/SphereScalarEnergyWriter.cpp
   ${GEOMHDISCC_SRC_DIR}/IoVariable/SphereTorPolEnergyWriter.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/RandomScalarState.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/RandomVectorState.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/SphereExactStateIds.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/SphereExactScalarState.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/SphereExactVectorState.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/Visualizers/ScalarFieldVisualizer.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/Visualizers/VectorFieldVisualizer.cpp
)
