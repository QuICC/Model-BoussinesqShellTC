set(MHDModelSources
   ${QUICC_SRC_DIR}/Equations/Sphere/Boussinesq/BoussinesqRTCSphereTransport.cpp
   ${QUICC_SRC_DIR}/Equations/Sphere/Boussinesq/BoussinesqRTCSphereMomentum.cpp
   ${QUICC_SRC_DIR}/Equations/Tools/EigenSHlTools.cpp
   ${QUICC_SRC_DIR}/Equations/Tools/EigenSHlmTools.cpp
   ${QUICC_SRC_DIR}/IoVariable/EnergyTags.cpp
   ${QUICC_SRC_DIR}/IoVariable/SphereScalarEnergyWriter.cpp
   ${QUICC_SRC_DIR}/IoVariable/SphereTorPolEnergyWriter.cpp
   ${QUICC_SRC_DIR}/Generator/States/RandomScalarState.cpp
   ${QUICC_SRC_DIR}/Generator/States/RandomVectorState.cpp
   ${QUICC_SRC_DIR}/Generator/States/SphereExactStateIds.cpp
   ${QUICC_SRC_DIR}/Generator/States/SphereExactScalarState.cpp
   ${QUICC_SRC_DIR}/Generator/States/SphereExactVectorState.cpp
   ${QUICC_SRC_DIR}/Generator/Visualizers/ScalarFieldVisualizer.cpp
   ${QUICC_SRC_DIR}/Generator/Visualizers/VectorFieldVisualizer.cpp
)
