set(MHDModelSources
   ${QUICC_SRC_DIR}/Equations/Shell/Boussinesq/BoussinesqCouetteShellMomentum.cpp
   ${QUICC_SRC_DIR}/Equations/Shell/Boussinesq/BoussinesqDynamoCouetteShellMomentum.cpp
   ${QUICC_SRC_DIR}/Equations/Shell/Boussinesq/BoussinesqDynamoCouetteShellInduction.cpp
   ${QUICC_SRC_DIR}/Equations/Tools/EigenSHlTools.cpp
   ${QUICC_SRC_DIR}/Equations/Tools/EigenSHlmTools.cpp
   ${QUICC_SRC_DIR}/IoVariable/EnergyTags.cpp
   ${QUICC_SRC_DIR}/IoVariable/ShellTorPolEnergyWriter.cpp
   ${QUICC_SRC_DIR}/Generator/States/RandomVectorState.cpp
   ${QUICC_SRC_DIR}/Generator/States/ShellExactStateIds.cpp
   ${QUICC_SRC_DIR}/Generator/States/ShellExactVectorState.cpp
   ${QUICC_SRC_DIR}/Generator/Visualizers/VectorFieldVisualizer.cpp
   ${QUICC_SRC_DIR}/Generator/Visualizers/SphericalVerticalFieldVisualizer.cpp
)
