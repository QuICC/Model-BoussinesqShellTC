set(MHDModelSources
   ${GEOMHDISCC_SRC_DIR}/Equations/Shell/Boussinesq/BoussinesqTCShellTransport.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Shell/Boussinesq/BoussinesqTCShellMomentum.cpp
   ${GEOMHDISCC_SRC_DIR}/Equations/Tools/EquationEigenSHlTools.cpp
   ${GEOMHDISCC_SRC_DIR}/IoVariable/EnergyTags.cpp
   ${GEOMHDISCC_SRC_DIR}/IoVariable/ShellScalarEnergyWriter.cpp
   ${GEOMHDISCC_SRC_DIR}/IoVariable/ShellTorPolEnergyWriter.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/RandomScalarState.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/RandomVectorState.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/ShellExactStateIds.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/ShellExactScalarState.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/States/ShellExactVectorState.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/Visualizers/ScalarFieldVisualizer.cpp
   ${GEOMHDISCC_SRC_DIR}/Generator/Visualizers/VectorFieldVisualizer.cpp
)
