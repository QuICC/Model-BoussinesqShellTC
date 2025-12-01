target_sources(${QUICC_CURRENT_MODEL_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  ITCState.cpp
  ITCVisualization.cpp
  ITCModel.cpp
  ITCBackend.cpp
  Momentum.cpp
  MomentumKernel.cpp
  Transport.cpp
  TransportKernel.cpp
  )

add_subdirectory(Explicit)
