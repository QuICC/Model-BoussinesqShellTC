include(ConfigureBenchmarks)
quicc_add_benchmark(Explicit
  MODEL "BoussinesqShellTC"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  )
