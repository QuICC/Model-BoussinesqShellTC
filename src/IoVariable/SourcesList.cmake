# Create list of sources
set(MHDSources
   StateFileTags.cpp
   StateFileWriter.cpp
   StateFileReader.cpp
   IVariableAsciiEWriter.cpp
   VariableHdf5Tags.cpp
   IVariableHdf5NWriter.cpp
   IVariableHdf5Reader.cpp
   VisualizationFileTags.cpp
   VisualizationFileWriter.cpp
   EnergyTags.cpp
   NusseltTags.cpp
   ContinuityTags.cpp
)

if(QUICC_SPATIALSCHEME STREQUAL "TTT")
   list(APPEND MHDSources
      Cartesian3DNusseltZWriter.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "TFT")
   list(APPEND MHDSources
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "TFF")
   list(APPEND MHDSources
      Cartesian1DMagneticEnergyWriter.cpp
      Cartesian1DNusseltDZWriter.cpp
      Cartesian1DNusseltXWriter.cpp
      Cartesian1DNusseltZWriter.cpp
      Cartesian1DPrimitiveEnergyWriter.cpp
      Cartesian1DScalarEnergyWriter.cpp
      Cartesian1DStreamEnergyWriter.cpp
      Cartesian1DTorPolEnergyWriter.cpp
      ContinuityWriter.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "FFF")
   list(APPEND MHDSources
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "AFT")
   list(APPEND MHDSources
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "CFT" OR QUICC_SPATIALSCHEME STREQUAL "WFT")
   list(APPEND MHDSources
      CylinderScalarEnergyWriter.cpp
      CylinderTorPolEnergyWriter.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "SLFL" OR QUICC_SPATIALSCHEME STREQUAL "SLFM")
   list(APPEND MHDSources
      ShellScalarEnergyWriter.cpp
      ShellTorPolEnergyWriter.cpp
      ShellTorPolEnergySpectraWriter.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "BLFL" OR QUICC_SPATIALSCHEME STREQUAL "BLFM" OR QUICC_SPATIALSCHEME STREQUAL "WLFL" OR QUICC_SPATIALSCHEME STREQUAL "WLFM")
   list(APPEND MHDSources
      SphereScalarEnergyWriter.cpp
      SphereTorPolEnergyWriter.cpp
      )
endif(QUICC_SPATIALSCHEME STREQUAL "TTT")
