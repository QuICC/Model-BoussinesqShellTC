set(MHDTestDefTag GEOMHDISCC_SPATIALSCHEME)
set(MHDTestDefList "TTT" "TFT" "TFF" "FFF" "CFT" "WFT" "SLF" "WLF")

# Create list of sources for test case
set(MHDTestSources
   Enums/DimensionTools.cpp
   IoTools/Formatter.cpp
   IoTools/Formatter.cpp
   IoTools/IdToHuman.cpp
   IoTools/HumanToId.cpp
   IoVariable/VisualizationFileWriter.cpp
   IoVariable/VisualizationFileTags.cpp
   IoVariable/IVariableHdf5NWriter.cpp
   IoVariable/VariableHdf5Tags.cpp
   IoHdf5/IHdf5NWriter.cpp
   IoHdf5/IHdf5Writer.cpp
   IoHdf5/Hdf5File.cpp
   Resolutions/Resolution.cpp
   Resolutions/SimulationResolution.cpp
   Resolutions/CoreResolution.cpp
   Resolutions/TransformResolution.cpp
   Resolutions/TransformSetup.cpp
   Resolutions/Tools/IndexCounter.cpp
   Resolutions/Tools/RegularIndexCounter.cpp
   LoadSplitter/LoadSplitter.cpp
   LoadSplitter/Algorithms/SplittingDescription.cpp
   LoadSplitter/Algorithms/SerialSplitting.cpp
   LoadSplitter/Algorithms/SplittingAlgorithm.cpp
   LoadSplitter/Algorithms/SplittingTools.cpp
   FastTransforms/FftwTools.cpp
   FastTransforms/FftSetup.cpp
   SpatialSchemes/3D/TTTScheme.cpp
   SpatialSchemes/3D/TFTScheme.cpp
   SpatialSchemes/3D/TFFScheme.cpp
   SpatialSchemes/3D/FFFScheme.cpp
   SpatialSchemes/3D/IRegular3DScheme.cpp
   SpatialSchemes/ISpatialScheme.cpp
   SpatialSchemes/ISchemeCosts.cpp
   Variables/VariableBase.cpp
   Variables/VariableRequirement.cpp
   Variables/FieldRequirement.cpp
   Variables/RequirementTools.cpp
   SpectralOperators/IOperator.cpp
   SpectralOperators/ChebyshevOperator.cpp
   Equations/CouplingInformation.cpp
   Equations/EquationData.cpp
   Equations/EquationParameters.cpp
   Equations/IEquation.cpp
   Equations/IScalarEquation.cpp
   Equations/IVectorEquation.cpp
   Simulation/SimulationBoundary.cpp
   Debug/SerialDebugger.cpp
   Timers/ITimer.cpp
   Timers/SerialTimer.cpp
)

function(geomhdiscc_append_test condition srcs)
   if(${condition} STREQUAL "TTT")
       set(${srcs}
          Equations/Tools/EquationNoEigenTools.cpp
          PARENT_SCOPE)
   endif(${condition} STREQUAL "TTT")
   if(${condition} STREQUAL "TFT")
      set(${srcs}
         Equations/Tools/EquationEigen1DTools.cpp
         PARENT_SCOPE)
   endif(${condition} STREQUAL "TFT")
   if(${condition} STREQUAL "TFF")
      set(${srcs}
         Equations/Tools/EquationEigen2DTools.cpp
         PARENT_SCOPE)
   endif(${condition} STREQUAL "TFF")
   if(${condition} STREQUAL "FFF")
      set(${srcs}
         Equations/Tools/EquationEigen3DTools.cpp
         PARENT_SCOPE)
   endif(${condition} STREQUAL "FFF")
endfunction(geomhdiscc_append_test condition srcs)

# Include all files for the framework
include(../src/Framework/SourcesList.cmake)
foreach(src ${MHDSources})
   list(APPEND MHDTestSources Framework/${src})
endforeach(src ${MHDSources})
set(MHDSources )

# Include all files for the exceptions
include(../src/Exceptions/SourcesList.cmake)
foreach(src ${MHDSources})
   list(APPEND MHDTestSources Exceptions/${src})
endforeach(src ${MHDSources})
set(MHDSources )

# Include all files for the storage profiler
include(../src/StorageProfiler/SourcesList.cmake)
foreach(src ${MHDSources})
   list(APPEND MHDTestSources StorageProfiler/${src})
endforeach(src ${MHDSources})
set(MHDSources )

# Include all files for the profiler
include(../src/Profiler/SourcesList.cmake)
foreach(src ${MHDSources})
   list(APPEND MHDTestSources Profiler/${src})
endforeach(src ${MHDSources})
set(MHDSources )
