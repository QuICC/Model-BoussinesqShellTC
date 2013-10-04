# Required additional definitions
set(MHDTestDefTag GEOMHDISCC_SPATIALSCHEME)
set(MHDTestDefList "TTT" "TFT" "TFF" "FFF" "CFT" "WFT" "SLF" "WLF")

# Create list of sources for test case
set(MHDTestSources
   Base/MathConstants.cpp
   Enums/DimensionTools.cpp
   IoTools/Formatter.cpp
   Resolutions/Resolution.cpp
   Resolutions/CoreResolution.cpp
   Resolutions/SimulationResolution.cpp
   Resolutions/TransformResolution.cpp
   Resolutions/TransformSetup.cpp
   Resolutions/Tools/IndexCounter.cpp
   Resolutions/Tools/RegularIndexCounter.cpp
   LoadSplitter/LoadSplitter.cpp
   LoadSplitter/Algorithms/SplittingTools.cpp
   LoadSplitter/Algorithms/SplittingDescription.cpp
   LoadSplitter/Algorithms/SplittingAlgorithm.cpp
   LoadSplitter/Algorithms/SerialSplitting.cpp
   SpatialSchemes/ISchemeCosts.cpp
   SpatialSchemes/ISpatialScheme.cpp
   SpatialSchemes/3D/IRegular3DScheme.cpp
   SpatialSchemes/3D/TTTScheme.cpp
   SpatialSchemes/3D/TFTScheme.cpp
   SpatialSchemes/3D/TFFScheme.cpp
   SpatialSchemes/3D/FFFScheme.cpp
   FastTransforms/FftSetup.cpp
   FastTransforms/FftwTools.cpp
   FastTransforms/FftwLibrary.cpp
   FastTransforms/FftwTransform.cpp
   FastTransforms/ChebyshevFftwTransform.cpp
   Variables/VariableRequirement.cpp
   Variables/FieldRequirement.cpp
   Variables/RequirementTools.cpp
   Variables/VariableBase.cpp
   TypeSelectors/ParallelSelector.cpp
   TransformGroupers/IForwardGrouper.cpp
   TransformGroupers/IBackwardGrouper.cpp
   TransformConfigurators/ForwardConfigurator.cpp
   TransformConfigurators/ForwardSerialConfigurator.cpp
   TransformConfigurators/BackwardConfigurator.cpp
   TransformConfigurators/BackwardSerialConfigurator.cpp
   TransformCoordinators/TransformCoordinatorTools.cpp
   SpectralOperators/IOperator.cpp
   SpectralOperators/ChebyshevOperator.cpp
   Equations/CouplingInformation.cpp
   Equations/EquationData.cpp
   Equations/EquationParameters.cpp
   Equations/IEquation.cpp
   Equations/IScalarEquation.cpp
   Equations/IVectorEquation.cpp
   Equations/Tests/TestSpatialSchemeForwardScalar.cpp
   Equations/Tests/TestSpatialSchemeBackwardScalar.cpp
   IoTools/IdToHuman.cpp
   IoTools/HumanToId.cpp
   Timers/ITimer.cpp
   Timers/SerialTimer.cpp
   Simulation/SimulationBoundary.cpp
   SpectralOperators/IOperator.cpp
   SpectralOperators/ChebyshevOperator.cpp
   Communicators/Converters/NoIndexConv.cpp
   Communicators/Converters/SHIndexConv.cpp
   Communicators/Converters/PMIndexConv.cpp
)

function(geomhdiscc_append_test condition srcs)
   if(${condition} STREQUAL "TFT")
      set(${srcs}
         Equations/Tools/Equation1DEigenTools.cpp
         PARENT_SCOPE)
   endif(${condition} STREQUAL "TFT")
   if(${condition} STREQUAL "TFF")
      set(${srcs}
         Equations/Tools/Equation2DEigenTools.cpp
         PARENT_SCOPE)
   endif(${condition} STREQUAL "TFF")
endfunction(geomhdiscc_append_test condition srcs)

if(GEOMHDISCC_MPI)
   list(APPEND MHDTestSources
      Timers/MpiTimer.cpp
      Communicators/CommunicationBuffer.cpp
      Communicators/Converters/MpiConverterTools.cpp
      LoadSplitter/Algorithms/SingleSplitting.cpp
      LoadSplitter/Algorithms/FixedSplitting.cpp
      LoadSplitter/Algorithms/TubularSplitting.cpp
   )
   if(GEOMHDISCC_MPIALGO STREQUAL "Single1D" OR GEOMHDISCC_MPIALGO STREQUAL "Fixed")
      list(APPEND MHDTestSources
         TransformConfigurators/BackwardSingle1DConfigurator.cpp
         TransformConfigurators/ForwardSingle1DConfigurator.cpp
      )
   endif(GEOMHDISCC_MPIALGO STREQUAL "Single1D" OR GEOMHDISCC_MPIALGO STREQUAL "Fixed")
   if(GEOMHDISCC_MPIALGO STREQUAL "Single2D")
      list(APPEND MHDTestSources
         TransformConfigurators/BackwardSingle2DConfigurator.cpp
         TransformConfigurators/ForwardSingle2DConfigurator.cpp
      )
   endif(GEOMHDISCC_MPIALGO STREQUAL "Single2D")
   if(GEOMHDISCC_MPIALGO STREQUAL "Tubular")
      list(APPEND MHDTestSources
         TransformConfigurators/BackwardTubularConfigurator.cpp
         TransformConfigurators/ForwardTubularConfigurator.cpp
      )
   endif(GEOMHDISCC_MPIALGO STREQUAL "Tubular")
endif(GEOMHDISCC_MPI)

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

# Include all files for the profiler
include(../src/Debug/SourcesList.cmake)
foreach(src ${MHDSources})
   list(APPEND MHDTestSources Debug/${src})
endforeach(src ${MHDSources})
set(MHDSources )
