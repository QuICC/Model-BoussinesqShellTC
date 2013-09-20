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
   Resolutions/TransformSetup.cpp
   Resolutions/Tools/IndexCounter.cpp
   Resolutions/Tools/RegularIndexCounter.cpp
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
   Equations/Tests/TestTFTForwardScalar.cpp
   Equations/Tests/TestTFTBackwardScalar.cpp
   Equations/Tests/TestTTTForwardScalar.cpp
   Equations/Tests/TestTTTBackwardScalar.cpp
   Equations/Tests/TestTFTForwardScalar.cpp
   Equations/Tests/TestTFTBackwardScalar.cpp
   Equations/Tests/TestTFFForwardScalar.cpp
   Equations/Tests/TestTFFBackwardScalar.cpp
   Equations/Tests/TestFFFForwardScalar.cpp
   Equations/Tests/TestFFFBackwardScalar.cpp
   Equations/Tests/TestCFTForwardScalar.cpp
   Equations/Tests/TestCFTBackwardScalar.cpp
   Equations/Tests/TestSLFForwardScalar.cpp
   Equations/Tests/TestSLFBackwardScalar.cpp
   Equations/Tests/TestWFTForwardScalar.cpp
   Equations/Tests/TestWFTBackwardScalar.cpp
   Equations/Tests/TestWLFForwardScalar.cpp
   Equations/Tests/TestWLFBackwardScalar.cpp
   IoTools/IdToHuman.cpp
   IoTools/HumanToId.cpp
   Timers/ITimer.cpp
   Timers/SerialTimer.cpp
)

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
