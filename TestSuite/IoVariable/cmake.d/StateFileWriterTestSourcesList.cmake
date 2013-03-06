# Create list of sources for test case
set(MHDTestSources
   IoTools/Formatter.cpp
   IoVariable/StateFileWriter.cpp
   IoVariable/StateFileTags.cpp
   IoVariable/IVariableHdf5NWriter.cpp
   IoVariable/VariableHdf5Tags.cpp
   IoHdf5/IHdf5NWriter.cpp
   IoHdf5/IHdf5Writer.cpp
   IoHdf5/Hdf5File.cpp
   Resolutions/Resolution.cpp
   Resolutions/SimulationResolution.cpp
   Resolutions/CoreResolution.cpp
   Resolutions/TransformResolution.cpp
   LoadSplitter/LoadSplitter.cpp
   LoadSplitter/Algorithms/SplittingDescription.cpp
   LoadSplitter/Algorithms/SerialSplitting.cpp
   LoadSplitter/Algorithms/SplittingAlgorithm.cpp
   LoadSplitter/Algorithms/SplittingTools.cpp
   FastTransforms/FftwTools.cpp
   FastTransforms/FftSetup.cpp
   SpatialSchemes/3D/TFTScheme.cpp
   SpatialSchemes/3D/IRegular3DScheme.cpp
   SpatialSchemes/ISpatialScheme.cpp
   SpatialSchemes/ISchemeCosts.cpp
   Variables/VariableBase.cpp
   IoTools/IdToHuman.cpp
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
