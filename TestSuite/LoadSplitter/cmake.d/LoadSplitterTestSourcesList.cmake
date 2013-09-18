# Required additional definitions
set(MHDTestDefTag GEOMHDISCC_SPATIALSCHEME)
set(MHDTestDefList "TTT" "TFT" "TFF" "FFF" "CFT" "WFT" "SLF" "WLF")

# Create list of sources for test case
set(MHDTestSources
   Exceptions/Exception.cpp
   Profiler/ProfilerBase.cpp
   Framework/FrameworkBase.cpp
   Enums/DimensionTools.cpp
   IoTools/Formatter.cpp
   IoXml/XmlFile.cpp
   IoXml/IXmlWriter.cpp
   IoXml/GxlWriter.cpp
   Resolutions/Resolution.cpp
   Resolutions/CoreResolution.cpp
   Resolutions/SimulationResolution.cpp
   Resolutions/TransformResolution.cpp
   LoadSplitter/LoadSplitter.cpp
   LoadSplitter/Algorithms/SplittingTools.cpp
   LoadSplitter/Algorithms/SplittingDescription.cpp
   LoadSplitter/Algorithms/SplittingAlgorithm.cpp
   SpatialSchemes/ISchemeCosts.cpp
   SpatialSchemes/ISpatialScheme.cpp
   SpatialSchemes/3D/IRegular3DScheme.cpp
   SpatialSchemes/3D/TFTScheme.cpp
   Resolutions/TransformSetup.cpp
   Resolutions/Tools/IndexCounter.cpp
   Resolutions/Tools/RegularIndexCounter.cpp
   FastTransforms/FftSetup.cpp
   FastTransforms/FftwTools.cpp
)

if(GEOMHDISCC_MPI)
   list(APPEND MHDTestSources
      Framework/MpiFramework.cpp
      Profiler/MpiProfiler.cpp
      LoadSplitter/Algorithms/SingleSplitting.cpp
      LoadSplitter/Algorithms/TubularSplitting.cpp
      )
else(GEOMHDISCC_MPI)
   list(APPEND MHDTestSources
      Framework/SerialFramework.cpp
      Profiler/SerialProfiler.cpp
      LoadSplitter/Algorithms/SerialSplitting.cpp
      )
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
