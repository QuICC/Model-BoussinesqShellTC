# Create list of sources for test case
set(MHDTestSources
   Exceptions/Exception.cpp
   Profiler/ProfilerBase.cpp
   Framework/FrameworkBase.cpp
   Enums/DimensionTools.cpp
   IoTools/Formatter.cpp
   IoAscii/AsciiFile.cpp
   IoAscii/IAsciiWriter.cpp
   IoAscii/DirectAsciiWriter.cpp
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
   FastTransforms/FftSetup.cpp
   FastTransforms/FftwTools.cpp
   IoTools/VisualizeResolution.cpp
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
