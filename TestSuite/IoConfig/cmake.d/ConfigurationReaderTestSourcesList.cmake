# Create list of sources for test case
set(MHDTestSources
   IoConfig/ConfigurationReader.cpp
   IoConfig/ConfigParts/IConfigurationPart.cpp
   IoConfig/ConfigParts/BoundaryPart.cpp
   IoConfig/ConfigParts/IoPart.cpp
   IoConfig/ConfigParts/ParallelPart.cpp
   IoConfig/ConfigParts/PhysicalPart.cpp
   IoConfig/ConfigParts/RunPart.cpp
   IoConfig/ConfigParts/TimesteppingPart.cpp
   IoConfig/ConfigParts/TruncationPart.cpp
   IoXml/IXmlReader.cpp
   IoXml/XmlFile.cpp
   IoTools/Formatter.cpp
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
