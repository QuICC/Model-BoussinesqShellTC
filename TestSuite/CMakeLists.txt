#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%% Create all the targets with a general implementation %%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Get List of directories
file(GLOB TestDirs "*")

# Loop over all test directories
foreach(TestDir ${TestDirs})

   # Make sure it is a directory
   if(IS_DIRECTORY ${TestDir})

      # Get all test executables in directory
      file(GLOB TestFiles ${TestDir}/*.cpp)

      # Create targets for each test executable
      foreach(TestFile ${TestFiles})

         # Extract filename from full path
         get_filename_component(TestName ${TestFile} NAME_WE)

         # Get test specific sources list
         include(${TestDir}/cmake.d/${TestName}SourcesList.cmake)

         # Initialise sources list
         set(TestSrcsList ${TestFile})

         # Loop over all sources and add full path source to list
         foreach(MHDSource ${MHDTestSources})
            list(APPEND TestSrcsList ${GEOMHDISCC_SRC_DIR}/${MHDSource})
         endforeach(MHDSource)

         # Add test executable to target list
         add_executable(${TestName} ${TestSrcsList})
         target_link_libraries(${TestName} "gtest")
         list(APPEND AllTests ${TestName})

         # Add message
         message(STATUS " --> Added ${TestName} test case")

      endforeach(TestFile ${TestFiles})

   endif(IS_DIRECTORY ${TestDir})

endforeach(TestDir ${TestDirs})

add_custom_target(all_tests)
add_dependencies(all_tests ${AllTests})
