#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%% Create all the targets with a general implementation %%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Get all possible test directories
file(GLOB TmpDirs "../include/*")
foreach(TestDir ${TmpDirs})
   # Make sure it only contains directories
   if(IS_DIRECTORY ${TestDir})
      get_filename_component(short ${TestDir} NAME_WE)
      list(APPEND PossibleTestDirs ${short})
   endif(IS_DIRECTORY ${TestDir})
endforeach(TestDir ${TmpDirs})
#message(STATUS ${PossibleTestDirs})

# Get available test directories
file(GLOB TmpDirs "*")
foreach(TestDir ${TmpDirs})
   # Make sure it only contains directories
   if(IS_DIRECTORY ${TestDir})
      get_filename_component(short ${TestDir} NAME_WE)
      list(APPEND TestDirs ${short})
   endif(IS_DIRECTORY ${TestDir})
endforeach(TestDir ${TmpDirs})
#message(STATUS ${TestDirs})

# Loop over all test directories
foreach(TestDir ${PossibleTestDirs})

   # Check if test directory is present in TestSuite
   list(FIND TestDirs ${TestDir} present)

   # If the test directory is present
   if(present EQUAL -1)

      # Output warning message
      MESSAGE(STATUS "")
      MESSAGE(STATUS " !!! No tests available for ${TestDir}")

   else(present EQUAL -1)

      MESSAGE(STATUS "")
      MESSAGE(STATUS " Tests available for ${TestDir}:")

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
         message(STATUS " --> Added ${TestName}")

      endforeach(TestFile ${TestFiles})

   endif(present EQUAL -1)

endforeach(TestDir ${TestDirs})

add_custom_target(all_tests)
add_dependencies(all_tests ${AllTests})
