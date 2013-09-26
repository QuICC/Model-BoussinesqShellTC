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

         # Create test executable name
         set(TestExecName ${TestName})
         if(NOT ${MHDTestBackendFlag} STREQUAL "")
            set(TestExecName ${${MHDTestBackendFlag}}${TestExecName})
         endif(NOT ${MHDTestBackendFlag} STREQUAL "")
         if(NOT ${MHDTestPrefix} STREQUAL "")
            set(TestExecName ${MHDTestPrefix}${TestExecName})
         endif(NOT ${MHDTestPrefix} STREQUAL "")
         if(NOT ${MHDTestDefTag} STREQUAL "")
            foreach(TestDef ${MHDTestDefList})
               set(TestExecName2 ${TestDef}${TestExecName})
               geomhdiscc_append_test(${TestDef} TestDefSrcs)
               # Add test executable to target list
               set(TestDefSrcsList ${TestSrcsList})
               foreach(MHDSource ${TestDefSrcs})
                  list(APPEND TestDefSrcsList ${GEOMHDISCC_SRC_DIR}/${MHDSource})
               endforeach(MHDSource)
               add_executable(${TestExecName2} ${TestDefSrcsList})
               target_link_libraries(${TestExecName2} "gtest")
               list(APPEND AllTests ${TestExecName2})
               set_target_properties(${TestExecName2} PROPERTIES COMPILE_FLAGS
                  "-D${MHDTestDefTag}_${TestDef}")
               # Add message
               message(STATUS " --> Added ${TestExecName2}")
            endforeach(TestDef ${MHDTestDefList})
         else(NOT ${MHDTestDefTag} STREQUAL "")
            # Add test executable to target list
            add_executable(${TestExecName} ${TestSrcsList})
            target_link_libraries(${TestExecName} "gtest")
            list(APPEND AllTests ${TestExecName})

            # Add message
            message(STATUS " --> Added ${TestExecName}")
         endif(NOT ${MHDTestDefTag} STREQUAL "")

         set(MHDTestBackendFlag "")
         set(MHDTestPrefix "")
         set(MHDTestDefTag "")
         set(MHDTestDefList "")

      endforeach(TestFile ${TestFiles})

   endif(present EQUAL -1)

endforeach(TestDir ${TestDirs})

add_custom_target(all_tests)
add_dependencies(all_tests ${AllTests})
