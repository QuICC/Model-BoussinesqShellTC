#
# Create a configuration selection entry
#
function (geomhdiscc_create_choice choiceList choiceName chosen)
   foreach(choice ${${choiceList}})
      if(DEFINED str)
         set(str "${str}, ${choice}")
      else(DEFINED str)
         set(str "${choice}")
      endif(DEFINED str)
   endforeach(choice  ${${choiceList}})
   set(str "Available ${choiceName}(s): ${str}")
   # Get Length of choice list
   list(LENGTH ${choiceList} len)
   if(${len} EQUAL 1)
      set(${chosen} ${${choiceList}} CACHE STRING ${str} FORCE)
   else(${len} EQUAL 1)
      set(${chosen} ${${chosen}} CACHE STRING ${str} FORCE)
   endif(${len} EQUAL 1)
endfunction(geomhdiscc_create_choice)

# 
# Verify that selected setup is valid
#
function (geomhdiscc_check_choice choiceList choiceName chosen out)
   set(${out} OFF PARENT_SCOPE)
   if(NOT ${chosen} STREQUAL "")
      list(FIND ${choiceList} ${${chosen}} test)
      if(test EQUAL -1)
         message(SEND_ERROR "------------>>> Unkown ${${choiceName}} <<<------------")
      else(test EQUAL -1)
         set(${out} ON PARENT_SCOPE)
         message(STATUS " --> ${choiceName}: ${${chosen}}")
      endif(test EQUAL -1)
   else(NOT ${chosen} STREQUAL "")
      message(SEND_ERROR "------------>>> ${choiceName} is required <<<------------")
   endif(NOT ${chosen} STREQUAL "")
endfunction(geomhdiscc_check_choice)

#
# Provide a configuration selection entry
#
function (geomhdiscc_provide_choice choiceList choiceName chosen out)
   # Create selection entry
   geomhdiscc_create_choice(${choice} ${choiceList} ${choiceName} ${chosen})
   # Validate selection
   geomhdiscc_check_choice(${choice} ${choiceList} ${choiceName} ${chosen} ${out})
   set(${out} ${${out}} PARENT_SCOPE)
endfunction(geomhdiscc_provide_choice)

#
# Load platform file and setup requested platform
#
function (geomhdiscc_load_platform platform)
   include("cmake.d/platforms/${${platform}}.cmake")
   string(TOUPPER "${${platform}}" plat)
   add_definitions("-DGEOMHDISCC_ON_${plat}")
endfunction (geomhdiscc_load_platform)

#
# Add definition to flags
#
function (geomhdiscc_add_definition base)
   string(TOUPPER "${base}_${${base}}" def)
   add_definitions("-D${def}")
endfunction (geomhdiscc_add_definition base)

#
# Link to external libraries
#
function (geomhdiscc_link_external varName)
   set(tmp ${${varName}})
   if(NOT ${tmp} STREQUAL "")
      string(TOUPPER ${tmp} libName)
   endif(NOT ${tmp} STREQUAL "")
   # Check if variable is defined
   if(DEFINED GEOMHDISCC_LIBRARIES_${libName})
      list(LENGTH GEOMHDISCC_LIBRARIES_${libName} len)
      # Check the number of elements in list and rule out if not equal to 1
      if(${len} EQUAL 1)
         # Check if library requires lookup through find_package
         if(${GEOMHDISCC_LIBRARIES_${libName}} STREQUAL "auto")
            find_package(${libName})
            if(${${libName}_FOUND})
               set(GEOMHDISCC_LIBRARIES_${libName} ${${libName}_LIBRARIES})
            endif(${${libName}_FOUND})
         endif(${GEOMHDISCC_LIBRARIES_${libName}} STREQUAL "auto") 
      endif(${len} EQUAL 1)
   endif(DEFINED GEOMHDISCC_LIBRARIES_${libName})
   # Loop over all library directories to add 
   foreach(inc ${GEOMHDISCC_LIBDIR_${libName}})
      link_directories(${inc})
   endforeach(inc)
   # Loop over all libraries to link to
   foreach(lib ${GEOMHDISCC_LIBRARIES_${libName}})
      link_libraries(${lib})
   endforeach(lib)
   # Loop over all include directories to add 
   foreach(inc ${GEOMHDISCC_INCLUDES_${libName}})
      include_directories(${inc})
   endforeach(inc)
   # Get upper case compiler
   set(tmp ${GEOMHDISCC_COMPILER})
   if(NOT ${tmp} STREQUAL "")
      string(TOUPPER ${tmp} compUP)
   endif(NOT ${tmp} STREQUAL "")
   # Loop over all library directories to add that are compiler specific
   foreach(inc ${GEOMHDISCC_LIBDIR_${libName}_${compUP}})
      link_directories(${inc})
   endforeach(inc)
   # Loop over all include directories to add that are compiler specific
   foreach(inc ${GEOMHDISCC_INCLUDES_${libName}_${compUP}})
      include_directories(${inc})
   endforeach(inc)
endfunction (geomhdiscc_link_external libName)

#
# Recursive crawler through sources
#
function (geomhdiscc_append_sources MHDAll MHDPath MHDDirs)
   # Loop over all the subdirectories
   foreach(MHDDir ${${MHDDirs}})
      # Set new path including subdirectory
      set(MHDNext ${MHDPath}/${MHDDir})
      # Unset list of subdirectories and sources
      set(MHDSources)
      set(MHDSrcSubDirs)
      # Include SourcesList.cmake file
      include(${MHDNext}/SourcesList.cmake)
      # Check if there are additional subdirectories
      if(DEFINED MHDSrcSubDirs)
         geomhdiscc_append_sources(${MHDAll} ${MHDNext} MHDSrcSubDirs)
      endif()
      # Loop over all sources and add full path source to list
      foreach(MHDSource ${MHDSources})
         list(APPEND ${MHDAll} ${MHDNext}/${MHDSource})
      endforeach(MHDSource)
      # Update the list of all sources in parent scope
      set(${MHDAll} ${${MHDAll}} PARENT_SCOPE)
   endforeach(MHDDir)
endfunction ()

#
# Append path to list of files
#
function (geomhdiscc_add_path MHDList MHDPath)
   set(MHDTmp )
   foreach(MHDFile ${${MHDList}})
      list(APPEND MHDTmp ${MHDPath}/${MHDFile})
   endforeach(MHDFile ${${MHDList}})
   set(${MHDList} ${MHDTmp} PARENT_SCOPE)
endfunction()

#
# Create executable
#
function (geomhdiscc_add_executable MHDModel MHDScheme MHDPostfix MHDExecSrc MHDModelSrcs MHDAllSrcs)
   # Create new name for executable
   STRING(REGEX REPLACE "Model" ${MHDPostfix} ExecName ${MHDModel})

   # Create list of source files
   set(SrcsList ${MHDExecSrc} ${GEOMHDISCC_SRC_DIR}/PhysicalModels/${MHDModel}.cpp ${${MHDAllSrcs}} ${${MHDModelSrcs}})

   # Add executable to target list
   add_executable(${ExecName} ${SrcsList})

   # Set special properties of target
   set_target_properties(${ExecName} PROPERTIES OUTPUT_NAME
      ${ExecName} COMPILE_FLAGS "-DGEOMHDISCC_SPATIALSCHEME_${MHDScheme} -DGEOMHDISCC_RUNSIM_MODEL=${MHDModel}")

   # Show message
   message(STATUS " --> Added ${ExecName} executable")
endfunction ()
