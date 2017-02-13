#
# Create a configuration selection entry
#
function (quicc_create_choice choiceList choiceName chosen)
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
endfunction(quicc_create_choice)

# 
# Verify that selected setup is valid
#
function (quicc_check_choice choiceList choiceName chosen out)
   set(${out} OFF PARENT_SCOPE)
   if(NOT ${chosen} STREQUAL "")
      list(FIND ${choiceList} ${${chosen}} test)
      if(test EQUAL -1)
         message(SEND_ERROR "------------>>> Unkown ${${choiceName}} <<<------------")
      else(test EQUAL -1)
         if(NOT ${chosen} STREQUAL "None")
            set(${out} ON PARENT_SCOPE)
         endif(NOT ${chosen} STREQUAL "None")
         message(STATUS " --> ${choiceName}: ${${chosen}}")
      endif(test EQUAL -1)
   else(NOT ${chosen} STREQUAL "")
      message(SEND_ERROR "------------>>> ${choiceName} is required <<<------------")
   endif(NOT ${chosen} STREQUAL "")
endfunction(quicc_check_choice)

#
# Provide a configuration selection entry
#
function (quicc_provide_choice choiceList choiceName chosen out)
   # Create selection entry
   quicc_create_choice(${choice} ${choiceList} ${choiceName} ${chosen})
   # Validate selection
   quicc_check_choice(${choice} ${choiceList} ${choiceName} ${chosen} ${out})
   set(${out} ${${out}} PARENT_SCOPE)
endfunction(quicc_provide_choice)

#
# Load platform file and setup requested platform
#
function (quicc_load_platform platform)
   include("cmake.d/platforms/${${platform}}.cmake")
   string(TOUPPER "${${platform}}" plat)
   add_definitions("-DQUICC_ON_${plat}")
endfunction (quicc_load_platform)

#
# Add definition to flags
#
function (quicc_add_definition base)
   string(TOUPPER "${base}_${${base}}" def)
   add_definitions("-D${def}")
endfunction (quicc_add_definition base)

#
# Link to external libraries
#
function (quicc_link_external varName)
   set(tmp ${${varName}})
   if(${ARGC} EQUAL 2 AND NOT ${ARGV1} STREQUAL "None")
      string(TOUPPER ${${ARGV1}} model)
      set(tmp "${tmp}_${model}")
   endif(${ARGC} EQUAL 2 AND NOT ${ARGV1} STREQUAL "None")
   if(NOT ${tmp} STREQUAL "")
      string(TOUPPER ${tmp} libName)
   endif(NOT ${tmp} STREQUAL "")
   # Check if variable is defined
   if(DEFINED QUICC_LIBRARIES_${libName})
      list(LENGTH QUICC_LIBRARIES_${libName} len)
      # Check the number of elements in list and rule out if not equal to 1
      if(${len} EQUAL 1)
         # Check if library requires lookup through find_package
         if(${QUICC_LIBRARIES_${libName}} STREQUAL "auto")
            find_package(${libName})
            if(${${libName}_FOUND})
               set(QUICC_LIBRARIES_${libName} ${${libName}_LIBRARIES})
            endif(${${libName}_FOUND})
         endif(${QUICC_LIBRARIES_${libName}} STREQUAL "auto") 
      endif(${len} EQUAL 1)
   endif(DEFINED QUICC_LIBRARIES_${libName})
   # Loop over all library directories to add 
   foreach(inc ${QUICC_LIBDIR_${libName}})
      link_directories(${inc})
   endforeach(inc)
   # Loop over all libraries to link to
   foreach(lib ${QUICC_LIBRARIES_${libName}})
      link_libraries(${lib})
   endforeach(lib)
   # Loop over all include directories to add 
   foreach(inc ${QUICC_INCLUDES_${libName}})
      include_directories(${inc})
   endforeach(inc)
   # Get upper case compiler
   set(tmp ${QUICC_COMPILER})
   if(NOT ${tmp} STREQUAL "")
      string(TOUPPER ${tmp} compUP)
   endif(NOT ${tmp} STREQUAL "")
   # Loop over all library directories to add that are compiler specific
   foreach(inc ${QUICC_LIBDIR_${libName}_${compUP}})
      link_directories(${inc})
   endforeach(inc)
   # Loop over all include directories to add that are compiler specific
   foreach(inc ${QUICC_INCLUDES_${libName}_${compUP}})
      include_directories(${inc})
   endforeach(inc)
endfunction (quicc_link_external libName)

#
# Recursive crawler through sources
#
function (quicc_append_sources MHDAll MHDPath MHDDirs)
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
         quicc_append_sources(${MHDAll} ${MHDNext} MHDSrcSubDirs)
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
function (quicc_add_path MHDList MHDPath)
   set(MHDTmp )
   foreach(MHDFile ${${MHDList}})
      list(APPEND MHDTmp ${MHDPath}/${MHDFile})
   endforeach(MHDFile ${${MHDList}})
   set(${MHDList} ${MHDTmp} PARENT_SCOPE)
endfunction()

#
# Create executable
#
function (quicc_add_executable ModelPath StateRef Scheme SchemeDim MHDForm Postfix MHDExecSrc MHDModelSrcs MHDAllSrcs)
   # Create simple model name
   string(REGEX REPLACE ${StateRef}/ "" ModelId ${ModelPath})
   string(REGEX REPLACE "/" "" ModelName ${ModelId})
   string(REGEX REPLACE "/" "::" CPPModel ${ModelId})

   # Create new name for executable
   set(ExecName ${ModelName}${Postfix})

   # Create upper case scheme name
   string(TOUPPER ${Scheme} UpScheme)

   # Create list of source files
   set(SrcsList ${MHDExecSrc} ${${MHDAllSrcs}} ${${MHDModelSrcs}})

   # Add executable to target list
   add_executable(${ExecName} ${SrcsList})

   # Set special properties of target
   if(${MHDForm} STREQUAL "DEFAULT")
      set_target_properties(${ExecName} PROPERTIES OUTPUT_NAME
         ${ExecName} COMPILE_FLAGS "-DQUICC_SPATIALSCHEME_${UpScheme} -DQUICC_SPATIALDIMENSION_${SchemeDim} -DQUICC_MODEL_PATH=${StateRef} -DQUICC_RUNSIM_PATH=${ModelId} -DQUICC_RUNSIM_CPPMODEL=${CPPModel}")
   else(${MHDForm} STREQUAL "DEFAULT")
      set_target_properties(${ExecName} PROPERTIES OUTPUT_NAME
         ${ExecName} COMPILE_FLAGS "-DQUICC_SPATIALSCHEME_${UpScheme} -DQUICC_SPATIALDIMENSION_${SchemeDim} -DQUICC_SPATIALSCHEME_${UpScheme}_${MHDForm} -DQUICC_MODEL_PATH=${StateRef} -DQUICC_RUNSIM_PATH=${ModelId} -DQUICC_RUNSIM_CPPMODEL=${CPPModel}")
   endif(${MHDForm} STREQUAL "DEFAULT")

   # Show message
   message(STATUS "    --> added ${ExecName}")
endfunction ()

#
# Create executable
#
function (quicc_add_test MHDModel MHDScheme MHDSchemeDim MHDForm MHDPostfix MHDExecSrc MHDModelSrcs MHDAllSrcs)
   # Create new name for executable
   STRING(REGEX REPLACE "Test" ${MHDPostfix} ExecName ${MHDModel})

   # Create upper case scheme name
   string(TOUPPER ${MHDScheme} UpMHDScheme)

   # Create list of source files
   set(SrcsList ${MHDExecSrc} ${${MHDAllSrcs}} ${${MHDModelSrcs}})

   # Add executable to target list
   add_executable(${ExecName} ${SrcsList})

   # Set special properties of target
   if(${MHDForm} STREQUAL "DEFAULT")
      set_target_properties(${ExecName} PROPERTIES OUTPUT_NAME
         ${ExecName} COMPILE_FLAGS "-DQUICC_SPATIALSCHEME_${UpMHDScheme} -DQUICC_SPATIALDIMENSION_${MHDSchemeDim} -DQUICC_RUNSIM_MODEL=${MHDModel}")
   else(${MHDForm} STREQUAL "DEFAULT")
      set_target_properties(${ExecName} PROPERTIES OUTPUT_NAME
         ${ExecName} COMPILE_FLAGS "-DQUICC_SPATIALSCHEME_${UpMHDScheme} -DQUICC_SPATIALDIMENSION_${MHDSchemeDim} -DQUICC_SPATIALSCHEME_${UpMHDScheme}_${MHDForm} -DQUICC_RUNSIM_MODEL=${MHDModel}")
   endif(${MHDForm} STREQUAL "DEFAULT")

   # Show message
   message(STATUS " --> Added ${ExecName} test executable")
endfunction ()

#
# Recursive crawler to generate model list
#
function (quicc_walk_models MHDAll MHDPath MHDDirs MHDSrcRef MHDStateRef)
   # Loop over all the subdirectories
   foreach(MHDDir ${MHDDirs})
      # Set new path including subdirectory
      set(MHDNext ${MHDPath}/${MHDDir})
      # Unset list of subdirectories and sources
      set(MHDSources)
      set(MHDSrcSubDirs)
      # Include SourcesList.cmake file
      include(${MHDNext}/SourcesList.cmake)
      # Check if there are additional subdirectories
      if(DEFINED MHDSrcSubDirs)
         quicc_walk_models("${MHDAll}" "${MHDNext}" "${MHDSrcSubDirs}" "${MHDSrcRef}" "${MHDStateRef}")
      endif()
      if(MHDSources MATCHES PhysicalModel.cpp)
         set(MHDModSrcs ${${MHDAll}})
         foreach(MHDSource ${MHDSources})
            list(APPEND MHDModSrcs ${MHDNext}/${MHDSource})
         endforeach(MHDSource)
         string(REGEX REPLACE ${MHDSrcRef} "" MHDModPath ${MHDNext})
         quicc_target_model("${MHDModPath}" "${MHDStateRef}" "${MHDModSrcs}")
      else()
         # Loop over all sources and add full path source to list
         foreach(MHDSource ${MHDSources})
            list(APPEND ${MHDAll} ${MHDNext}/${MHDSource})
         endforeach(MHDSource)
         # Update the list of all sources in parent scope
         set(${MHDAll} ${${MHDAll}} PARENT_SCOPE)
      endif()
   endforeach(MHDDir)
endfunction ()

function (quicc_find_models Path ModelDir ModelState)
   set(MHDTmp ${ModelDir} ${ModelDir}/${ModelState})
   set(MHDLocal )
   quicc_walk_models(MHDLocal "${Path}" "${MHDTmp}" "${Path}" "${ModelDir}/${ModelState}")
endfunction ()

function (quicc_target_model MHDModPath MHDStateRef MHDModSrcs)
   string(REGEX REPLACE "^/" "" MHDModPath ${MHDModPath})
   # Set path to model header
   set(ModelFile ${QUICC_INCLUDE_DIR}/${MHDModPath}/PhysicalModel.hpp)

   # Extract the spatial scheme from model file
   file(STRINGS ${ModelFile} Scheme REGEX "typedef .* SchemeType;")
   STRING(REGEX REPLACE "typedef Schemes::([a-z,A-Z]*)Scheme SchemeType;" "\\1" Scheme ${Scheme})
   STRING(REGEX REPLACE " " "" Scheme ${Scheme})
   STRING(TOUPPER ${Scheme} Scheme) 
   set(QUICC_SPATIALSCHEME ${Scheme})

   # Extract the spatial scheme dimension from model file
   file(STRINGS ${ModelFile} SchemeDim REGEX "#include \"SpatialSchemes/.*/")
   STRING(REGEX REPLACE "#include \"SpatialSchemes/([1-9]*D)/.*" "\\1" SchemeDim ${SchemeDim})
   STRING(REGEX REPLACE " " "" SchemeDim ${SchemeDim})
   set(QUICC_SPATIALDIMENSION ${SchemeDim})

   # Extract the spatial scheme formulation from model file
   file(STRINGS ${ModelFile} Formulation REGEX "// QUICC_SPATIALSCHEME_FORMULATION = .*;")
   if(NOT ${Formulation} STREQUAL "")
      STRING(REGEX REPLACE "// QUICC_SPATIALSCHEME_FORMULATION = ([a-z,A-Z]*);" "\\1" Formulation ${Formulation})
      STRING(REGEX REPLACE " " "" Formulation ${Formulation})
      set(QUICC_SPATIALSCHEME_FORMULATION ${Formulation})
   else(NOT ${Formulation} STREQUAL "")
      set(Formulation "DEFAULT")
      set(QUICC_SPATIALSCHEME_FORMULATION "DEFAULT")
   endif(NOT ${Formulation} STREQUAL "")

   # Generate sources list
   set(All_Srcs )
   include(../cmake.d/setup/AllSrc.cmake)

   # Create executable name
   string(REGEX REPLACE ${MHDStateRef}/ "" PrettyName ${MHDModPath})
   message(STATUS "Generating executables for ${PrettyName}")
   # Create simulation target
   quicc_add_executable("${MHDModPath}" "${MHDStateRef}" "${Scheme}" "${SchemeDim}" "${Formulation}" "Model" "RunSimulation.cpp"
      MHDModSrcs All_Srcs)

   # Create configuration file target
   quicc_add_executable("${MHDModPath}" "${MHDStateRef}" "${Scheme}" "${SchemeDim}" "${Formulation}" "Config" "WriteConfig.cpp"
      MHDModSrcs All_Srcs)

   # Create state file target
   quicc_add_executable("${MHDModPath}" "${MHDStateRef}" "${Scheme}" "${SchemeDim}" "${Formulation}" "State" "GenerateState.cpp"
      MHDModSrcs All_Srcs)

   # Create visualization file target
   quicc_add_executable("${MHDModPath}" "${MHDStateRef}" "${Scheme}" "${SchemeDim}" "${Formulation}" "Visu" "VisualizeState.cpp"
      MHDModSrcs All_Srcs)
endfunction ()
