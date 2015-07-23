#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%% Create all the targets with a general implementation %%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Get list of physical models and test models
include(${GEOMHDISCC_INCLUDE_DIR}/PhysicalModels/cmake.d/ModelsList.cmake)
geomhdiscc_add_path(MHDPhysicalModels ${GEOMHDISCC_INCLUDE_DIR}/PhysicalModels)
geomhdiscc_add_path(MHDTestModels ${GEOMHDISCC_INCLUDE_DIR}/PhysicalModels)

# Loop over all physical models
foreach(ModelFile ${MHDPhysicalModels})

   # Extract filename from full path
   get_filename_component(ModelName ${ModelFile} NAME_WE)

   # Extract the spatial scheme from model file
   file(STRINGS ${ModelFile} Scheme REGEX "typedef .* SchemeType;")
   STRING(REGEX REPLACE "typedef Schemes::([a-z,A-Z]*)Scheme SchemeType;" "\\1" Scheme ${Scheme})
   STRING(REGEX REPLACE " " "" Scheme ${Scheme})
   STRING(TOUPPER ${Scheme} Scheme) 
   set(GEOMHDISCC_SPATIALSCHEME ${Scheme})

   # Extract the spatial scheme dimension from model file
   file(STRINGS ${ModelFile} SchemeDim REGEX "#include \"SpatialSchemes/.*/")
   STRING(REGEX REPLACE "#include \"SpatialSchemes/([1-9]*D)/.*" "\\1" SchemeDim ${SchemeDim})
   STRING(REGEX REPLACE " " "" SchemeDim ${SchemeDim})
   set(GEOMHDISCC_SPATIALDIMENSION ${SchemeDim})

   # Extract the spatial scheme formulation from model file
   file(STRINGS ${ModelFile} Formulation REGEX "// GEOMHDISCC_SPATIALSCHEME_FORMULATION = .*;")
   if(NOT ${Formulation} STREQUAL "")
      STRING(REGEX REPLACE "// GEOMHDISCC_SPATIALSCHEME_FORMULATION = ([a-z,A-Z]*);" "\\1" Formulation ${Formulation})
      STRING(REGEX REPLACE " " "" Formulation ${Formulation})
      set(GEOMHDISCC_SPATIALSCHEME_FORMULATION ${Formulation})
   else(NOT ${Formulation} STREQUAL "")
      set(Formulation "DEFAULT")
      set(GEOMHDISCC_SPATIALSCHEME_FORMULATION "DEFAULT")
   endif(NOT ${Formulation} STREQUAL "")

   # Generate sources list
   set(All_Srcs )
   include(../cmake.d/setup/AllSrc.cmake)

   # Get model specific sources
   include(${GEOMHDISCC_INCLUDE_DIR}/PhysicalModels/cmake.d/${ModelName}Sources.cmake)

   # Create simulation target
   geomhdiscc_add_executable(${ModelName} ${Scheme} ${SchemeDim} ${Formulation} "Model" "RunSimulation.cpp"
      MHDModelSources All_Srcs)

   # Create configuration file target
   geomhdiscc_add_executable(${ModelName} ${Scheme} ${SchemeDim} ${Formulation} "Config" "WriteConfig.cpp"
      MHDModelSources All_Srcs)

   # Create state file target
   geomhdiscc_add_executable(${ModelName} ${Scheme} ${SchemeDim} ${Formulation} "State" "GenerateState.cpp"
      MHDModelSources All_Srcs)

   # Create visualization file target
   geomhdiscc_add_executable(${ModelName} ${Scheme} ${SchemeDim} ${Formulation} "Visu" "VisualizeState.cpp"
      MHDModelSources All_Srcs)
endforeach(ModelFile)

# Loop over all test models
message(STATUS "***********************************************")
message(STATUS "***************** Test models *****************")
message(STATUS "***********************************************")
foreach(ModelFile ${MHDTestModels})

   # Extract filename from full path
   get_filename_component(ModelName ${ModelFile} NAME_WE)

   # Extract the spatial scheme from model file
   file(STRINGS ${ModelFile} Scheme REGEX "typedef .* SchemeType;")
   STRING(REGEX REPLACE "typedef Schemes::([a-z,A-Z]*)Scheme SchemeType;" "\\1" Scheme ${Scheme})
   STRING(REGEX REPLACE " " "" Scheme ${Scheme})
   STRING(TOUPPER ${Scheme} Scheme) 
   set(GEOMHDISCC_SPATIALSCHEME ${Scheme})

   # Extract the spatial scheme dimension from model file
   file(STRINGS ${ModelFile} SchemeDim REGEX "#include \"SpatialSchemes/.*/")
   STRING(REGEX REPLACE "#include \"SpatialSchemes/([1-9]*D)/.*" "\\1" SchemeDim ${SchemeDim})
   STRING(REGEX REPLACE " " "" SchemeDim ${SchemeDim})
   set(GEOMHDISCC_SPATIALDIMENSION ${SchemeDim})

   # Extract the spatial scheme formulation from model file
   file(STRINGS ${ModelFile} Formulation REGEX "// GEOMHDISCC_SPATIALSCHEME_FORMULATION = .*;")
   if(NOT ${Formulation} STREQUAL "")
      STRING(REGEX REPLACE "// GEOMHDISCC_SPATIALSCHEME_FORMULATION = ([a-z,A-Z]*);" "\\1" Formulation ${Formulation})
      STRING(REGEX REPLACE " " "" Formulation ${Formulation})
      set(GEOMHDISCC_SPATIALSCHEME_FORMULATION ${Formulation})
   else(NOT ${Formulation} STREQUAL "")
      set(Formulation "DEFAULT")
      set(GEOMHDISCC_SPATIALSCHEME_FORMULATION "DEFAULT")
   endif(NOT ${Formulation} STREQUAL "")

   # Generate sources list
   set(All_Srcs )
   include(../cmake.d/setup/AllSrc.cmake)

   # Get model specific sources
   include(${GEOMHDISCC_INCLUDE_DIR}/PhysicalModels/cmake.d/${ModelName}Sources.cmake)

   # Create simulation target
   geomhdiscc_add_test(${ModelName} ${Scheme} ${SchemeDim} ${Formulation} "Test" "RunSimulation.cpp"
      MHDModelSources All_Srcs)

   # Create configuration file target
   geomhdiscc_add_test(${ModelName} ${Scheme} ${SchemeDim} ${Formulation} "Config" "WriteConfig.cpp"
      MHDModelSources All_Srcs)

   # Create state file target
   geomhdiscc_add_test(${ModelName} ${Scheme} ${SchemeDim} ${Formulation} "State" "GenerateState.cpp"
      MHDModelSources All_Srcs)

   # Create visualization file target
   geomhdiscc_add_test(${ModelName} ${Scheme} ${SchemeDim} ${Formulation} "Visu" "VisualizeState.cpp"
      MHDModelSources All_Srcs)
endforeach(ModelFile)
