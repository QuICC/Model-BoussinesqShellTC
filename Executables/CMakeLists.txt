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
   set(GEOMHDISCC_SPATIALSCHEME ${Scheme})

   # Generate sources list
   set(All_Srcs )
   include(../cmake.d/setup/AllSrc.cmake)

   # Get model specific sources
   include(${GEOMHDISCC_INCLUDE_DIR}/PhysicalModels/cmake.d/${ModelName}Sources.cmake)

   # Create simulation target
   geomhdiscc_add_executable(${ModelName} ${Scheme} "Model" "RunSimulation.cpp"
      MHDModelSources All_Srcs)

   # Create configuration file target
   geomhdiscc_add_executable(${ModelName} ${Scheme} "Config" "WriteConfig.cpp"
      MHDModelSources All_Srcs)

   # Create state file target
   geomhdiscc_add_executable(${ModelName} ${Scheme} "State" "GenerateState.cpp"
      MHDModelSources All_Srcs)

   # Create visualization file target
   geomhdiscc_add_executable(${ModelName} ${Scheme} "Visu" "VisualizeState.cpp"
      MHDModelSources All_Srcs)
endforeach(ModelFile)

# Loop over all test models
foreach(ModelFile ${MHDTestModels})

   # Extract filename from full path
   get_filename_component(ModelName ${ModelFile} NAME_WE)

   # Extract the spatial scheme from model file
   file(STRINGS ${ModelFile} Scheme REGEX "typedef .* SchemeType;")
   STRING(REGEX REPLACE "typedef Schemes::([a-z,A-Z]*)Scheme SchemeType;" "\\1" Scheme ${Scheme})
   STRING(REGEX REPLACE " " "" Scheme ${Scheme})

   # Get model specific sources
   include(${GEOMHDISCC_INCLUDE_DIR}/PhysicalModels/cmake.d/${ModelName}Sources.cmake)

   # Create simulation target
   geomhdiscc_add_executable(${ModelName} ${Scheme} "Model" "RunSimulation.cpp"
      MHDModelSources All_Srcs)

   # Create configuration file target
   geomhdiscc_add_executable(${ModelName} ${Scheme} "Config" "WriteConfig.cpp"
      MHDModelSources All_Srcs)

   # Create state file target
   geomhdiscc_add_executable(${ModelName} ${Scheme} "State" "GenerateState.cpp"
      MHDModelSources All_Srcs)

   # Create visualization file target
   geomhdiscc_add_executable(${ModelName} ${Scheme} "Visu" "VisualizeState.cpp"
      MHDModelSources All_Srcs)
endforeach(ModelFile)
