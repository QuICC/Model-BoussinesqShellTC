#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%% Create all the targets with a general implementation %%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Get list of physical models
file(GLOB PhysicalModels
   ${GEOMHDISCC_INCLUDE_DIR}/PhysicalModels/*Model.hpp)

# Loop over all physical models
foreach(ModelFile ${PhysicalModels})

   # Extract filename from full path
   get_filename_component(ModelName ${ModelFile} NAME_WE)

   # Extract the spatial scheme from model file
   file(STRINGS ${ModelFile} Scheme REGEX "typedef .* SchemeType;")
   STRING(REGEX REPLACE "typedef Schemes::([a-z,A-Z]*)Scheme SchemeType;" "\\1" Scheme ${Scheme})
   STRING(REGEX REPLACE " " "" Scheme ${Scheme})

   # Get model specific sources
   include(${GEOMHDISCC_INCLUDE_DIR}/PhysicalModels/cmake.d/${ModelName}Sources.cmake)

   # Create list of source files
   set(SrcsList RunSimulation.cpp ${GEOMHDISCC_SRC_DIR}/PhysicalModels/${ModelName}.cpp ${All_Srcs} ${MHDModelSources})

   # Add executable to target list
   add_executable(${ModelName} ${SrcsList})

   # Set special properties of target
   set_target_properties(${ModelName} PROPERTIES OUTPUT_NAME ${ModelName}
      COMPILE_FLAGS "-DGEOMHDISCC_SPATIALSCHEME_${Scheme} -DGEOMHDISCC_RUNSIM_MODEL=${ModelName}")

   # Show message
   message(STATUS " --> Added ${ModelName} executable")

   # Create list of source files
   set(SrcsList WriteConfig.cpp ${GEOMHDISCC_SRC_DIR}/PhysicalModels/${ModelName}.cpp ${All_Srcs} ${MHDModelSources})

   # Add executable to target list
   add_executable(${ModelName}Config ${SrcsList})

   # Set special properties of target
   set_target_properties(${ModelName}Config PROPERTIES OUTPUT_NAME
      ${ModelName}Config COMPILE_FLAGS "-DGEOMHDISCC_SPATIALSCHEME_${Scheme} -DGEOMHDISCC_RUNSIM_MODEL=${ModelName}")

   # Show message
   message(STATUS " --> Added ${ModelName}Config executable")

   # Create list of source files
   set(SrcsList GenerateState.cpp ${GEOMHDISCC_SRC_DIR}/PhysicalModels/${ModelName}.cpp ${All_Srcs} ${MHDModelSources})

   # Add executable to target list
   add_executable(${ModelName}State ${SrcsList})

   # Set special properties of target
   set_target_properties(${ModelName}State PROPERTIES OUTPUT_NAME
      ${ModelName}State COMPILE_FLAGS "-DGEOMHDISCC_SPATIALSCHEME_${Scheme} -DGEOMHDISCC_RUNSIM_MODEL=${ModelName}")

   # Show message
   message(STATUS " --> Added ${ModelName}State executable")

   # Create list of source files
   set(SrcsList VisualizeState.cpp ${GEOMHDISCC_SRC_DIR}/PhysicalModels/${ModelName}.cpp ${All_Srcs} ${MHDModelSources})

   # Add executable to target list
   add_executable(${ModelName}Visu ${SrcsList})

   # Set special properties of target
   set_target_properties(${ModelName}Visu PROPERTIES OUTPUT_NAME
      ${ModelName}Visu COMPILE_FLAGS "-DGEOMHDISCC_SPATIALSCHEME_${Scheme} -DGEOMHDISCC_RUNSIM_MODEL=${ModelName}")

   # Show message
   message(STATUS " --> Added ${ModelName}Visu executable")
endforeach(ModelFile)
