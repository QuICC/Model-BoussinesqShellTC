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

   # Create list of source files
   set(SrcsList RunSimulation.cpp ${GEOMHDISCC_SRC_DIR}/PhysicalModels/${ModelName}.cpp ${All_Srcs})

   # Add executable to target list
   add_executable(${ModelName} ${SrcsList})

   # Set special properties of target
   set_target_properties(${ModelName} PROPERTIES OUTPUT_NAME ${ModelName}
      COMPILE_FLAGS "-D GEOMHDISCC_RUNSIM_MODEL=${ModelName}")

   # Create list of source files
   set(SrcsList WriteConfig.cpp ${GEOMHDISCC_SRC_DIR}/PhysicalModels/${ModelName}.cpp ${All_Srcs})

   # Add executable to target list
   add_executable(${ModelName}Config ${SrcsList})

   # Set special properties of target
   set_target_properties(${ModelName}Config PROPERTIES OUTPUT_NAME
      ${ModelName}Config COMPILE_FLAGS "-D GEOMHDISCC_RUNSIM_MODEL=${ModelName}")

   # Create list of source files
   set(SrcsList GenerateState.cpp ${GEOMHDISCC_SRC_DIR}/PhysicalModels/${ModelName}.cpp ${All_Srcs})

   # Add executable to target list
   add_executable(${ModelName}State ${SrcsList})

   # Set special properties of target
   set_target_properties(${ModelName}State PROPERTIES OUTPUT_NAME
      ${ModelName}State COMPILE_FLAGS "-D GEOMHDISCC_RUNSIM_MODEL=${ModelName}")

   # Show message
   message(STATUS " --> Added ${ModelName}Config executable")
endforeach(ModelFile)
