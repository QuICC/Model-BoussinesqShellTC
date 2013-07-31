###################################################
#-------------- AVAILABLE PLATFORMS --------------#
###################################################

# Get list of platforms
file(GLOB Platforms
   ${CMAKE_SOURCE_DIR}/cmake.d/platforms/[A-Z]*.cmake)

# Cleanup list
foreach(Platform ${Platforms})
   # Extract platform name
   get_filename_component(PlatformName ${Platform} NAME_WE)

   set(GEOMHDISCC_PLATFORMS ${PlatformName} ${GEOMHDISCC_PLATFORMS})
endforeach(Platform ${Platforms})
