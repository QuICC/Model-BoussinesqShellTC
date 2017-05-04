###################################################
#-------------- AVAILABLE PLATFORMS --------------#
###################################################

# Get list of platforms
IF(APPLE)
    file(GLOB Platforms
        ${CMAKE_SOURCE_DIR}/cmake.d/platforms/*.cmake)
ELSE()
    file(GLOB Platforms
        ${CMAKE_SOURCE_DIR}/cmake.d/platforms/[A-Z]*.cmake)
ENDIF()

# Cleanup list
foreach(Platform ${Platforms})
   # Extract platform name
   get_filename_component(PlatformName ${Platform} NAME_WE)

   set(QUICC_PLATFORMS ${PlatformName} ${QUICC_PLATFORMS})
endforeach(Platform ${Platforms})
