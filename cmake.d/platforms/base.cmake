###################################################
#-------------- AVAILABLE PLATFORMS --------------#
###################################################

# Get list of platforms
# Leo: this doesn't work 

IF(APPLE)
    file(GLOB Platforms
        ${CMAKE_SOURCE_DIR}/cmake.d/platforms/*.cmake)
ELSE()
    file(GLOB Platforms
        ${CMAKE_SOURCE_DIR}/cmake.d/platforms/[A-Z]*.cmake)
ENDIF()

    #IF (APPLE)
    #ENDIF (APPLE)


# Cleanup list
foreach(Platform ${Platforms})
   # Extract platform name
   get_filename_component(PlatformName ${Platform} NAME_WE)

   set(QUICC_PLATFORMS ${PlatformName} ${QUICC_PLATFORMS})
endforeach(Platform ${Platforms})
