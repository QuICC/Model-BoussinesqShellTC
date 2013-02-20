###################################################
#--------------- COMPILER SETTINGS ---------------#
###################################################

#
# General compiler Settings
#
set(CMAKE_CXX_FLAGS_RELEASE "-DEIGEN_NO_DEBUG -DNDEBUG" CACHE STRING "" FORCE)

if(GEOMHDISCC_MPI)
   set(CMAKE_CXX_COMPILER ${GEOMHDISCC_CC_MPI_${GEOMHDISCC_COMPILER}})
else(GEOMHDISCC_MPI)
   set(CMAKE_CXX_COMPILER ${GEOMHDISCC_CC_SERIAL_${GEOMHDISCC_COMPILER}})
endif(GEOMHDISCC_MPI)

set(CMAKE_CXX_FLAGS "${GEOMHDISCC_CC_ARCH_${GEOMHDISCC_COMPILER}}" CACHE STRING "" FORCE)

if(GEOMHDISCC_MPI)
   foreach(inc ${GEOMHDISCC_CC_INC_MPI_${GEOMHDISCC_COMPILER}})
      include_directories(${inc})
   endforeach(inc)
   foreach(lib ${GEOMHDISCC_CC_LIB_MPI_${GEOMHDISCC_COMPILER}})
      link_libraries(${lib})
   endforeach(lib)
else(GEOMHDISCC_MPI)
   foreach(inc ${GEOMHDISCC_CC_INC_${GEOMHDISCC_COMPILER}})
      include_directories(${inc})
   endforeach(inc)
   foreach(lib ${GEOMHDISCC_CC_LIB_${GEOMHDISCC_COMPILER}})
      link_libraries(${lib})
   endforeach(lib)
endif(GEOMHDISCC_MPI)

# General libraries and includes
foreach(lib ${GEOMHDISCC_LIBRARIES})
   link_libraries(${lib})
endforeach(lib)
foreach(inc ${GEOMHDISCC_INCLUDES})
   include_directories(${inc})
endforeach(inc)

# Smart pointers libraries and includes
foreach(lib ${GEOMHDISCC_LIBRARIES_${GEOMHDISCC_SMARTPTR}})
   link_libraries(${lib})
endforeach(lib)
foreach(inc ${GEOMHDISCC_INCLUDES_${GEOMHDISCC_SMARTPTR}})
   include_directories(${inc})
endforeach(inc)

# FFT implementation libraries and includes
foreach(lib ${GEOMHDISCC_LIBRARIES_${GEOMHDISCC_FFT}})
   link_libraries(${lib})
endforeach(lib)
foreach(inc ${GEOMHDISCC_INCLUDES_${GEOMHDISCC_FFT}})
   include_directories(${inc})
endforeach(inc)

# Linear algebra libraries and includes
foreach(lib ${GEOMHDISCC_LIBRARIES_${GEOMHDISCC_LINALG}})
   link_libraries(${lib})
endforeach(lib)
foreach(inc ${GEOMHDISCC_INCLUDES_${GEOMHDISCC_LINALG}})
   include_directories(${inc})
endforeach(inc)

# Sparse linear algebra libraries and includes
foreach(lib ${GEOMHDISCC_LIBRARIES_${GEOMHDISCC_SPLINALG}})
   link_libraries(${lib})
endforeach(lib)
foreach(inc ${GEOMHDISCC_INCLUDES_${GEOMHDISCC_SPLINALG}})
   include_directories(${inc})
endforeach(inc)

# Large IO format libraries and includes
foreach(lib ${GEOMHDISCC_LIBRARIES_${GEOMHDISCC_LARGEIO}})
   link_libraries(${lib})
endforeach(lib)
foreach(inc ${GEOMHDISCC_INCLUDES_${GEOMHDISCC_LARGEIO}})
   include_directories(${inc})
endforeach(inc)

# Multiple precision libraries and includes
foreach(lib ${GEOMHDISCC_LIBRARIES_${GEOMHDISCC_MPLIB}})
   link_libraries(${lib})
endforeach(lib)
foreach(inc ${GEOMHDISCC_INCLUDES_${GEOMHDISCC_MPLIB}})
   include_directories(${inc})
endforeach(inc)
