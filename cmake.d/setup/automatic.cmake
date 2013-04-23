###################################################
#--------------- COMPILER SETTINGS ---------------#
###################################################

#
# General compiler Settings
#
set(CMAKE_CXX_FLAGS_RELEASE "-DEIGEN_NO_DEBUG -DNDEBUG" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall" CACHE STRING "" FORCE)

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
geomhdiscc_link_external(GEOMHDISCC_SMARTPTR)

# FFT implementation libraries and includes
geomhdiscc_link_external(GEOMHDISCC_FFT)

# Linear algebra libraries and includes
geomhdiscc_link_external(GEOMHDISCC_LINALG)

# Sparse linear algebra libraries and includes
geomhdiscc_link_external(GEOMHDISCC_SPLINALG)

# Sparse eigen solver libraries and includes
geomhdiscc_link_external(GEOMHDISCC_SPEIGSOLVER)

# Large IO format libraries and includes
geomhdiscc_link_external(GEOMHDISCC_LARGEIO)

# Multiple precision libraries and includes
geomhdiscc_link_external(GEOMHDISCC_MPLIB)
