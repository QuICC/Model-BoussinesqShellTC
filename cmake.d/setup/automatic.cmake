###################################################
#--------------- COMPILER SETTINGS ---------------#
###################################################

#
# General compiler Settings
#
set(CMAKE_CXX_FLAGS_RELEASE "-DEIGEN_NO_DEBUG -DNDEBUG" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-DEIGEN_NO_DEBUG -DNDEBUG -g" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall" CACHE STRING "" FORCE)
if(GEOMHDISCC_DISABLE_RDYNAMIC)
   set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "" CACHE STRING "" FORCE)
   mark_as_advanced(FORCE CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
endif(GEOMHDISCC_DISABLE_RDYNAMIC)
if(GEOMHDISCC_ENABLE_DYNAMIC)
   set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "-dynamic ${CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS}" CACHE STRING "" FORCE)
   mark_as_advanced(FORCE CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
endif(GEOMHDISCC_ENABLE_DYNAMIC)

string(TOUPPER "${GEOMHDISCC_COMPILER}" GEOMHDISCC_COMPILER)

if(GEOMHDISCC_COMPILER STREQUAL "SCALASCA")
   set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE scorep)
   set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK scorep)
endif(GEOMHDISCC_COMPILER STREQUAL "SCALASCA")

if(GEOMHDISCC_MPI)
   set(CMAKE_CXX_COMPILER ${GEOMHDISCC_CC_MPI_${GEOMHDISCC_COMPILER}})
else(GEOMHDISCC_MPI)
   set(CMAKE_CXX_COMPILER ${GEOMHDISCC_CC_SERIAL_${GEOMHDISCC_COMPILER}})
endif(GEOMHDISCC_MPI)

set(CMAKE_CXX_FLAGS "${GEOMHDISCC_CC_ARCH_${GEOMHDISCC_COMPILER}}" CACHE STRING "" FORCE)
if(NOT GEOMHDISCC_THREADS STREQUAL "None")
   string(TOUPPER "${GEOMHDISCC_THREADS}" model)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GEOMHDISCC_CC_${model}_${GEOMHDISCC_COMPILER}}" CACHE STRING "" FORCE)
endif(NOT GEOMHDISCC_THREADS STREQUAL "None")

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
quicc_link_external(GEOMHDISCC_SMARTPTR)

# FFT implementation libraries and includes
quicc_link_external(GEOMHDISCC_FFT GEOMHDISCC_THREADS)

# Linear algebra libraries and includes
quicc_link_external(GEOMHDISCC_LINALG)

# Sparse linear algebra libraries and includes
quicc_link_external(GEOMHDISCC_SPLINALG)

# Sparse SPD linear algebra libraries and includes
quicc_link_external(GEOMHDISCC_SPSPDLINALG)

# Large IO format libraries and includes
quicc_link_external(GEOMHDISCC_LARGEIO)

# Multiple precision libraries and includes
quicc_link_external(GEOMHDISCC_MPLIB)

# Multiple precision libraries and includes
quicc_link_external(GEOMHDISCC_PYTHON)
