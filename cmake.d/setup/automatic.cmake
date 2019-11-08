###################################################
#--------------- COMPILER SETTINGS ---------------#
###################################################

#
# General compiler Settings
#
set(CMAKE_CXX_FLAGS_RELEASE "-DEIGEN_NO_DEBUG -DNDEBUG" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-DEIGEN_NO_DEBUG -DNDEBUG -g" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall" CACHE STRING "" FORCE)
if(QUICC_DISABLE_RDYNAMIC)
   set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "" CACHE STRING "" FORCE)
   mark_as_advanced(FORCE CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
endif(QUICC_DISABLE_RDYNAMIC)
if(QUICC_ENABLE_DYNAMIC)
   set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "-dynamic ${CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS}" CACHE STRING "" FORCE)
   mark_as_advanced(FORCE CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
endif(QUICC_ENABLE_DYNAMIC)

string(TOUPPER "${QUICC_COMPILER}" QUICC_COMPILER)

if(QUICC_COMPILER STREQUAL "SCALASCA")
   set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE scorep)
   set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK scorep)
endif(QUICC_COMPILER STREQUAL "SCALASCA")

if(QUICC_MPI)
   set(CMAKE_CXX_COMPILER ${QUICC_CC_MPI_${QUICC_COMPILER}})
else(QUICC_MPI)
   set(CMAKE_CXX_COMPILER ${QUICC_CC_SERIAL_${QUICC_COMPILER}})
endif(QUICC_MPI)

set(CMAKE_CXX_FLAGS "${QUICC_CC_ARCH_${QUICC_COMPILER}}" CACHE STRING "" FORCE)
if(NOT QUICC_THREADS STREQUAL "None")
   string(TOUPPER "${QUICC_THREADS}" model)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${QUICC_CC_${model}_${QUICC_COMPILER}}" CACHE STRING "" FORCE)
endif(NOT QUICC_THREADS STREQUAL "None")

if(QUICC_MPI)
   foreach(inc ${QUICC_CC_INC_MPI_${QUICC_COMPILER}})
      include_directories(${inc})
   endforeach(inc)
   foreach(lib ${QUICC_CC_LIB_MPI_${QUICC_COMPILER}})
      link_libraries(${lib})
   endforeach(lib)
else(QUICC_MPI)
   foreach(inc ${QUICC_CC_INC_${QUICC_COMPILER}})
      include_directories(${inc})
   endforeach(inc)
   foreach(lib ${QUICC_CC_LIB_${QUICC_COMPILER}})
      link_libraries(${lib})
   endforeach(lib)
endif(QUICC_MPI)

# General libraries and includes
foreach(lib ${QUICC_LIBRARIES})
   link_libraries(${lib})
endforeach(lib)
foreach(inc ${QUICC_INCLUDES})
   include_directories(${inc})
endforeach(inc)

# Smart pointers libraries and includes
quicc_link_external(QUICC_SMARTPTR)

# FFT implementation libraries and includes
quicc_link_external(QUICC_FFT QUICC_THREADS)

# Linear algebra libraries and includes
quicc_link_external(QUICC_LINALG)

# Sparse linear algebra libraries and includes
quicc_link_external(QUICC_SPLINALG)

# Sparse SPD linear algebra libraries and includes
quicc_link_external(QUICC_SPSPDLINALG)

# Large IO format libraries and includes
quicc_link_external(QUICC_LARGEIO)

# Multiple precision libraries and includes
quicc_link_external(QUICC_MPBACKEND)

# Multiple precision libraries and includes
quicc_link_external(QUICC_PYTHON)
