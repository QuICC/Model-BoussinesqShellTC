###################################################
#----------- TESTSUITE CONFIGURATION -------------#
###################################################

#
# Provide tests
#
option(GEOMHDISCC_TEST "Provide tests?" OFF)

if(GEOMHDISCC_TEST)
   message(STATUS "***********************************************")
   message(STATUS "*************** TestSuite setup ***************")
   message(STATUS "***********************************************")

   add_subdirectory("TestSuite" EXCLUDE_FROM_ALL)
endif(GEOMHDISCC_TEST)
