###################################################
#----------- TESTSUITE CONFIGURATION -------------#
###################################################

#
# Provide tests
#
option(QUICC_TEST "Provide tests?" OFF)

if(QUICC_TEST)
   message(STATUS "***********************************************")
   message(STATUS "*************** TestSuite setup ***************")
   message(STATUS "***********************************************")

   add_subdirectory("TestSuite" EXCLUDE_FROM_ALL)
endif(QUICC_TEST)
