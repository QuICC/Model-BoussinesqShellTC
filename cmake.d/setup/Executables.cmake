###################################################
#----------------- Executables -------------------#
###################################################

message(STATUS "***********************************************")
message(STATUS "************** Executables setup **************")
message(STATUS "***********************************************")


###################################################
#---------- ACTIVATE DEVLOPMENT MODELS -----------#
###################################################

option(QUICC_DEVEL_MODEL "Enable development models?" OFF)
mark_as_advanced(FORCE QUICC_DEVEL_MODEL)

add_subdirectory("Executables" EXCLUDE_FROM_ALL)
