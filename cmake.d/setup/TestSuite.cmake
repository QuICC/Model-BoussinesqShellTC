###################################################
#----------- TESTSUITE CONFIGURATION -------------#
###################################################

#
# Install tests
#
INSTALL(DIRECTORY TestSuite/Bash DESTINATION ${PROJECT_BINARY_DIR}/TestSuite FILES_MATCHING PATTERN "*.sh")
