###################################################
#----------------- PYTHON SCRIPTS ----------------#
###################################################

#
# Install Scripts package
#
INSTALL(DIRECTORY Scripts/Bash DESTINATION ${PROJECT_BINARY_DIR}/Scripts FILES_MATCHING PATTERN "*.sh")
INSTALL(DIRECTORY Scripts/Python DESTINATION ${PROJECT_BINARY_DIR}/Scripts FILES_MATCHING PATTERN "*.py")
INSTALL(DIRECTORY Scripts/Slurm DESTINATION ${PROJECT_BINARY_DIR}/Scripts FILES_MATCHING PATTERN "*.slurm")
