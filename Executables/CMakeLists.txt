#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%% Create all the targets with a general implementation %%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Get the stable models
quicc_find_models(${QUICC_SRC_DIR} "Model" "Stable")
# Get the staging models
quicc_find_models(${QUICC_SRC_DIR} "Model" "Staging")

if(QUICC_DEVEL_MODEL)
   # Get the devel models
   quicc_find_models(${QUICC_SRC_DIR} "Model" "Devel")
endif(QUICC_DEVEL_MODEL)
