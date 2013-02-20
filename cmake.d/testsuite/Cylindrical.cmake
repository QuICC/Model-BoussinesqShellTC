#####################################
######## SpatialSimulation ##########
#####################################

#
# Add SpatialSimulation unit test
#
set(CylUnitTests ${CylUnitTests} unit_SpatialSimulation)

#
# Create all the tests
#
foreach(UnitTest ${CylUnitTests})
   set(SrcsList TestSuite/Cylindrical/${UnitTest} ${All_Srcs})
   get_filename_component(unit ${UnitTest} NAME)
   add_executable(${unit} ${SrcsList})
endforeach(UnitTest)
