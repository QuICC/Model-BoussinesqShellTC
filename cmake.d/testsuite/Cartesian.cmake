#####################################
######## SpatialSimulation ##########
#####################################

#
# Add SpatialSimulation unit test
#
set(CarUnitTests ${CarUnitTests} unit_SpatialSimulation)

#
# Create all the tests
#
foreach(UnitTest ${CarUnitTests})
   set(SrcsList TestSuite/Cartesian/${UnitTest} ${All_Srcs})
   get_filename_component(unit ${UnitTest} NAME)
   add_executable(${unit} ${SrcsList})
endforeach(UnitTest)
