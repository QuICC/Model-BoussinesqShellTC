#####################################
######## SpatialSimulation ##########
#####################################

#
# Add SpatialSimulation unit test
#
set(SphUnitTests ${SphUnitTests} unit_SpatialSimulation)

#
# Create all the tests
#
foreach(UnitTest ${SphUnitTests})
   set(SrcsList TestSuite/Spherical/${UnitTest} ${All_Srcs})
   get_filename_component(unit ${UnitTest} NAME)
   add_executable(${unit} ${SrcsList})
endforeach(UnitTest)
