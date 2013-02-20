#
# Add ScalarField1D unit test
#
set(UnitTests ${UnitTests} unit_ScalarField1D)

#
# Add ScalarField2D unit test
#
set(UnitTests ${UnitTests} unit_ScalarField2D)

#
# Add ScalarField3D unit test
#
set(UnitTests ${UnitTests} unit_ScalarField3D)

#
# Add Transform1DCoordinator unit test
#
set(UnitTests ${UnitTests} unit_Transform1DCoordinator)

#
# Add Transform2DCoordinator unit test
#
set(UnitTests ${UnitTests} unit_Transform2DCoordinator)

#
# Add Transform3DCoordinator unit test
#
set(UnitTests ${UnitTests} unit_Transform3DCoordinator)

#
# Add ConfigurationFile unit test
#
set(UnitTests ${UnitTests} unit_ConfigurationFile)

#
# Add QuasiInverse unit test
#
set(UnitTests ${UnitTests} unit_QuasiInverse)


#####################################
########## LoadSplitters ############
#####################################

#
# Add LoadSplitter unit test
#
set(UnitTests ${UnitTests} "LoadSplitters/unit_LoadSplitter")


#####################################
############ Polynomials ############
#####################################

#
# Add AssociatedLegendreBasis unit test
#
set(UnitTests ${UnitTests} "Polynomials/unit_AssociatedLegendreBasis")

#
# Add ChebyshevBasis unit test
#
set(UnitTests ${UnitTests} "Polynomials/unit_ChebyshevBasis")

#
# Add JacobiBasis unit test
#
set(UnitTests ${UnitTests} "Polynomials/unit_JacobiBasis")

#
# Add OnesidedChebyshevBasis unit test
#
set(UnitTests ${UnitTests} "Polynomials/unit_OnesidedChebyshevBasis")

#
# Add OnesidedJacobiBasis unit test
#
set(UnitTests ${UnitTests} "Polynomials/unit_OnesidedJacobiBasis")

#
# Add WorlandBasis unit test
#
set(UnitTests ${UnitTests} "Polynomials/unit_WorlandBasis")


#####################################
############ Quadratures ############
#####################################

#
# Add ChebyshevRule unit test
#
set(UnitTests ${UnitTests} "Quadratures/unit_ChebyshevRule")

#
# Add LegendreRule unit test
#
set(UnitTests ${UnitTests} "Quadratures/unit_LegendreRule")

#
# Add JacobiRule unit test
#
set(UnitTests ${UnitTests} "Quadratures/unit_JacobiRule")

#
# Create all the tests
#
foreach(UnitTest ${UnitTests})
   set(SrcsList TestSuite/Base/${UnitTest} ${All_Srcs})
   get_filename_component(unit ${UnitTest} NAME)
   add_executable(${unit} ${SrcsList})
endforeach(UnitTest)
