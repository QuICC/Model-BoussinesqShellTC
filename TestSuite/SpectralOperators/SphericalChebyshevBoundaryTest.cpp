/** 
 * @file SphericalChebyshevBoundaryTest.cpp
 * @brief Implementation of test case for SphericalChebyshevBoundary
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

#include "gtest/gtest.h"

#include "Framework/FrameworkMacro.h"
#include "SpectralOperators/IBoundary.hpp"
#include "SpectralOperators/SphericalChebyshevBoundary.hpp"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the SphericalChebyshevBoundary implementation
    */
   class SphericalChebyshevBoundaryTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         SphericalChebyshevBoundaryTest();

         /**
          * @brief Destructor
          */
         virtual ~SphericalChebyshevBoundaryTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   SphericalChebyshevBoundaryTest::SphericalChebyshevBoundaryTest()
   {
      // Initilise framework
      FrameworkMacro::init();

      // Set nCpu for serial run
      int nCpu = 1;

      // Set ID and nCpu in MPI case
      #ifdef GEOMHDISCC_MPI
         // Get MPI size
         int size;
         MPI_Comm_size(MPI_COMM_WORLD, &nCpu);
      #endif //GEOMHDISCC_MPI

      // Setup framework
      FrameworkMacro::setup(nCpu);
   }

   SphericalChebyshevBoundaryTest::~SphericalChebyshevBoundaryTest()
   {
      // Finalise framework
      FrameworkMacro::finalize();
   }

//   void SphericalChebyshevBoundaryTest::SetUp()
//   {
//   }

//   void SphericalChebyshevBoundaryTest::TearDown()
//   {
//   }

   /**
    * @brief Test value boundary condition
    */
   TEST_F(SphericalChebyshevBoundaryTest, Value)
   {
      // basis size
      int nN = 16;

      // Create Chebyshev boundary
      Spectral::SphericalChebyshevBoundary bound(nN);

      // Create left value boundary condition (i.e x = -1)
      Array bc = bound.value(Spectral::IBoundary::LEFT);

      // Check boundary values
      EXPECT_EQ(bc(0), 0.5*static_cast<MHDFloat>(std::pow(-1,0)));
      for(int i = 1; i < nN; i++)
      {
         EXPECT_EQ(bc(i), static_cast<MHDFloat>(std::pow(-1,i)));
      }

      // Create right value boundary condition (i.e x = 1)
      bc = bound.value(Spectral::IBoundary::RIGHT);

      // Check boundary values
      EXPECT_EQ(bc(0), 0.5);
      for(int i = 1; i < nN; i++)
      {
         EXPECT_EQ(bc(i), 1.0);
      }
   }

   /**
    * @brief Test first derivative boundary condition
    */
   TEST_F(SphericalChebyshevBoundaryTest, FirstDerivative)
   {
      // basis size
      int nN = 16;

      // Create Chebyshev boundary
      Spectral::SphericalChebyshevBoundary bound(nN);

      // Create left first derivative boundary condition (i.e x = -1)
      Array bc = bound.firstDerivative(Spectral::IBoundary::LEFT);

      // Check boundary values
      for(int i = 0; i < nN; i++)
      {
         EXPECT_EQ(bc(i), static_cast<MHDFloat>(std::pow(-1,i-1)*i*i));
      }

      // Create right first derivative boundary condition (i.e x = 1)
      bc = bound.firstDerivative(Spectral::IBoundary::RIGHT);

      // Check boundary values
      for(int i = 0; i < nN; i++)
      {
         EXPECT_EQ(bc(i), static_cast<MHDFloat>(i*i));
      }
   }

   /**
    * @brief Test second derivative boundary condition
    */
   TEST_F(SphericalChebyshevBoundaryTest, SecondDerivative)
   {
      // basis size
      int nN = 16;

      // Create Chebyshev boundary
      Spectral::SphericalChebyshevBoundary bound(nN);

      // Create left first derivative boundary condition (i.e x = -1)
      Array bc = bound.secondDerivative(Spectral::IBoundary::LEFT);

      // Check boundary values
      for(int i = 0; i < nN; i++)
      {
         EXPECT_EQ(bc(i), static_cast<MHDFloat>(std::pow(-1,i)*(std::pow(i,4)-std::pow(i,2)))/3.);
      }

      // Create right first derivative boundary condition (i.e x = 1)
      bc = bound.secondDerivative(Spectral::IBoundary::RIGHT);

      // Check boundary values
      for(int i = 0; i < nN; i++)
      {
         EXPECT_EQ(bc(i), static_cast<MHDFloat>(std::pow(i,4)-std::pow(i,2))/3.);
      }
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}