/** 
 * @file AssociatedLegendreTransformTest.cpp
 * @brief Implementation of test cases for Associated Legendre transforms
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#include "gtest/gtest.h"

#include "Framework/FrameworkMacro.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Googletest fixture for the AssociatedLegendreTransform implementation
    */
   class AssociatedLengedreTransformTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         AssociatedLegendreTransformTest();

         /**
          * @brief Destructor
          */
         virtual ~AssociatedLegenreTransformTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
         
         /// Acceptable error
         double mError;
   };

   AssociatedLegendreTransformTest::AssociatedLegendreTransformTest()
      : mError(1e-13)
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

   AssociatedLegendreTransformTest::~AssociatedLegendreTransformTest()
   {
      // Finalise framework
      FrameworkMacro::finalize();
   }

   /**
    * @brief Test mesh grid
    *
    * @param AssociatedLegendreTransformTest Test fixture ID
    * @param MeshGrid                        Test ID
    */
   TEST_F(AssociatedLegendreTransformTest, MeshGrid)
   {
   }

}
}

/**
 * @brief Main function to execute all test cases
 *
 * @param argc Number of arguments
 * @param argv Arguments
 */
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
