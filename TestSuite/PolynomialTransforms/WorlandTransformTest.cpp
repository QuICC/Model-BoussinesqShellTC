/** 
 * @file WorlandTransformTest.cpp
 * @brief Implementation of test cases for Worland transforms
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

#include "gtest/gtest.h"

#include "Framework/FrameworkMacro.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Googletest fixture for the WorlandTransform implementation
    */
   class WorlandTransformTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         WorlandTransformTest();

         /**
          * @brief Destructor
          */
         virtual ~WorlandTransformTest();

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

   WorlandTransformTest::WorlandTransformTest()
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

   WorlandTransformTest::~WorlandTransformTest()
   {
      // Finalise framework
      FrameworkMacro::finalize();
   }

   /**
    * @brief Test mesh grid
    *
    * @param WorlandTransformTest Test fixture ID
    * @param MeshGrid                        Test ID
    */
   TEST_F(WorlandTransformTest, MeshGrid)
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
