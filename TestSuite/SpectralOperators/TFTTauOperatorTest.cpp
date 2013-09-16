/** 
 * @file TFTTauOperatorTest.cpp
 * @brief Implementation of test cases for TFTTauOperator
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

#include "gtest/gtest.h"

#include "Framework/FrameworkMacro.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the TFTTauOperator implementation
    */
   class TFTTauOperatorTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         TFTTauOperatorTest();

         /**
          * @brief Destructor
          */
         virtual ~TFTTauOperatorTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   TFTTauOperatorTest::TFTTauOperatorTest()
   {
   }

   TFTTauOperatorTest::~TFTTauOperatorTest()
   {
   }

//   void TFTTauOperatorTest::SetUp()
//   {
//   }

//   void TFTTauOperatorTest::TearDown()
//   {
//   }

   /**
    * @brief Dummy placeholder test
    *
    * @param TFTTauOperatorTest Test fixture ID
    * @param Placeholder      Test ID
    */
   TEST_F(TFTTauOperatorTest, Placeholder)
   {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/**
 * @brief Main function to execute all test cases
 *
 * @param argc Number of arguments
 * @param argv Arguments
 */
int main(int argc, char **argv)
{
   // Initilise framework
   GeoMHDiSCC::FrameworkMacro::init();

   // Set nCpu for serial run
   int nCpu = 1;

   // Set ID and nCpu in MPI case
   #ifdef GEOMHDISCC_MPI
      // Get MPI size
      MPI_Comm_size(MPI_COMM_WORLD, &nCpu);
   #endif //GEOMHDISCC_MPI

   // Setup framework
   GeoMHDiSCC::FrameworkMacro::setup(nCpu);

   ::testing::InitGoogleTest(&argc, argv);
   int status = RUN_ALL_TESTS();

   // Finalise framework
   GeoMHDiSCC::FrameworkMacro::finalize();

   return status;
}
