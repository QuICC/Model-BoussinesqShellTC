/** 
 * @file MpiTimerTest.cpp
 * @brief Implementation of test cases for the MPI timer
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#include "gtest/gtest.h"

#include "Framework/FrameworkMacro.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the MpiTimer implementation
    */
   class MpiTimerTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         MpiTimerTest();

         /**
          * @brief Destructor
          */
         virtual ~MpiTimerTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   MpiTimerTest::MpiTimerTest()
   {
   }

   MpiTimerTest::~MpiTimerTest()
   {
   }

//   void MpiTimerTest::SetUp()
//   {
//   }

//   void MpiTimerTest::TearDown()
//   {
//   }

   /**
    * @brief Dummy placeholder test
    *
    * @param MpiTimerTest Test fixture ID
    * @param Placeholder      Test ID
    */
   TEST_F(MpiTimerTest, Placeholder)
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
