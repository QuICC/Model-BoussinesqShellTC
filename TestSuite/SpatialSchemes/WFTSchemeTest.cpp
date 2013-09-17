/** 
 * @file WFTSchemeTest.cpp
 * @brief Implementation of test cases for the Worland/Fourier/Chebyshev spatial scheme
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#include "gtest/gtest.h"

#include "Framework/FrameworkMacro.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the WFTScheme implementation
    */
   class WFTSchemeTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         WFTSchemeTest();

         /**
          * @brief Destructor
          */
         virtual ~WFTSchemeTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   WFTSchemeTest::WFTSchemeTest()
   {
   }

   WFTSchemeTest::~WFTSchemeTest()
   {
   }

//   void WFTSchemeTest::SetUp()
//   {
//   }

//   void WFTSchemeTest::TearDown()
//   {
//   }

   /**
    * @brief Dummy placeholder test
    *
    * @param WFTSchemeTest Test fixture ID
    * @param Placeholder      Test ID
    */
   TEST_F(WFTSchemeTest, Placeholder)
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
