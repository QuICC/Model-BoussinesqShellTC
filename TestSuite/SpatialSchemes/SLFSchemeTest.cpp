/** 
 * @file SLFSchemeTest.cpp
 * @brief Implementation of test cases for the spherical Chebyshev/Associated Legendre/Fourier spatial scheme
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#include "gtest/gtest.h"

#include "Framework/FrameworkMacro.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the SLFScheme implementation
    */
   class SLFSchemeTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         SLFSchemeTest();

         /**
          * @brief Destructor
          */
         virtual ~SLFSchemeTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   SLFSchemeTest::SLFSchemeTest()
   {
   }

   SLFSchemeTest::~SLFSchemeTest()
   {
   }

//   void SLFSchemeTest::SetUp()
//   {
//   }

//   void SLFSchemeTest::TearDown()
//   {
//   }

   /**
    * @brief Dummy placeholder test
    *
    * @param SLFSchemeTest Test fixture ID
    * @param Placeholder      Test ID
    */
   TEST_F(SLFSchemeTest, Placeholder)
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
