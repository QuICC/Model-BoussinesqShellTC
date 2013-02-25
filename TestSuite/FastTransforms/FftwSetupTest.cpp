/** \file FftwSetupTest.cpp
 *  \brief Implementation of test case for FftwSetup
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the FftwSetup implementation
    */
   class FftwSetupTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         FftwSetupTest();

         /**
          * @brief Destructor
          */
         virtual ~FftwSetupTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   FftwSetupTest::FftwSetupTest()
   {
   }

   FftwSetupTest::~FftwSetupTest()
   {
   }

//   void FftwSetupTest::SetUp()
//   {
//   }

//   void FftwSetupTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(FftwSetupTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
