/** \file FftSetupTest.cpp
 *  \brief Implementation of test case for FftSetup
 */

#include "FastTransforms/FftSetup.hpp"
#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the FftSetup implementation
    */
   class FftSetupTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         FftSetupTest();

         /**
          * @brief Destructor
          */
         virtual ~FftSetupTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   FftSetupTest::FftSetupTest()
   {
   }

   FftSetupTest::~FftSetupTest()
   {
   }

//   void FftSetupTest::SetUp()
//   {
//   }

//   void FftSetupTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(FftSetupTest, Constructor) {
      Transform::FftSetup setup(11,11,11,true);

      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
