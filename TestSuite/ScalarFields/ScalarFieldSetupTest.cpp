/** \file ScalarFieldSetupTest.cpp
 *  \brief Implementation of test case for ScalarFieldSetup
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the ScalarFieldSetup implementation
    */
   class ScalarFieldSetupTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         ScalarFieldSetupTest();

         /**
          * @brief Destructor
          */
         virtual ~ScalarFieldSetupTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   ScalarFieldSetupTest::ScalarFieldSetupTest()
   {
   }

   ScalarFieldSetupTest::~ScalarFieldSetupTest()
   {
   }

//   void ScalarFieldSetupTest::SetUp()
//   {
//   }

//   void ScalarFieldSetupTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(ScalarFieldSetupTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
