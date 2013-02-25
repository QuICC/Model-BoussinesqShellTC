/** \file CoreResolutionTest.cpp
 *  \brief Implementation of test case for CoreResolution
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the CoreResolution implementation
    */
   class CoreResolutionTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         CoreResolutionTest();

         /**
          * @brief Destructor
          */
         virtual ~CoreResolutionTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   CoreResolutionTest::CoreResolutionTest()
   {
   }

   CoreResolutionTest::~CoreResolutionTest()
   {
   }

//   void CoreResolutionTest::SetUp()
//   {
//   }

//   void CoreResolutionTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(CoreResolutionTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
