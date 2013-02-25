/** \file ProfilerToolsTest.cpp
 *  \brief Implementation of test case for ProfilerTools
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the ProfilerTools implementation
    */
   class ProfilerToolsTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         ProfilerToolsTest();

         /**
          * @brief Destructor
          */
         virtual ~ProfilerToolsTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   ProfilerToolsTest::ProfilerToolsTest()
   {
   }

   ProfilerToolsTest::~ProfilerToolsTest()
   {
   }

//   void ProfilerToolsTest::SetUp()
//   {
//   }

//   void ProfilerToolsTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(ProfilerToolsTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
