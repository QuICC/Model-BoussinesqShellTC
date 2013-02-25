/** \file TimestepperTest.cpp
 *  \brief Implementation of test case for Timestepper
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the Timestepper implementation
    */
   class TimestepperTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         TimestepperTest();

         /**
          * @brief Destructor
          */
         virtual ~TimestepperTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   TimestepperTest::TimestepperTest()
   {
   }

   TimestepperTest::~TimestepperTest()
   {
   }

//   void TimestepperTest::SetUp()
//   {
//   }

//   void TimestepperTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(TimestepperTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
