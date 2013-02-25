/** \file ExecutionTimerTest.cpp
 *  \brief Implementation of test case for ExecutionTimer
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the ExecutionTimer implementation
    */
   class ExecutionTimerTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         ExecutionTimerTest();

         /**
          * @brief Destructor
          */
         virtual ~ExecutionTimerTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   ExecutionTimerTest::ExecutionTimerTest()
   {
   }

   ExecutionTimerTest::~ExecutionTimerTest()
   {
   }

//   void ExecutionTimerTest::SetUp()
//   {
//   }

//   void ExecutionTimerTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(ExecutionTimerTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
