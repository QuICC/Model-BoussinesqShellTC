/** \file PeriodicOperatorTest.cpp
 *  \brief Implementation of test case for PeriodicOperator
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the PeriodicOperator implementation
    */
   class PeriodicOperatorTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         PeriodicOperatorTest();

         /**
          * @brief Destructor
          */
         virtual ~PeriodicOperatorTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   PeriodicOperatorTest::PeriodicOperatorTest()
   {
   }

   PeriodicOperatorTest::~PeriodicOperatorTest()
   {
   }

//   void PeriodicOperatorTest::SetUp()
//   {
//   }

//   void PeriodicOperatorTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(PeriodicOperatorTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
