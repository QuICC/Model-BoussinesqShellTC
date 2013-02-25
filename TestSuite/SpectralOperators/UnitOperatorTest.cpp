/** \file UnitOperatorTest.cpp
 *  \brief Implementation of test case for UnitOperator
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the UnitOperator implementation
    */
   class UnitOperatorTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         UnitOperatorTest();

         /**
          * @brief Destructor
          */
         virtual ~UnitOperatorTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   UnitOperatorTest::UnitOperatorTest()
   {
   }

   UnitOperatorTest::~UnitOperatorTest()
   {
   }

//   void UnitOperatorTest::SetUp()
//   {
//   }

//   void UnitOperatorTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(UnitOperatorTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
