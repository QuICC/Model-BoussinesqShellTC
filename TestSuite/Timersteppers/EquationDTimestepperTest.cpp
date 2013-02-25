/** \file EquationDTimestepperTest.cpp
 *  \brief Implementation of test case for EquationDTimestepper
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the EquationDTimestepper implementation
    */
   class EquationDTimestepperTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         EquationDTimestepperTest();

         /**
          * @brief Destructor
          */
         virtual ~EquationDTimestepperTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   EquationDTimestepperTest::EquationDTimestepperTest()
   {
   }

   EquationDTimestepperTest::~EquationDTimestepperTest()
   {
   }

//   void EquationDTimestepperTest::SetUp()
//   {
//   }

//   void EquationDTimestepperTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(EquationDTimestepperTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
