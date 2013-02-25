/** \file EquationZTimestepperTest.cpp
 *  \brief Implementation of test case for EquationZTimestepper
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the EquationZTimestepper implementation
    */
   class EquationZTimestepperTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         EquationZTimestepperTest();

         /**
          * @brief Destructor
          */
         virtual ~EquationZTimestepperTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   EquationZTimestepperTest::EquationZTimestepperTest()
   {
   }

   EquationZTimestepperTest::~EquationZTimestepperTest()
   {
   }

//   void EquationZTimestepperTest::SetUp()
//   {
//   }

//   void EquationZTimestepperTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(EquationZTimestepperTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
