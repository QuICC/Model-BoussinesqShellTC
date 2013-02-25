/** \file VectorVariableTest.cpp
 *  \brief Implementation of test case for VectorVariable
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the VectorVariable implementation
    */
   class VectorVariableTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         VectorVariableTest();

         /**
          * @brief Destructor
          */
         virtual ~VectorVariableTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   VectorVariableTest::VectorVariableTest()
   {
   }

   VectorVariableTest::~VectorVariableTest()
   {
   }

//   void VectorVariableTest::SetUp()
//   {
//   }

//   void VectorVariableTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(VectorVariableTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
