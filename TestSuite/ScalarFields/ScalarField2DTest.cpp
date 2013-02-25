/** \file ScalarField2DTest.cpp
 *  \brief Implementation of test case for ScalarField2D
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the ScalarField2D implementation
    */
   class ScalarField2DTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         ScalarField2DTest();

         /**
          * @brief Destructor
          */
         virtual ~ScalarField2DTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   ScalarField2DTest::ScalarField2DTest()
   {
   }

   ScalarField2DTest::~ScalarField2DTest()
   {
   }

//   void ScalarField2DTest::SetUp()
//   {
//   }

//   void ScalarField2DTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(ScalarField2DTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
