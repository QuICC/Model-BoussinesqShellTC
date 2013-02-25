/** \file ScalarField1DTest.cpp
 *  \brief Implementation of test case for ScalarField1D
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the ScalarField1D implementation
    */
   class ScalarField1DTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         ScalarField1DTest();

         /**
          * @brief Destructor
          */
         virtual ~ScalarField1DTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   ScalarField1DTest::ScalarField1DTest()
   {
   }

   ScalarField1DTest::~ScalarField1DTest()
   {
   }

//   void ScalarField1DTest::SetUp()
//   {
//   }

//   void ScalarField1DTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(ScalarField1DTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
