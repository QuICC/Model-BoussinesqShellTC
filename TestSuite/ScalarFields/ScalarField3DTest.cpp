/** \file ScalarField3DTest.cpp
 *  \brief Implementation of test case for ScalarField3D
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the ScalarField3D implementation
    */
   class ScalarField3DTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         ScalarField3DTest();

         /**
          * @brief Destructor
          */
         virtual ~ScalarField3DTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   ScalarField3DTest::ScalarField3DTest()
   {
   }

   ScalarField3DTest::~ScalarField3DTest()
   {
   }

//   void ScalarField3DTest::SetUp()
//   {
//   }

//   void ScalarField3DTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(ScalarField3DTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
