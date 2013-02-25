/** \file TransformResolutionTest.cpp
 *  \brief Implementation of test case for TransformResolution
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the TransformResolution implementation
    */
   class TransformResolutionTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         TransformResolutionTest();

         /**
          * @brief Destructor
          */
         virtual ~TransformResolutionTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   TransformResolutionTest::TransformResolutionTest()
   {
   }

   TransformResolutionTest::~TransformResolutionTest()
   {
   }

//   void TransformResolutionTest::SetUp()
//   {
//   }

//   void TransformResolutionTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(TransformResolutionTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
