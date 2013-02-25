/** \file FftwToolsTest.cpp
 *  \brief Implementation of test case for FftwTools
 */

#include "FastTransforms/FftwTools.hpp"
#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the FftwTools implementation
    */
   class FftwToolsTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         FftwToolsTest();

         /**
          * @brief Destructor
          */
         virtual ~FftwToolsTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   FftwToolsTest::FftwToolsTest()
   {
   }

   FftwToolsTest::~FftwToolsTest()
   {
   }

//   void FftwToolsTest::SetUp()
//   {
//   }

//   void FftwToolsTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(FftwToolsTest, Constructor)
   {
      int size;
      size = Transform::FftwTools::dealiasFft(10);
      size = Transform::FftwTools::dealiasMixedFft(10);
      size = Transform::FftwTools::optimizeFft(10);

      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
