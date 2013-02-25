/** \file FlatLayoutTest.cpp
 *  \brief Implementation of test case for FlatLayout
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the FlatLayout implementation
    */
   class FlatLayoutTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         FlatLayoutTest();

         /**
          * @brief Destructor
          */
         virtual ~FlatLayoutTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   FlatLayoutTest::FlatLayoutTest()
   {
   }

   FlatLayoutTest::~FlatLayoutTest()
   {
   }

//   void FlatLayoutTest::SetUp()
//   {
//   }

//   void FlatLayoutTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(FlatLayoutTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
