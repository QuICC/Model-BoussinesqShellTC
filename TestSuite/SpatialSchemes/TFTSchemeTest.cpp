/** \file TFTSchemeTest.cpp
 *  \brief Implementation of test case for TFTScheme
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the TFTScheme implementation
    */
   class TFTSchemeTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         TFTSchemeTest();

         /**
          * @brief Destructor
          */
         virtual ~TFTSchemeTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   TFTSchemeTest::TFTSchemeTest()
   {
   }

   TFTSchemeTest::~TFTSchemeTest()
   {
   }

//   void TFTSchemeTest::SetUp()
//   {
//   }

//   void TFTSchemeTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(TFTSchemeTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
