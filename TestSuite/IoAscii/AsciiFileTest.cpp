/** \file AsciiFileTest.cpp
 *  \brief Implementation of test case for AsciiFile
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the AsciiFile implementation
    */
   class AsciiFileTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         AsciiFileTest();

         /**
          * @brief Destructor
          */
         virtual ~AsciiFileTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   AsciiFileTest::AsciiFileTest()
   {
   }

   AsciiFileTest::~AsciiFileTest()
   {
   }

//   void AsciiFileTest::SetUp()
//   {
//   }

//   void AsciiFileTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(AsciiFileTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
