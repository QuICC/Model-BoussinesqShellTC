/** \file BinaryFileTest.cpp
 *  \brief Implementation of test case for BinaryFile
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the BinaryFile implementation
    */
   class BinaryFileTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         BinaryFileTest();

         /**
          * @brief Destructor
          */
         virtual ~BinaryFileTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   BinaryFileTest::BinaryFileTest()
   {
   }

   BinaryFileTest::~BinaryFileTest()
   {
   }

//   void BinaryFileTest::SetUp()
//   {
//   }

//   void BinaryFileTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(BinaryFileTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
