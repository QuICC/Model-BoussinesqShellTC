/** \file PrecisionTest.cpp
 *  \brief Implementation of test case for Precision
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the Precision implementation
    */
   class PrecisionTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         PrecisionTest();

         /**
          * @brief Destructor
          */
         virtual ~PrecisionTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   PrecisionTest::PrecisionTest()
   {
   }

   PrecisionTest::~PrecisionTest()
   {
   }

//   void PrecisionTest::SetUp()
//   {
//   }

//   void PrecisionTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(PrecisionTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
