/** \file StdOutPipeTest.cpp
 *  \brief Implementation of test case for StdOutPipe
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the StdOutPipe implementation
    */
   class StdOutPipeTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         StdOutPipeTest();

         /**
          * @brief Destructor
          */
         virtual ~StdOutPipeTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   StdOutPipeTest::StdOutPipeTest()
   {
   }

   StdOutPipeTest::~StdOutPipeTest()
   {
   }

//   void StdOutPipeTest::SetUp()
//   {
//   }

//   void StdOutPipeTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(StdOutPipeTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
