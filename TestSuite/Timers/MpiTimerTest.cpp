/** \file MpiTimerTest.cpp
 *  \brief Implementation of test case for MpiTimer
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the MpiTimer implementation
    */
   class MpiTimerTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         MpiTimerTest();

         /**
          * @brief Destructor
          */
         virtual ~MpiTimerTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   MpiTimerTest::MpiTimerTest()
   {
   }

   MpiTimerTest::~MpiTimerTest()
   {
   }

//   void MpiTimerTest::SetUp()
//   {
//   }

//   void MpiTimerTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(MpiTimerTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
