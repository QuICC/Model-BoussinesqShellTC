/** \file ResolutionTest.cpp
 *  \brief Implementation of test case for Resolution
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the Resolution implementation
    */
   class ResolutionTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         ResolutionTest();

         /**
          * @brief Destructor
          */
         virtual ~ResolutionTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   ResolutionTest::ResolutionTest()
   {
   }

   ResolutionTest::~ResolutionTest()
   {
   }

//   void ResolutionTest::SetUp()
//   {
//   }

//   void ResolutionTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(ResolutionTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
