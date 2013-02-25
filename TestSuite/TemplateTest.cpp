/** \file ++TESTNAME++Test.cpp
 *  \brief Implementation of test case for ++TESTNAME++
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the ++TESTNAME++ implementation
    */
   class ++TESTNAME++Test : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         ++TESTNAME++Test();

         /**
          * @brief Destructor
          */
         virtual ~++TESTNAME++Test();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   ++TESTNAME++Test::++TESTNAME++Test()
   {
   }

   ++TESTNAME++Test::~++TESTNAME++Test()
   {
   }

//   void ++TESTNAME++Test::SetUp()
//   {
//   }

//   void ++TESTNAME++Test::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(++TESTNAME++Test, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
