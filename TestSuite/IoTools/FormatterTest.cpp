/** \file FormatterTest.cpp
 *  \brief Implementation of test case for Formatter
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the Formatter implementation
    */
   class FormatterTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         FormatterTest();

         /**
          * @brief Destructor
          */
         virtual ~FormatterTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   FormatterTest::FormatterTest()
   {
   }

   FormatterTest::~FormatterTest()
   {
   }

//   void FormatterTest::SetUp()
//   {
//   }

//   void FormatterTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(FormatterTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
