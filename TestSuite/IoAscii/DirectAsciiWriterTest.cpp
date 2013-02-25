/** \file DirectAsciiWriterTest.cpp
 *  \brief Implementation of test case for DirectAsciiWriter
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the DirectAsciiWriter implementation
    */
   class DirectAsciiWriterTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         DirectAsciiWriterTest();

         /**
          * @brief Destructor
          */
         virtual ~DirectAsciiWriterTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   DirectAsciiWriterTest::DirectAsciiWriterTest()
   {
   }

   DirectAsciiWriterTest::~DirectAsciiWriterTest()
   {
   }

//   void DirectAsciiWriterTest::SetUp()
//   {
//   }

//   void DirectAsciiWriterTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(DirectAsciiWriterTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
