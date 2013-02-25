/** \file DirectBinaryWriterTest.cpp
 *  \brief Implementation of test case for DirectBinaryWriter
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the DirectBinaryWriter implementation
    */
   class DirectBinaryWriterTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         DirectBinaryWriterTest();

         /**
          * @brief Destructor
          */
         virtual ~DirectBinaryWriterTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   DirectBinaryWriterTest::DirectBinaryWriterTest()
   {
   }

   DirectBinaryWriterTest::~DirectBinaryWriterTest()
   {
   }

//   void DirectBinaryWriterTest::SetUp()
//   {
//   }

//   void DirectBinaryWriterTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(DirectBinaryWriterTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
