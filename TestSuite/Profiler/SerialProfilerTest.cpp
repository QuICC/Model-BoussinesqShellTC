/** \file SerialProfilerTest.cpp
 *  \brief Implementation of test case for SerialProfiler
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the SerialProfiler implementation
    */
   class SerialProfilerTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         SerialProfilerTest();

         /**
          * @brief Destructor
          */
         virtual ~SerialProfilerTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   SerialProfilerTest::SerialProfilerTest()
   {
   }

   SerialProfilerTest::~SerialProfilerTest()
   {
   }

//   void SerialProfilerTest::SetUp()
//   {
//   }

//   void SerialProfilerTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(SerialProfilerTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
