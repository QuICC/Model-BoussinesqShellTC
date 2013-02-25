/** \file ConfigurationReaderTest.cpp
 *  \brief Implementation of test case for ConfigurationReader
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the ConfigurationReader implementation
    */
   class ConfigurationReaderTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         ConfigurationReaderTest();

         /**
          * @brief Destructor
          */
         virtual ~ConfigurationReaderTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   ConfigurationReaderTest::ConfigurationReaderTest()
   {
   }

   ConfigurationReaderTest::~ConfigurationReaderTest()
   {
   }

//   void ConfigurationReaderTest::SetUp()
//   {
//   }

//   void ConfigurationReaderTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(ConfigurationReaderTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
