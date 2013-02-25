/** \file ConfigurationWriterTest.cpp
 *  \brief Implementation of test case for ConfigurationWriter
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the ConfigurationWriter implementation
    */
   class ConfigurationWriterTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         ConfigurationWriterTest();

         /**
          * @brief Destructor
          */
         virtual ~ConfigurationWriterTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   ConfigurationWriterTest::ConfigurationWriterTest()
   {
   }

   ConfigurationWriterTest::~ConfigurationWriterTest()
   {
   }

//   void ConfigurationWriterTest::SetUp()
//   {
//   }

//   void ConfigurationWriterTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(ConfigurationWriterTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
