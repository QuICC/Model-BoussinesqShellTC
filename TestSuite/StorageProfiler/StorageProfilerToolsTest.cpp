/** \file StorageProfilerToolsTest.cpp
 *  \brief Implementation of test case for StorageProfilerTools
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the StorageProfilerTools implementation
    */
   class StorageProfilerToolsTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         StorageProfilerToolsTest();

         /**
          * @brief Destructor
          */
         virtual ~StorageProfilerToolsTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   StorageProfilerToolsTest::StorageProfilerToolsTest()
   {
   }

   StorageProfilerToolsTest::~StorageProfilerToolsTest()
   {
   }

//   void StorageProfilerToolsTest::SetUp()
//   {
//   }

//   void StorageProfilerToolsTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(StorageProfilerToolsTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
