/** \file StorageProfilerTest.cpp
 *  \brief Implementation of test case for StorageProfiler
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the StorageProfiler implementation
    */
   class StorageProfilerTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         StorageProfilerTest();

         /**
          * @brief Destructor
          */
         virtual ~StorageProfilerTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   StorageProfilerTest::StorageProfilerTest()
   {
   }

   StorageProfilerTest::~StorageProfilerTest()
   {
   }

//   void StorageProfilerTest::SetUp()
//   {
//   }

//   void StorageProfilerTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(StorageProfilerTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
