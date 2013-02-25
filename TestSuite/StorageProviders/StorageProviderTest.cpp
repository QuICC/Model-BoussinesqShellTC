/** \file StorageProviderTest.cpp
 *  \brief Implementation of test case for StorageProvider
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the StorageProvider implementation
    */
   class StorageProviderTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         StorageProviderTest();

         /**
          * @brief Destructor
          */
         virtual ~StorageProviderTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   StorageProviderTest::StorageProviderTest()
   {
   }

   StorageProviderTest::~StorageProviderTest()
   {
   }

//   void StorageProviderTest::SetUp()
//   {
//   }

//   void StorageProviderTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(StorageProviderTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
