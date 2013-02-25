/** \file StoragePairProviderTest.cpp
 *  \brief Implementation of test case for StoragePairProvider
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the StoragePairProvider implementation
    */
   class StoragePairProviderTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         StoragePairProviderTest();

         /**
          * @brief Destructor
          */
         virtual ~StoragePairProviderTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   StoragePairProviderTest::StoragePairProviderTest()
   {
   }

   StoragePairProviderTest::~StoragePairProviderTest()
   {
   }

//   void StoragePairProviderTest::SetUp()
//   {
//   }

//   void StoragePairProviderTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(StoragePairProviderTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
