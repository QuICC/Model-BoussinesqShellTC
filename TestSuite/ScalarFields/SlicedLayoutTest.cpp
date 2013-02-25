/** \file SlicedLayoutTest.cpp
 *  \brief Implementation of test case for SlicedLayout
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the SlicedLayout implementation
    */
   class SlicedLayoutTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         SlicedLayoutTest();

         /**
          * @brief Destructor
          */
         virtual ~SlicedLayoutTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   SlicedLayoutTest::SlicedLayoutTest()
   {
   }

   SlicedLayoutTest::~SlicedLayoutTest()
   {
   }

//   void SlicedLayoutTest::SetUp()
//   {
//   }

//   void SlicedLayoutTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(SlicedLayoutTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
