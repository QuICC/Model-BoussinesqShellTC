/** \file ControlFileTest.cpp
 *  \brief Implementation of test case for ControlFile
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the ControlFile implementation
    */
   class ControlFileTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         ControlFileTest();

         /**
          * @brief Destructor
          */
         virtual ~ControlFileTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   ControlFileTest::ControlFileTest()
   {
   }

   ControlFileTest::~ControlFileTest()
   {
   }

//   void ControlFileTest::SetUp()
//   {
//   }

//   void ControlFileTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(ControlFileTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
