/** \file VectorFieldTest.cpp
 *  \brief Implementation of test case for VectorField
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the VectorField implementation
    */
   class VectorFieldTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         VectorFieldTest();

         /**
          * @brief Destructor
          */
         virtual ~VectorFieldTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   VectorFieldTest::VectorFieldTest()
   {
   }

   VectorFieldTest::~VectorFieldTest()
   {
   }

//   void VectorFieldTest::SetUp()
//   {
//   }

//   void VectorFieldTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(VectorFieldTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
