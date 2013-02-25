/** \file Transform2DCoordinatorTest.cpp
 *  \brief Implementation of test case for Transform2DCoordinator
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the Transform2DCoordinator implementation
    */
   class Transform2DCoordinatorTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         Transform2DCoordinatorTest();

         /**
          * @brief Destructor
          */
         virtual ~Transform2DCoordinatorTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   Transform2DCoordinatorTest::Transform2DCoordinatorTest()
   {
   }

   Transform2DCoordinatorTest::~Transform2DCoordinatorTest()
   {
   }

//   void Transform2DCoordinatorTest::SetUp()
//   {
//   }

//   void Transform2DCoordinatorTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(Transform2DCoordinatorTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
