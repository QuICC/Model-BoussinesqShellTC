/** \file Transform3DCoordinatorTest.cpp
 *  \brief Implementation of test case for Transform3DCoordinator
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the Transform3DCoordinator implementation
    */
   class Transform3DCoordinatorTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         Transform3DCoordinatorTest();

         /**
          * @brief Destructor
          */
         virtual ~Transform3DCoordinatorTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   Transform3DCoordinatorTest::Transform3DCoordinatorTest()
   {
   }

   Transform3DCoordinatorTest::~Transform3DCoordinatorTest()
   {
   }

//   void Transform3DCoordinatorTest::SetUp()
//   {
//   }

//   void Transform3DCoordinatorTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(Transform3DCoordinatorTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
