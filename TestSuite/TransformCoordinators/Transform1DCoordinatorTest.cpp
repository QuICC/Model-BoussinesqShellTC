/** \file Transform1DCoordinatorTest.cpp
 *  \brief Implementation of test case for Transform1DCoordinator
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the Transform1DCoordinator implementation
    */
   class Transform1DCoordinatorTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         Transform1DCoordinatorTest();

         /**
          * @brief Destructor
          */
         virtual ~Transform1DCoordinatorTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   enum DummyEnum {TRA1D, TRA2D, TRA3D};

   template <DummyEnum ID> class Nothing;

   struct DummyConf {
      template<DummyEnum  ID1, template<DummyEnum ID2> class STEP> static void  project(int dummy);
   };

   template<> void DummyConf::project<TRA1D, Nothing>(int dummy){};
   template<> void DummyConf::project<TRA2D, Nothing>(int dummy){};
   template<> void DummyConf::project<TRA3D, Nothing>(int dummy){};

   template<> void DummyConf::project<TRA1D, Step>(int dummy)
   {

   }

   Transform1DCoordinatorTest::Transform1DCoordinatorTest()
   {
   }

   Transform1DCoordinatorTest::~Transform1DCoordinatorTest()
   {
   }

//   void Transform1DCoordinatorTest::SetUp()
//   {
//   }

//   void Transform1DCoordinatorTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(Transform1DCoordinatorTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
