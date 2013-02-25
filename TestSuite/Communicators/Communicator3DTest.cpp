/** \file Communicator3DTest.cpp
 *  \brief Implementation of test case for Communicator3D
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the Communicator3D implementation
    */
   class Communicator3DTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         Communicator3DTest();

         /**
          * @brief Destructor
          */
         virtual ~Communicator3DTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   Communicator3DTest::Communicator3DTest()
   {
   }

   Communicator3DTest::~Communicator3DTest()
   {
   }

//   void Communicator3DTest::SetUp()
//   {
//   }

//   void Communicator3DTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(Communicator3DTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
