/** \file Communicator2DTest.cpp
 *  \brief Implementation of test case for Communicator2D
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the Communicator2D implementation
    */
   class Communicator2DTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         Communicator2DTest();

         /**
          * @brief Destructor
          */
         virtual ~Communicator2DTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   Communicator2DTest::Communicator2DTest()
   {
   }

   Communicator2DTest::~Communicator2DTest()
   {
   }

//   void Communicator2DTest::SetUp()
//   {
//   }

//   void Communicator2DTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(Communicator2DTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
