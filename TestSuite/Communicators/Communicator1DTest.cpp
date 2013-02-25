/** \file Communicator1DTest.cpp
 *  \brief Implementation of test case for Communicator1D
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the Communicator1D implementation
    */
   class Communicator1DTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         Communicator1DTest();

         /**
          * @brief Destructor
          */
         virtual ~Communicator1DTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   Communicator1DTest::Communicator1DTest()
   {
   }

   Communicator1DTest::~Communicator1DTest()
   {
   }

//   void Communicator1DTest::SetUp()
//   {
//   }

//   void Communicator1DTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(Communicator1DTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
