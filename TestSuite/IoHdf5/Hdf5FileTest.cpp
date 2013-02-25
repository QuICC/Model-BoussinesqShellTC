/** \file Hdf5FileTest.cpp
 *  \brief Implementation of test case for Hdf5File
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the Hdf5File implementation
    */
   class Hdf5FileTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         Hdf5FileTest();

         /**
          * @brief Destructor
          */
         virtual ~Hdf5FileTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   Hdf5FileTest::Hdf5FileTest()
   {
   }

   Hdf5FileTest::~Hdf5FileTest()
   {
   }

//   void Hdf5FileTest::SetUp()
//   {
//   }

//   void Hdf5FileTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(Hdf5FileTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
