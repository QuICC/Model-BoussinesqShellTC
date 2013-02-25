/** \file MpiProfilerTest.cpp
 *  \brief Implementation of test case for MpiProfiler
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the MpiProfiler implementation
    */
   class MpiProfilerTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         MpiProfilerTest();

         /**
          * @brief Destructor
          */
         virtual ~MpiProfilerTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   MpiProfilerTest::MpiProfilerTest()
   {
   }

   MpiProfilerTest::~MpiProfilerTest()
   {
   }

//   void MpiProfilerTest::SetUp()
//   {
//   }

//   void MpiProfilerTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(MpiProfilerTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
