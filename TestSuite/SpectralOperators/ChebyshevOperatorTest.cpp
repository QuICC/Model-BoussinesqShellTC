/** \file ChebyshevOperatorTest.cpp
 *  \brief Implementation of test case for ChebyshevOperator
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the ChebyshevOperator implementation
    */
   class ChebyshevOperatorTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         ChebyshevOperatorTest();

         /**
          * @brief Destructor
          */
         virtual ~ChebyshevOperatorTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   ChebyshevOperatorTest::ChebyshevOperatorTest()
   {
   }

   ChebyshevOperatorTest::~ChebyshevOperatorTest()
   {
   }

//   void ChebyshevOperatorTest::SetUp()
//   {
//   }

//   void ChebyshevOperatorTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(ChebyshevOperatorTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
