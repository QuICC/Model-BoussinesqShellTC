/** \file XmlFileTest.cpp
 *  \brief Implementation of test case for XmlFile
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the XmlFile implementation
    */
   class XmlFileTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         XmlFileTest();

         /**
          * @brief Destructor
          */
         virtual ~XmlFileTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   XmlFileTest::XmlFileTest()
   {
   }

   XmlFileTest::~XmlFileTest()
   {
   }

//   void XmlFileTest::SetUp()
//   {
//   }

//   void XmlFileTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(XmlFileTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
