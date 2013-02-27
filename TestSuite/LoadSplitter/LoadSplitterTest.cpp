/** \file LoadSplitterTest.cpp
 *  \brief Implementation of test case for LoadSplitter
 */

#include "LoadSplitter/LoadSplitter.hpp"
#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the LoadSplitter implementation
    */
   class LoadSplitterTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         LoadSplitterTest();

         /**
          * @brief Destructor
          */
         virtual ~LoadSplitterTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   LoadSplitterTest::LoadSplitterTest()
   {
   }

   LoadSplitterTest::~LoadSplitterTest()
   {
   }

//   void LoadSplitterTest::SetUp()
//   {
//   }

//   void LoadSplitterTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(LoadSplitterTest, Constructor)
   {  
      int id = 0;
      int nCpu = 16;
      Parallel::LoadSplitter   splitter(id, nCpu);

      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
