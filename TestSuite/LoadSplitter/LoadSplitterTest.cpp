/** \file LoadSplitterTest.cpp
 *  \brief Implementation of test case for LoadSplitter
 */

#include "Framework/FrameworkMacro.h"
#include "LoadSplitter/LoadSplitter.hpp"
#include "SpatialSchemes/3D/TFTScheme.hpp"
#include "IoTools/VisualizeResolution.hpp"
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
      FrameworkMacro::init();

      int id = 0;
      int nCpu = 1;

      FrameworkMacro::setup(nCpu);

      Parallel::LoadSplitter   splitter(id, nCpu);

      ArrayI dims(3);
      dims(0) = 10; dims(1) = 12; dims(2) = 15;
      splitter.init<TFTScheme>(dims);

      splitter.showSplittings(10);

      IoTools::VisualizeResolution::show(std::cout, splitter.bestSplitting().first);

      FrameworkMacro::finalize();

      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
