/** \file LoadSplitterTest.cpp
 *  \brief Implementation of test case for LoadSplitter
 */

#include "Framework/FrameworkMacro.h"
#include "LoadSplitter/LoadSplitter.hpp"
#include "SpatialSchemes/3D/TFTScheme.hpp"
#include "IoTools/VisualizeResolution.hpp"
#include "IoAscii/DirectAsciiWriter.hpp"
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
    * @brief Test TFTScheme
    */
   TEST_F(LoadSplitterTest, TFTScheme)
   {  
      if(FrameworkMacro::id() == 0)
      {
         Parallel::LoadSplitter   splitter(0, FrameworkMacro::nCpu());

         ArrayI dims(3);
         dims(0) = 64; dims(1) = 64; dims(2) = 64;
         splitter.init<TFTScheme>(dims);

         splitter.showSplittings(10);

         std::string name = "TFT";

         IoAscii::DirectAsciiWriter outFile(name, ".dot", "", "", "");

         outFile.init();

         IoTools::VisualizeResolution::show(outFile.file(), name, splitter.bestSplitting().first);

         outFile.finalize();

         FrameworkMacro::finalize();

         ASSERT_TRUE(true);
      }
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv)
{
   // Initilise framework
   GeoMHDiSCC::FrameworkMacro::init();

   // Set nCpu for serial run
   int nCpu = 1;

   // Set ID and nCpu in MPI case
   #ifdef GEOMHDISCC_MPI
      // Get MPI size
      MPI_Comm_size(MPI_COMM_WORLD, &nCpu);
   #endif //GEOMHDISCC_MPI

   // Setup framework
   GeoMHDiSCC::FrameworkMacro::setup(nCpu);

   ::testing::InitGoogleTest(&argc, argv);
   int status = RUN_ALL_TESTS();

   // Finalise framework
   GeoMHDiSCC::FrameworkMacro::finalize();

   return status;
}
