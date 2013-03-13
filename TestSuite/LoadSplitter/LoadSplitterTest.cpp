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
      // Only one cPU needs to work
      if(FrameworkMacro::id() == 0)
      {
         // Create load splitter
         Parallel::LoadSplitter   splitter(0, FrameworkMacro::nCpu());

         // Create test resolution
         ArrayI dims(3);
         dims(0) = 64; dims(1) = 64; dims(2) = 64;

         // Initialise the load splitter with TFT scheme
         splitter.init<TFTScheme>(dims);

         // Show the 5 best splittings' description
         splitter.showSplittings(5);

         // Set dot file name
         std::string name = "TFT";

         // Create direct access ASCII file
         IoAscii::DirectAsciiWriter outFile(name, ".dot", "Transpose communication graph", "Graphviz dot file", "1.0");

         // Initialise output file
         outFile.init();

         // Write communication graph structure to file
         IoTools::VisualizeResolution::show(outFile.file(), name, splitter.bestSplitting().first);

         // Finalise output file
         outFile.finalize();

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
