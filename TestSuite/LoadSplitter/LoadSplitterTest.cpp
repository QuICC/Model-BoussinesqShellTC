/** 
 * @file LoadSplitterTest.cpp
 * @brief Implementation of test cases for LoadSplitter
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#include "Framework/FrameworkMacro.h"
#include "LoadSplitter/LoadSplitter.hpp"
#include "TypeSelectors/SpatialSchemeSelector.hpp"
#include "IoXml/GxlWriter.hpp"
#include "gtest/gtest.h"

#ifdef GEOMHDISCC_SPATIALSCHEME_TTT
   #define SCHEMENAME "TTT"
#endif //GEOMHDISCC_SPATIALSCHEME_TTT
#ifdef GEOMHDISCC_SPATIALSCHEME_TFT
   #define SCHEMENAME "TFT"
#endif //GEOMHDISCC_SPATIALSCHEME_TFT
#ifdef GEOMHDISCC_SPATIALSCHEME_TFF
   #define SCHEMENAME "TFF"
#endif //GEOMHDISCC_SPATIALSCHEME_TFF
#ifdef GEOMHDISCC_SPATIALSCHEME_FFF
   #define SCHEMENAME "FFF"
#endif //GEOMHDISCC_SPATIALSCHEME_FFF
#ifdef GEOMHDISCC_SPATIALSCHEME_CFT
   #define SCHEMENAME "CFT"
#endif //GEOMHDISCC_SPATIALSCHEME_CFT
#ifdef GEOMHDISCC_SPATIALSCHEME_WFT
   #define SCHEMENAME "WFT"
#endif //GEOMHDISCC_SPATIALSCHEME_WFT
#ifdef GEOMHDISCC_SPATIALSCHEME_SLF
   #define SCHEMENAME "SLF"
#endif //GEOMHDISCC_SPATIALSCHEME_SLF
#ifdef GEOMHDISCC_SPATIALSCHEME_WLF
   #define SCHEMENAME "WLF"
#endif //GEOMHDISCC_SPATIALSCHEME_WLF


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
    * @brief Test spatial scheme load splitting
    */
   TEST_F(LoadSplitterTest, SpatialScheme)
   {  
      // Only one CPU needs to work
      if(FrameworkMacro::id() == 0)
      {
         // Create load splitter
         Parallel::LoadSplitter   splitter(0, FrameworkMacro::nCpu());

         // Create test resolution
         ArrayI dims(3);
         dims(0) = 64; dims(1) = 64; dims(2) = 64;

         // Initialise the load splitter with spatial scheme
         splitter.init<Schemes::SpatialType>(dims);

         // Show the 5 best splittings' description
         splitter.showSplittings(5);

         // Set the scheme name
         std::string name = SCHEMENAME;

         // Create GXL graph format file
         IoXml::GxlWriter gxl(name + "_transpose_graph");

         // Initialise and write graph for resolution
         gxl.init();
         gxl.graphResolution(splitter.bestSplitting().first);
         gxl.write();
         gxl.finalize();

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
