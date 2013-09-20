/** 
 * @file LoadSplitterTest.cpp
 * @brief Implementation of test cases for LoadSplitter
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#include "gtest/gtest.h"

#include "Framework/FrameworkMacro.h"
#include "LoadSplitter/LoadSplitter.hpp"
#include "TypeSelectors/SpatialSchemeSelector.hpp"
#include "IoXml/GxlWriter.hpp"

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

         // Create GXL graph format file
         IoXml::GxlWriter gxl(Schemes::SpatialType::type() + "_transpose_graph");

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

      if(argc == 2)
      {
         if(nCpu == 1)
         {
            std::istringstream   iss;
            iss.str(argv[1]);
            iss >> nCpu;
            std::cerr << nCpu << std::endl;
         } else
         {
            throw GeoMHDiSCC::Exception("Can only emulate large run with single CPU");
         }
      }
   #endif //GEOMHDISCC_MPI

   // Setup framework
   try
   {
      GeoMHDiSCC::FrameworkMacro::setup(nCpu);
   }
   catch(GeoMHDiSCC::Exception& e)
   {
   }

   ::testing::InitGoogleTest(&argc, argv);
   int status = RUN_ALL_TESTS();

   // Finalise framework
   GeoMHDiSCC::FrameworkMacro::finalize();

   return status;
}
