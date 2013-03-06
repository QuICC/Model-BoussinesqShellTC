/** \file VisualizationFileWriterTest.cpp
 *  \brief Implementation of test case for VisualizationFileWriter
 */

#include "gtest/gtest.h"

#include "Framework/FrameworkMacro.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the VisualizationFileWriter implementation
    */
   class VisualizationFileWriterTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         VisualizationFileWriterTest();

         /**
          * @brief Destructor
          */
         virtual ~VisualizationFileWriterTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   VisualizationFileWriterTest::VisualizationFileWriterTest()
   {
      // Initilise framework
      FrameworkMacro::init();

      // Set nCpu for serial run
      int nCpu = 1;

      // Set ID and nCpu in MPI case
      #ifdef GEOMHDISCC_MPI
         // Get MPI size
         int size;
         MPI_Comm_size(MPI_COMM_WORLD, &nCpu);
      #endif //GEOMHDISCC_MPI

      // Setup framework
      FrameworkMacro::setup(nCpu);
   }

   VisualizationFileWriterTest::~VisualizationFileWriterTest()
   {
      // Finalise framework
      FrameworkMacro::finalize();
   }

//   void VisualizationFileWriterTest::SetUp()
//   {
//   }

//   void VisualizationFileWriterTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(VisualizationFileWriterTest, Constructor)
   {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
