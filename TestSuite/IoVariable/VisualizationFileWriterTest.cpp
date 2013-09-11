/** 
 * @file VisualizationFileWriterTest.cpp
 * @brief Implementation of test case for VisualizationFileWriter
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

#include "gtest/gtest.h"

#include "Framework/FrameworkMacro.h"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "LoadSplitter/LoadSplitter.hpp"
#include "SpatialSchemes/3D/TFTScheme.hpp"

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
      // Set type string
      std::string type = TFTScheme::type();

      // Set spectral and physical dimensions
      ArrayI dim(3); dim(0) = 11; dim(1) = 13; dim(2) = 12;
     
      // Create the load splitter
      Parallel::LoadSplitter splitter(FrameworkMacro::id(), FrameworkMacro::nCpu());

      // Initialise the load splitter
      splitter.init<TFTScheme>(dim);

      // Get best splitting resolution object
      std::pair<SharedResolution, Parallel::SplittingDescription>  best = splitter.bestSplitting();

      // Store the shared resolution object
      SharedResolution spRes = best.first;

      // Create scalar variable
      Datatypes::SharedScalarVariableType   spScalar(new Datatypes::ScalarVariableType(spRes));

      // Initialise physical component
      spScalar->initPhysical();

      // Create test data
      for(int k = 0; k < spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>(); ++k)
      {
         for(int j = 0; j < spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(); ++j)
         {
            for(int i = 0; i < spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DATF1D>(); ++i)
            {
               spScalar->rDom(0).rPhys().setPoint(static_cast<MHDFloat>(i) + static_cast<MHDFloat>(j)*0.0001 + static_cast<MHDFloat>(k)*0.01,i,j,k);
            }
         }
      }

      // Create state file
      IoVariable::VisualizationFileWriter viz(type);

      // Set expected variable
      viz.expect(PhysicalNames::TEMPERATURE);
      viz.addScalar(std::make_pair(PhysicalNames::TEMPERATURE, spScalar));

      // Set fake mesh
      std::vector<Array>   meshs;
      int n;
      for(int i = 0; i < dim.size(); ++i)
      {
         n = spRes->sim()->dim(static_cast<Dimensions::Simulation::Id>(i),Dimensions::Space::PHYSICAL);
         meshs.push_back(Array(n));
         meshs.back().setLinSpaced(n,0, n - 1);
      }
      viz.setMesh(meshs);

      // Make sure all information was provided
      if(viz.isFull())
      {
         // Initialise state file
         viz.init();

         // Write state file
         viz.write();

         // Finalize state file
         viz.finalize();
      } else
      {
         ASSERT_TRUE(false);
      }

      ASSERT_TRUE(true);

   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
