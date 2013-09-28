/** 
 * @file VisualizationFileWriterTest.cpp
 * @brief Implementation of test case for VisualizationFileWriter
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#include "gtest/gtest.h"

#include "Framework/FrameworkMacro.h"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "LoadSplitter/LoadSplitter.hpp"
#include "TypeSelectors/SpatialSchemeSelector.hpp"

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
         virtual void SetUp();

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};

         /**
          * @brief Maximum order in 1D
          */
         int mFile1D;

         /**
          * @brief Maximum order in 2D
          */
         int mFile2D;

         /**
          * @brief Maximum order in 3D
          */
         int mFile3D;
         
         /**
          * @brief Shared resolution 
          */
         SharedResolution mspRes;
   };

   VisualizationFileWriterTest::VisualizationFileWriterTest()
      : mFile1D(23), mFile2D(25), mFile3D(27)
   {
   }

   VisualizationFileWriterTest::~VisualizationFileWriterTest()
   {
   }

   void VisualizationFileWriterTest::SetUp()
   {
      // Create load splitter
      Parallel::LoadSplitter   splitter(FrameworkMacro::id(), FrameworkMacro::nCpu());

      // Create test resolution
      ArrayI dims(3);
      dims(0) = this->mFile1D; dims(1) = this->mFile2D; dims(2) = this->mFile3D;

      // Initialise the load splitter with spatial scheme
      splitter.init<Schemes::SpatialType>(dims);

      // Get best splitting
      std::pair<SharedResolution, Parallel::SplittingDescription>  best = splitter.bestSplitting();

      // Get resolution
      this->mspRes = best.first;

      // Set additional options on final resolution object
      Schemes::SpatialType::tuneResolution(this->mspRes);
   }

//   void VisualizationFileWriterTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(VisualizationFileWriterTest, WriteScalar)
   {
      // Create scalar variable
      Datatypes::SharedScalarVariableType   spScalar(new Datatypes::ScalarVariableType(this->mspRes));

      // Initialise physical component
      spScalar->initPhysical();
      spScalar->initPhysicalDiff();

      // Create test data
      for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>(); ++k)
      {
         for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(); ++j)
         {
            for(int i = 0; i < this->mspRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DATF1D>(); ++i)
            {
               spScalar->rDom(0).rPhys().setPoint(static_cast<MHDFloat>(i) + static_cast<MHDFloat>(j)*0.0001 + static_cast<MHDFloat>(k)*0.01,i,j,k);
               spScalar->rDom(0).rGrad().rComp(FieldComponents::Physical::ONE).setPoint(static_cast<MHDFloat>(i) + static_cast<MHDFloat>(j)*0.0001 + static_cast<MHDFloat>(k)*0.01,i,j,k);
               spScalar->rDom(0).rGrad().rComp(FieldComponents::Physical::ONE).setPoint(-(static_cast<MHDFloat>(i) + static_cast<MHDFloat>(j)*0.0001 + static_cast<MHDFloat>(k)*0.01),i,j,k);
               spScalar->rDom(0).rGrad().rComp(FieldComponents::Physical::ONE).setPoint(static_cast<MHDFloat>(i) + static_cast<MHDFloat>(j)*0.0001 + static_cast<MHDFloat>(k)*0.01,i,j,k);
            }
         }
      }

      // Set type string
      std::string type = Schemes::SpatialType::type();

      // Create state file
      IoVariable::VisualizationFileWriter viz(type);

      // Set expected variable
      viz.expect(PhysicalNames::TEMPERATURE);
      viz.addScalar(std::make_pair(PhysicalNames::TEMPERATURE, spScalar));

      // Set fake mesh
      std::vector<Array>   meshs;
      int n;
      for(int i = 0; i < this->mspRes->cpu()->nDim(); ++i)
      {
         n = this->mspRes->sim()->dim(static_cast<Dimensions::Simulation::Id>(i),Dimensions::Space::PHYSICAL);
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
