/** 
 * @file StateFileWriterTest.cpp
 * @brief Implementation of test case for StateFileWriter
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#include "gtest/gtest.h"

#include "Framework/FrameworkMacro.h"
#include "IoVariable/StateFileWriter.hpp"
#include "LoadSplitter/LoadSplitter.hpp"
#include "TypeSelectors/SpatialSchemeSelector.hpp"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the StateFileWriter implementation
    */
   class StateFileWriterTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         StateFileWriterTest();

         /**
          * @brief Destructor
          */
         virtual ~StateFileWriterTest();

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

   StateFileWriterTest::StateFileWriterTest()
      : mFile1D(13), mFile2D(15), mFile3D(17)
   {
   }

   StateFileWriterTest::~StateFileWriterTest()
   {
   }

   void StateFileWriterTest::SetUp()
   {
      // Create load splitter
      Parallel::LoadSplitter   splitter(FrameworkMacro::id(), FrameworkMacro::nCpu());

      // Create test resolution
      ArrayI dims(3);
      dims(0) = this->mFile1D; dims(1) = this->mFile2D; dims(2) = this->mFile3D;

      // Initialise the load splitter with spatial scheme
      splitter.init<Schemes::SpatialSelector>(dims);

      // Get best splitting
      std::pair<SharedResolution, Parallel::SplittingDescription>  best = splitter.bestSplitting();

      // Get resolution
      this->mspRes = best.first;

      // Set additional options on final resolution object
      Schemes::SpatialSelector::tuneResolution(this->mspRes);
   }

//   void StateFileWriterTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(StateFileWriterTest, WriteState)
   {
      // Create scalar variable
      Datatypes::SharedScalarVariableType   spScalar(new Datatypes::ScalarVariableType(this->mspRes));

      // Create test data
      int i_, j_, k_;
      for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
      {
         k_ = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
         for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); ++j)
         {
            j_ = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k);
            for(int i = 0; i < this->mspRes->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::TRANSFORM); ++i)
            {
               i_ = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DATB1D>(i,k);
               spScalar->rDom(0).rPerturbation().setPoint(MHDComplex(static_cast<MHDFloat>(-k_),static_cast<MHDFloat>(i_)*0.01+static_cast<MHDFloat>(j_)),i,j,k);
            }
         }
      }

      // Set type string
      std::string type = Schemes::SpatialSelector::type();
      bool isRegular = Schemes::SpatialSelector::isRegular();

      // Create state file
      IoVariable::StateFileWriter   state(type, isRegular);

      state.expect(PhysicalNames::TEMPERATURE);
      state.addScalar(std::make_pair(PhysicalNames::TEMPERATURE, spScalar));

      // Make sure all information was provided
      if(state.isFull())
      {
         // Initialise state file
         state.init();

         // Write state file
         state.write();

         // Finalize state file
         state.finalize();
      } else
      {
         ASSERT_TRUE(false);
      }

      ASSERT_TRUE(true);

      // Synchronize over CPUs
      FrameworkMacro::synchronize();
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
