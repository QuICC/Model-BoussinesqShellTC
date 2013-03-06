/** \file StateFileWriterTest.cpp
 *  \brief Implementation of test case for StateFileWriter
 */

#include "gtest/gtest.h"

#include "Framework/FrameworkMacro.h"
#include "IoVariable/StateFileWriter.hpp"
#include "LoadSplitter/LoadSplitter.hpp"
#include "SpatialSchemes/3D/TFTScheme.hpp"

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
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   StateFileWriterTest::StateFileWriterTest()
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

   StateFileWriterTest::~StateFileWriterTest()
   {
      // Finalise framework
      FrameworkMacro::finalize();
   }

//   void StateFileWriterTest::SetUp()
//   {
//   }

//   void StateFileWriterTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(StateFileWriterTest, TFTState)
   {
      // Set type string
      std::string type = TFTScheme::type();
      bool isRegular = TFTScheme::isRegular();

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
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
