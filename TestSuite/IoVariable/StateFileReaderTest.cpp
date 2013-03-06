/** \file StateFileReaderTest.cpp
 *  \brief Implementation of test case for StateFileReader
 */

#include "gtest/gtest.h"

#include "Framework/FrameworkMacro.h"
#include "IoVariable/StateFileReader.hpp"
#include "LoadSplitter/LoadSplitter.hpp"
#include "SpatialSchemes/3D/TFTScheme.hpp"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the StateFileReader implementation
    */
   class StateFileReaderTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         StateFileReaderTest();

         /**
          * @brief Destructor
          */
         virtual ~StateFileReaderTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   StateFileReaderTest::StateFileReaderTest()
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

   StateFileReaderTest::~StateFileReaderTest()
   {
      // Finalise framework
      FrameworkMacro::finalize();
   }

//   void StateFileReaderTest::SetUp()
//   {
//   }

//   void StateFileReaderTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(StateFileReaderTest, Constructor)
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
      IoVariable::StateFileReader   state("0000", type, isRegular);

      state.expect(PhysicalNames::TEMPERATURE);
      state.addScalar(std::make_pair(PhysicalNames::TEMPERATURE, spScalar));

      // Make sure all information was provided
      if(state.isFull())
      {
         // Initialise state file
         state.init();

         // Write state file
         state.read();

         // Finalize state file
         state.finalize();
      } else
      {
         ASSERT_TRUE(false);
      }

      // Create test data
      for(int k = 0; k <= dim(1); ++k)
      {
         for(int j = 0; j <= dim(2); ++j)
         {
            for(int i = 0; i <= dim(0); ++i)
            {
               EXPECT_EQ(spScalar->dom(0).perturbation().point(i,j,k), MHDComplex(static_cast<MHDFloat>(-k),static_cast<MHDFloat>(i)*0.01+static_cast<MHDFloat>(j)));
            }
         }
      }
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
