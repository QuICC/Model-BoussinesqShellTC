/** \file SparseSolverTest.cpp
 *  \brief Implementation of test case for the sparse solvers
 */

#include "gtest/gtest.h"

#include "Framework/FrameworkMacro.h"
#include "TypeSelectors/SparseSolverSelector.hpp"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the sparse solvers
    */
   class SparseSolverTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         SparseSolverTest();

         /**
          * @brief Destructor
          */
         virtual ~SparseSolverTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   SparseSolverTest::SparseSolverTest()
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

   SparseSolverTest::~SparseSolverTest()
   {
      // Finalise framework
      FrameworkMacro::finalize();
   }

//   void SparseSolverTest::SetUp()
//   {
//   }

//   void SparseSolverTest::TearDown()
//   {
//   }

   /**
    * @brief Test case 1 for sparse solvers
    */
   TEST_F(SparseSolverTest, SparseSolver1)
   {
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
