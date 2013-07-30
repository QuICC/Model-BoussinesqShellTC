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
    * @brief Small and simple "analytical" solution
    */
   TEST_F(SparseSolverTest, Simple)
   {
   }

   /**
    * @brief Cartesian 1D \f$\nabla^2 x = b\f$ with Dirichlet boundary conditions
    */
   TEST_F(SparseSolverTest, HarmonicDirichlet1D)
   {
   }

   /**
    * @brief Cartesian 1D \f$\nabla^2 x = b\f$ with Neumann boundary conditions
    */
   TEST_F(SparseSolverTest, HarmonicNeumann1D)
   {
   }

   /**
    * @brief Cartesian 1D \f$\nabla^4 x = b\f$ with Dirichlet boundary conditions
    */
   TEST_F(SparseSolverTest, BiHarmonicDirichlet1D)
   {
   }

   /**
    * @brief Cartesian 1D \f$\nabla^4 x = b\f$ with Neumann boundary conditions
    */
   TEST_F(SparseSolverTest, BiHarmonicNeumann1D)
   {
   }

   /**
    * @brief Cartesian 2D \f$\nabla_x^2 \otimes\nabla_y^2 x = b\f$ with Dirichlet/Dirichlet boundary conditions
    */
   TEST_F(SparseSolverTest, HarmonicDirichlet2D)
   {
   }

   /**
    * @brief Cartesian 2D \f$\nabla_x^2 \otimes\nabla_y^2 x = b\f$ with Neumann/Neumann boundary conditions
    */
   TEST_F(SparseSolverTest, HarmonicNeumann2D)
   {
   }

   /**
    * @brief Cartesian 2D \f$\nabla_x^2 \otimes\nabla_y^2 x = b\f$ with mixed Dirichlet/Neumann boundary conditions
    */
   TEST_F(SparseSolverTest, HarmonicMixed2D)
   {
   }

   /**
    * @brief Cartesian 2D \f$\nabla_x^4 \otimes\nabla_y^4 x = b\f$ with Dirichlet/Dirichlet boundary conditions
    */
   TEST_F(SparseSolverTest, BiHarmonicDirichlet2D)
   {
   }

   /**
    * @brief Cartesian 2D \f$\nabla_x^4 \otimes\nabla_y^4 x = b\f$ with Neumann/Neumann boundary conditions
    */
   TEST_F(SparseSolverTest, BiHarmonicNeumann2D)
   {
   }

   /**
    * @brief Cartesian 2D \f$\nabla_x^4 \otimes\nabla_y^4 x = b\f$ with mixed Dirichlet/Neumann boundary conditions
    */
   TEST_F(SparseSolverTest, BiHarmonicMixed2D)
   {
   }

   /**
    * @brief Cartesian 2D \f$\nabla_x^4 \otimes\nabla_y^2 x = b\f$ with mixed Dirichlet/Neumann boundary conditions
    */
   TEST_F(SparseSolverTest, MixedHarmonicMixed2D)
   {
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
