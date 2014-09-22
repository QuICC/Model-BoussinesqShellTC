/**
 * @file SparseLinearSolverTest.cpp
 * @brief Simple tests for the different linear solvers
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//
#include "Framework/FrameworkMacro.h"

// System includes
//
#include <Eigen/SparseExtra>
#include <iostream>

// Project includes
//
#include "Base/Typedefs.hpp"
#include "TypeSelectors/SparseSolverSelector.hpp"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Solve real problem
    */
   int solveProblem(const SparseMatrix& mat, const Array& sol, const Array& rhs)
   {
      // Create solver
      Solver::SparseSelector<SparseMatrix>::Type   solver;
   
      // Factorize
      solver.compute(mat);

      Array lhs(sol.rows());
      lhs = solver.solve(rhs);

      std::cout << "Max absolute error: " << (sol - lhs).array().abs().maxCoeff() << std::endl;

      return 1;
   }

   /**
    * @brief Solve complex problem
    */
   int solveProblem(const SparseMatrixZ& mat, const ArrayZ& sol, const ArrayZ& rhs)
   {
      // Create solver
      Solver::SparseSelector<SparseMatrixZ>::Type   solver;
   
      // Factorize
      solver.compute(mat);

      ArrayZ lhs(sol.rows());
      lhs = solver.solve(rhs);

      MHDFloat absError = (sol - lhs).array().abs().maxCoeff();

      std::cout << "Max absolute error: " << absError << std::endl;

      return 1;
   }

   /**
    * @brief Setup problem 1
    * 
    * Simple real problem take from Pardiso example
    */
   void setupProblem1(SparseMatrix& mat, Array& sol, Array& rhs)
   {
      // Resize provided storage
      mat.resize(8,8);
      sol.resize(8);
      rhs.resize(8);

      // Build test matrix from triplets
      typedef Eigen::Triplet<MHDFloat> T;
      std::vector<T> tripletList;
      tripletList.reserve(20);
      tripletList.push_back(T(0,0,7.0));
      tripletList.push_back(T(0,2,1.0));
      tripletList.push_back(T(0,5,2.0));
      tripletList.push_back(T(0,6,7.0));
      tripletList.push_back(T(1,1,-4.0));
      tripletList.push_back(T(1,2,8.0));
      tripletList.push_back(T(1,4,2.0));
      tripletList.push_back(T(2,2,1.0));
      tripletList.push_back(T(2,7,5.0));
      tripletList.push_back(T(3,3,7.0));
      tripletList.push_back(T(3,6,9.0));
      tripletList.push_back(T(4,1,-4.0));
      tripletList.push_back(T(5,2,7.0));
      tripletList.push_back(T(5,5,3.0));
      tripletList.push_back(T(5,7,8.0));
      tripletList.push_back(T(6,1,1.0));
      tripletList.push_back(T(6,6,11.0));
      tripletList.push_back(T(7,2,-3.0));
      tripletList.push_back(T(7,6,2.0));
      tripletList.push_back(T(7,7,5.0));
      mat.setFromTriplets(tripletList.begin(), tripletList.end());

      // Build RHS
      rhs.setLinSpaced(8,0.,7.);

      // Build exact solution
      sol(0) = -1.153896103896104;sol(1) = -1;sol(2) = -0.9318181818181819;sol(3) = -0.3896103896103896;sol(4) = 2.227272727272728;sol(5) = 2.277272727272727;sol(6) = 0.6363636363636364;sol(7) = 0.5863636363636366;
   }

   /**
    * @brief Load problem from Matrix Market files
    */
   void loadProblem(SparseMatrix& mat, Array& sol, Array& rhs, const std::string& name)
   {
      mat.resize(0,0);
      rhs.resize(0);
      sol.resize(0);
      Eigen::loadMarket(mat, name+"_A.mtx");
      Eigen::loadMarketVector(sol, name + "_sol.mtx");
      Eigen::loadMarketVector(rhs, name + "_rhs.mtx");
   }

   /**
    * @brief Setup and run the real tests
    */
   int runRealTest()
   {  
      // Test counter
      int passedTests = 0;

      // Storage
      SparseMatrix mat;
      Array sol;
      Array rhs;

      setupProblem1(mat, sol, rhs);
      passedTests = passedTests + solveProblem(mat, sol, rhs);

      loadProblem(mat, sol, rhs, "laplacian1D");
      passedTests = passedTests + solveProblem(mat, sol, rhs);

      loadProblem(mat, sol, rhs, "bilaplacian1D");
      passedTests = passedTests + solveProblem(mat, sol, rhs);

      loadProblem(mat, sol, rhs, "laplacian2D");
      passedTests = passedTests + solveProblem(mat, sol, rhs);

      loadProblem(mat, sol, rhs, "bilaplacian2D");
      passedTests = passedTests + solveProblem(mat, sol, rhs);

      loadProblem(mat, sol, rhs, "bilaplacian3D");
      passedTests = passedTests + solveProblem(mat, sol, rhs);

      return passedTests;
   }

   /**
    * @brief Setup and run the complex tests
    */
   int runCplxTest()
   {  
      // Test counter
      int passedTests = 0;

      // Storage
      SparseMatrixZ mat;
      MatrixZ sol;
      MatrixZ rhs;

      return passedTests;
   }
}
}

/**
 * @brief General main, setting up MPI if required
 *
 * The actual program is in run to make sure MPI initialisations
 * are called before anything else and finalisation after destruction
 */
int main(int argc, char* argv[])
{
   // Initilise everything that can't be done inside a class
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

   // Compute simulation
   int realPassed = GeoMHDiSCC::TestSuite::runRealTest();
   int cplxPassed = GeoMHDiSCC::TestSuite::runCplxTest();

   std::cout << "Real tests passed: " << realPassed << std::endl;
   std::cout << "Complex tests passed: " << cplxPassed << std::endl;

   // Finalise everything that can't be done inside a class
   GeoMHDiSCC::FrameworkMacro::finalize();

   return realPassed+cplxPassed;
}
