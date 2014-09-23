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
#include <sstream>

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Timers/TimerMacro.h"
#include "TypeSelectors/SparseSolverSelector.hpp"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Solve real problem
    */
   Array solveProblem(const SparseMatrix& mat, const Array& sol, const Array& rhs)
   {
      TimerMacro  timer;
      Array timing = Array::Zero(2);

      // Create solver
      Solver::SparseSelector<SparseMatrix>::Type   solver;
   
      // Factorize
      timer.start();
      solver.compute(mat);
      timer.stop();
      timing(0) = timer.time();

      Array lhs(sol.rows());
      timer.start();
      lhs = solver.solve(rhs);
      timer.stop();
      timing(1) = timer.time();

      MHDFloat absError = (sol - lhs).array().abs().maxCoeff();

      return timing;
   }

   /**
    * @brief Solve complex problem
    */
   Array solveProblem(const SparseMatrixZ& mat, const ArrayZ& sol, const ArrayZ& rhs)
   {
      TimerMacro  timer;
      Array timing = Array::Zero(2);

      // Create solver
      Solver::SparseSelector<SparseMatrixZ>::Type   solver;
   
      // Factorize
      solver.compute(mat);

      ArrayZ lhs(sol.rows());
      lhs = solver.solve(rhs);

      MHDFloat absError = (sol - lhs).array().abs().maxCoeff();

      return timing;
   }

   /**
    * @brief Load problem from Matrix Market files
    */
   void loadProblem(SparseMatrix& mat, Array& sol, Array& rhs, const std::string& name)
   {
      rhs.resize(mat.rows());
      sol.resize(mat.cols());
      std::ostringstream   oss;
      oss << mat.rows();
      Eigen::loadMarket(mat, name +"_A_" + oss.str() + ".mtx");
      Eigen::loadMarketVector(sol, name + "_sol_" + oss.str() + ".mtx");
      Eigen::loadMarketVector(rhs, name + "_rhs_" + oss.str() + ".mtx");
   }

   /**
    * @brief Setup and run the real tests
    */
   int runRealTest()
   {  
      // Test counter
      int passedTests = 0;
      int timeLoop = 10;

      // Storage
      SparseMatrix mat;
      Array sol;
      Array rhs;

      std::vector<int> ns1D;
      ns1D.push_back(16); ns1D.push_back(32); ns1D.push_back(64); ns1D.push_back(128); ns1D.push_back(256); ns1D.push_back(512);

      std::vector<int> ns2D;
      ns2D.push_back(8*8); ns2D.push_back(16*16); ns2D.push_back(24*24); ns2D.push_back(32*32); ns2D.push_back(48*48);

      std::vector<int> ns3D;
      ns3D.push_back(8*8*8); ns3D.push_back(10*10*10); ns3D.push_back(12*12*12); ns3D.push_back(14*14*14); ns3D.push_back(16*16*16);
      std::vector<int>::iterator it;

      std::cerr << "#Laplacian1D: " << std::endl;
      for(it = ns1D.begin(); it != ns1D.end(); ++it)
      {
         mat.resize(*it,*it);

         loadProblem(mat, sol, rhs, "laplacian1D");
         Array timing = Array::Zero(2);
         for(int k = 0; k < timeLoop; k++)
         {
            timing = timing + solveProblem(mat, sol, rhs);
         }
         std::cerr << *it << "\t" << timing(0)/timeLoop << "\t" << timing(1)/timeLoop << std::endl;
      }

      std::cerr << "#Bilaplacian1D: " << std::endl;
      for(it = ns1D.begin(); it != ns1D.end(); ++it)
      {
         mat.resize(*it,*it);

         loadProblem(mat, sol, rhs, "bilaplacian1D");
         Array timing = Array::Zero(2);
         for(int k = 0; k < timeLoop; k++)
         {
            timing = timing + solveProblem(mat, sol, rhs);
         }
         std::cerr << *it << "\t" << timing(0)/timeLoop << "\t" << timing(1)/timeLoop << std::endl;
      }

      std::cerr << "#Laplacian2D: " << std::endl;
      for(it = ns2D.begin(); it != ns2D.end(); ++it)
      {
         mat.resize(*it,*it);

         loadProblem(mat, sol, rhs, "laplacian2D");
         Array timing = Array::Zero(2);
         for(int k = 0; k < timeLoop; k++)
         {
            timing = timing + solveProblem(mat, sol, rhs);
         }
         std::cerr << *it << "\t" << timing(0)/timeLoop << "\t" << timing(1)/timeLoop << std::endl;
      }

      std::cerr << "#Bilaplacian2D: " << std::endl;
      for(it = ns2D.begin(); it != ns2D.end(); ++it)
      {
         mat.resize(*it,*it);
         loadProblem(mat, sol, rhs, "bilaplacian2D");
         Array timing = Array::Zero(2);
         for(int k = 0; k < timeLoop; k++)
         {
            timing = timing + solveProblem(mat, sol, rhs);
         }
         std::cerr << *it << "\t" << timing(0)/timeLoop << "\t" << timing(1)/timeLoop << std::endl;
      }

      std::cerr << "#Laplacian3D: " << std::endl;
      for(it = ns3D.begin(); it != ns3D.end(); ++it)
      {
         mat.resize(*it,*it);

         loadProblem(mat, sol, rhs, "laplacian3D");
         Array timing = Array::Zero(2);
         for(int k = 0; k < timeLoop; k++)
         {
            timing = timing + solveProblem(mat, sol, rhs);
         }
         std::cerr << *it << "\t" << timing(0)/timeLoop << "\t" << timing(1)/timeLoop << std::endl;
      }

      std::cerr << "#Bilaplacian3D: " << std::endl;
      for(it = ns3D.begin(); it != ns3D.end(); ++it)
      {
         mat.resize(*it,*it);

         loadProblem(mat, sol, rhs, "bilaplacian3D");
         Array timing = Array::Zero(2);
         for(int k = 0; k < timeLoop; k++)
         {
            timing = timing + solveProblem(mat, sol, rhs);
         }
         std::cerr << *it << "\t" << timing(0)/timeLoop << "\t" << timing(1)/timeLoop << std::endl;
      }

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
