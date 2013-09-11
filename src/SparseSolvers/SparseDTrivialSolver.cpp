/** 
 * @file SparseDTrivialSolver.cpp
 * @brief Implementation of a general real trivial solver structure
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

// System includes
//

// External includes
//

// Class include
//
#include "SparseSolvers/SparseDTrivialSolver.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Solver {

   SparseDTrivialSolver::SparseDTrivialSolver(const int start)
      : SparseSolverBase(start)
   {
   }

   SparseDTrivialSolver::~SparseDTrivialSolver()
   {
   }

   void SparseDTrivialSolver::addStorage(const int rows, const int cols)
   {
      // Assert for non zero rows and columns
      assert(rows > 0);
      assert(cols > 0);

      // Add storage for solution
      this->mSolution.push_back(std::make_pair(Matrix(rows,cols),Matrix(rows,cols)));
      this->mSolution.back().first.setZero();
      this->mSolution.back().second.setZero();
   }

   int SparseDTrivialSolver::nSystem() const
   {
      return this->mSolution.size();
   }

   const DecoupledZMatrix& SparseDTrivialSolver::solution(const int idx) const
   {
      return this->mSolution.at(idx);
   }

   DecoupledZMatrix& SparseDTrivialSolver::rSolution(const int idx)
   {
      return this->mSolution.at(idx);
   }

   DecoupledZMatrix& SparseDTrivialSolver::rRHSData(const int idx)
   {
      // WARNING: this is the same as rSolution. It's used to simplify some implementations
      return this->mSolution.at(idx);
   }
}
}
