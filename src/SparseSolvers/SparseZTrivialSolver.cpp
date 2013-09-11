/** 
 * @file SparseZTrivialSolver.cpp
 * @brief Implementation of a general complex trivial solver structure
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
#include "SparseSolvers/SparseZTrivialSolver.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"

namespace GeoMHDiSCC {

namespace Solver {

   SparseZTrivialSolver::SparseZTrivialSolver(const int start)
      : SparseSolverBase(start)
   {
   }

   SparseZTrivialSolver::~SparseZTrivialSolver()
   {
   }

   void SparseZTrivialSolver::addStorage(const int rows, const int cols)
   {
      // Assert for non zero rows and columns
      assert(rows > 0);
      assert(cols > 0);

      // Add storage for solution
      this->mSolution.push_back(MatrixZ(rows,cols));
      this->mSolution.back().setZero();
   }

   int SparseZTrivialSolver::nSystem() const
   {
      return this->mSolution.size();
   }

   const MatrixZ& SparseZTrivialSolver::solution(const int idx) const
   {
      return this->mSolution.at(idx);
   }

   MatrixZ& SparseZTrivialSolver::rSolution(const int idx)
   {
      return this->mSolution.at(idx);
   }

   MatrixZ& SparseZTrivialSolver::rRHSData(const int idx)
   {
      // WARNING: this is the same as rSolution. It's used to simplify some implementations
      return this->mSolution.at(idx);
   }
}
}
