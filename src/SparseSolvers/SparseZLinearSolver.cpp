/** 
 * @file SparseZLinearSolver.cpp
 * @brief Implementation of a general complex linear solver structure
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
#include "SparseSolvers/SparseZLinearSolver.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"

namespace GeoMHDiSCC {

namespace Solver {

   SparseZLinearSolver::SparseZLinearSolver(const int start)
      : SparseSolverBase(start)
   {
   }

   SparseZLinearSolver::~SparseZLinearSolver()
   {
   }

   void SparseZLinearSolver::solve(const int step)
   {
      int start = step*this->nSystem();

      // Set unused modes to zero
      for(int i = 0; i < this->mZeroIdx; ++i)
      {
         this->mSolution.at(i).setZero();
      }

      // Solve other modes
      for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
      {
         this->mSolution.at(i) = this->mSolver.at(i+start)->solve(this->mRHSData.at(i));

         // Safety assert for successful solve
         assert(this->mSolver.at(i+start)->info() == Eigen::Success);
      }
   }

   void SparseZLinearSolver::initSolver()
   {
      // Initialise solver
      this->mSolver.reserve(this->mLHSMatrix.size());
      for(size_t i = 0; i < this->mLHSMatrix.size(); i++)
      {
         SharedPtrMacro<SparseSolverSelector<SparseMatrixZ>::SolverType >  solver(new SparseSolverSelector<SparseMatrixZ>::SolverType());

         this->mSolver.push_back(solver);
      }

      // Compute pattern and factorisation
      this->updateSolver();
   }

   void SparseZLinearSolver::updateSolver()
   {
      // Compute factorisation
      for(size_t i = 0; i < this->mLHSMatrix.size(); i++)
      {
         if(static_cast<int>(i) % this->nSystem() >= this->mZeroIdx)
         {
            // Safety assert to make sur matrix is compressed
            assert(this->mLHSMatrix.at(i).isCompressed());

            this->mSolver.at(i)->compute(this->mLHSMatrix.at(i));

            // Safety assert for successful factorisation
            assert(this->mSolver.at(i)->info() == Eigen::Success);
         }
      }
   }

   void SparseZLinearSolver::initMatrices(const int n)
   {
      // Do not reinitialise if work already done by other field
      if(this->mLHSMatrix.size() == 0)
      {
         // Reserve space for the LHS matrices
         this->mLHSMatrix.reserve(n);

         // Initialise storage
         for(int i = 0; i < n; ++i)
         {
            // Create storage for LHS matrices
            this->mLHSMatrix.push_back(SparseMatrixZ());
         }
      }
   }

   void SparseZLinearSolver::addStorage(const int rows, const int cols)
   {
      // Assert for non zero rows and columns
      assert(rows > 0);
      assert(cols > 0);

      // Add storage for RHS data
      this->mRHSData.push_back(MatrixZ(rows,cols));
      this->mRHSData.back().setZero();

      // Add storage for solution
      this->mSolution.push_back(MatrixZ(rows,cols));
      this->mSolution.back().setZero();
   }

   int SparseZLinearSolver::nSystem() const
   {
      return this->mRHSData.size();
   }

   SparseMatrixZ& SparseZLinearSolver::rLHSMatrix(const int idx)
   {
      return this->mLHSMatrix.at(idx);
   }

   MatrixZ& SparseZLinearSolver::rRHSData(const int idx)
   {
      return this->mRHSData.at(idx);
   }

   const MatrixZ& SparseZLinearSolver::solution(const int idx) const
   {
      return this->mSolution.at(idx);
   }

   MatrixZ& SparseZLinearSolver::rSolution(const int idx)
   {
      return this->mSolution.at(idx);
   }
}
}
