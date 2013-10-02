/** 
 * @file SparseDLinearSolver.cpp
 * @brief Implementation of a general real linear solver structure
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// External includes
//

// Class include
//
#include "SparseSolvers/SparseDLinearSolver.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Solver {

   SparseDLinearSolver::SparseDLinearSolver(const int start)
      : SparseSolverBase(start)
   {
   }

   SparseDLinearSolver::~SparseDLinearSolver()
   {
   }

   void SparseDLinearSolver::solve(const int step)
   {
      int start = step*this->nSystem();

      // Set unused modes to zero
      for(int i = 0; i < this->mZeroIdx; ++i)
      {
         this->mSolution.at(i).real().setZero();
         this->mSolution.at(i).imag().setZero();
      }

      for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
      {
         // Solve for the real component
         this->mSolution.at(i).real() = this->mSolver.at(i+start)->solve(this->mRHSData.at(i).real());

         // Safety assert for successful solve
         assert(this->mSolver.at(i+start)->info() == Eigen::Success);

         // Solve for the imaginary component
         this->mSolution.at(i).imag() = this->mSolver.at(i+start)->solve(this->mRHSData.at(i).imag());

         // Safety assert for successful solve
         assert(this->mSolver.at(i+start)->info() == Eigen::Success);
      }
   }

   void SparseDLinearSolver::initSolver()
   {
      // Initialise the solver
      this->mSolver.reserve(this->mLHSMatrix.size());
      for(size_t i = 0; i < this->mLHSMatrix.size(); i++)
      {
         SharedPtrMacro<SparseSolverSelector<SparseMatrix>::SolverType >  solver(new SparseSolverSelector<SparseMatrix>::SolverType());

         this->mSolver.push_back(solver);
      }

      // Compute the pattern and the factorisations
      this->updateSolver();
   }

   void SparseDLinearSolver::updateSolver()
   {
      // Compute the factorisations
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

   void SparseDLinearSolver::initMatrices(const int n)
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
            this->mLHSMatrix.push_back(SparseMatrix());
         }
      }
   }

   void SparseDLinearSolver::addStorage(const int rows, const int cols)
   {
      // Assert for non zero rows and columns
      assert(rows > 0);
      assert(cols > 0);

      // Add storage for RHS data
      this->mRHSData.push_back(DecoupledZMatrix(rows,cols));
      this->mRHSData.back().real().setZero();
      this->mRHSData.back().imag().setZero();

      // Add storage for solution
      this->mSolution.push_back(DecoupledZMatrix(rows,cols));
      this->mSolution.back().real().setZero();
      this->mSolution.back().imag().setZero();
   }

   int SparseDLinearSolver::nSystem() const
   {
      return this->mRHSData.size();
   }

   SparseMatrix& SparseDLinearSolver::rLHSMatrix(const int idx)
   {
      return this->mLHSMatrix.at(idx);
   }

   DecoupledZMatrix& SparseDLinearSolver::rRHSData(const int idx)
   {
      return this->mRHSData.at(idx);
   }

   const DecoupledZMatrix& SparseDLinearSolver::solution(const int idx) const
   {
      return this->mSolution.at(idx);
   }

   DecoupledZMatrix& SparseDLinearSolver::rSolution(const int idx)
   {
      return this->mSolution.at(idx);
   }
}
}
