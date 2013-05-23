/** \file EquationZTimestepper.cpp
 *  \brief Implementation of a general timestepper structure
 */

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Class include
//
#include "Timesteppers/EquationZTimestepper.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"
#include "Timesteppers/ImExRK3.hpp"

namespace GeoMHDiSCC {

namespace Timestep {

   EquationZTimestepper::EquationZTimestepper(const int start)
      : EquationTimestepperBase(start)
   {
   }

   EquationZTimestepper::~EquationZTimestepper()
   {
   }

   void EquationZTimestepper::computeRHS(const int step)
   {
      int start = step*this->nSystem();

      MatrixZ   tmp;

      if(ImExRK3::rhsNN(step) == 0.0)
      {
         for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
         {
            this->mRHSOld.at(i) = this->mRHSData.at(i);
            this->mRHSData.at(i) = this->mRHSMatrix.at(i+start)*this->mSolution.at(i) + ImExRK3::rhsN(step)*this->mRHSData.at(i);
         }
      } else
      {
         for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
         {
            tmp = this->mRHSData.at(i);
            this->mRHSData.at(i) = this->mRHSMatrix.at(i+start)*this->mSolution.at(i) + ImExRK3::rhsN(step)*this->mRHSData.at(i) + ImExRK3::rhsNN(step)*this->mRHSOld.at(i);
            this->mRHSOld.at(i) = tmp;
         }
      }
   }

   void EquationZTimestepper::solve(const int step)
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

   void  EquationZTimestepper::updateTimeMatrix(const MHDFloat lhsCoeff, const MHDFloat rhsCoeff, const int step)
   {
      // Set start offset
      int start = step*this->nSystem();

      // Loop over matrices within same step
      for(int i = 0; i < this->nSystem(); ++i)
      {
         // Get the number of nonzero elements in time dependence
         size_t nnz = this->mTMatrix.at(i).nonZeros();

         // Update LHS and RHS matrices
         size_t lhsJ = 0;
         size_t rhsJ = 0;
         for (size_t k=0; k< static_cast<size_t>(this->mTMatrix.at(i).outerSize()); ++k)
         {
            SparseMatrixZ::InnerIterator lhsIt(this->mLHSMatrix.at(start+i),lhsJ);
            SparseMatrixZ::InnerIterator rhsIt(this->mRHSMatrix.at(start+i),rhsJ);
            for (SparseMatrixZ::InnerIterator timeIt(this->mTMatrix.at(i),k); timeIt; ++timeIt)
            {
               // Only keep going if nonzero elements are left
               if(nnz > 0)
               {
                  // Update LHS matrix
                  if(timeIt.col() == lhsIt.col())
                  {
                     if(timeIt.row() == lhsIt.row())
                     {
                        // Update values
                        lhsIt.valueRef() += lhsCoeff*timeIt.value();

                        // Update LHS iterators and counters
                        ++lhsIt;
                        if(!lhsIt)
                        {
                           lhsJ++;
                        }
                     }
                  }

                  // Update LHS matrix
                  if(timeIt.col() == rhsIt.col())
                  {
                     if(timeIt.row() == rhsIt.row())
                     {
                        // Update values
                        rhsIt.valueRef() += rhsCoeff*timeIt.value();

                        // Update RHS iterators and counters
                        ++rhsIt;
                        if(!rhsIt)
                        {
                           rhsJ++;
                        }
                     }
                  }

                  // Update nonzero counter
                  nnz--;
               } else
               {
                  break;
               }
            }
         }

         // Safety assert to make sure all values have been updated
         assert(nnz == 0);
      }
   }

   void EquationZTimestepper::initSolver()
   {
      // Initialise solver
      this->mSolver.reserve(this->mLHSMatrix.size());
      for(size_t i = 0; i < this->mLHSMatrix.size(); i++)
      {
         SharedPtrMacro<SparseSolverMacro<SparseMatrixZ> >  solver(new SparseSolverMacro<SparseMatrixZ>());

         this->mSolver.push_back(solver);
      }

      // Compute pattern and factorisation
      this->updateSolver();
   }

   void EquationZTimestepper::updateSolver()
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

   void EquationZTimestepper::initMatrices(const int n)
   {
      // Reserve space for the LHS matrices
      this->mLHSMatrix.reserve(n);

      // Reserve space for the RHS matrices
      this->mRHSMatrix.reserve(n);

      // Initialise storage
      for(int i = 0; i < n; ++i)
      {
         // Create storage for LHS matrices
         this->mLHSMatrix.push_back(SparseMatrixZ());

         // Create storage for LHS matrices
         this->mRHSMatrix.push_back(SparseMatrixZ());
      }
   }

   SparseMatrixZ& EquationZTimestepper::rLHSMatrix(const int idx)
   {
      return this->mLHSMatrix.at(idx);
   }

   SparseMatrixZ& EquationZTimestepper::rRHSMatrix(const int idx)
   {
      return this->mRHSMatrix.at(idx);
   }

   SparseMatrixZ& EquationZTimestepper::rTMatrix(const int idx)
   {
      return this->mTMatrix.at(idx);
   }

   void EquationZTimestepper::addStorage(const int rows, const int cols)
   {
      // Assert for non zero rows and columns
      assert(rows > 0);
      assert(cols > 0);

      // Add storage for RHS data
      this->mRHSData.push_back(MatrixZ(rows,cols));
      this->mRHSData.back().setZero();

      // Add storage for old RHS
      this->mRHSOld.push_back(MatrixZ(rows,cols));
      this->mRHSOld.back().setZero();

      // Add storage for solution
      this->mSolution.push_back(MatrixZ(rows,cols));
      this->mSolution.back().setZero();

      // Add storage for the time matrix
      this->mTMatrix.push_back(SparseMatrixZ());
   }
}
}
