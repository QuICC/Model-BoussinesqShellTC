/** \file SparseDTimestepper.cpp
 *  \brief Implementation of a general real timestepper structure
 */

// System includes
//

// External includes
//

// Class include
//
#include "Timesteppers/SparseDTimestepper.hpp"

// Project includes
//
#include "Timesteppers/ImExRK3.hpp"

namespace GeoMHDiSCC {

namespace Timestep {

   SparseDTimestepper::SparseDTimestepper(const int start)
      : Solver::SparseDLinearSolver(start)
   {
   }

   SparseDTimestepper::~SparseDTimestepper()
   {
   }

   void SparseDTimestepper::computeRHS(const int step)
   {
      int start = step*this->nSystem();

      Matrix   tmp;

      if(ImExRK3::rhsNN(step) == 0.0)
      {
         for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
         {
            this->mRHSOld.at(i).first = this->mRHSData.at(i).first;
            this->mRHSData.at(i).first = this->mRHSMatrix.at(i+start)*this->mSolution.at(i).first + ImExRK3::rhsN(step)*this->mRHSData.at(i).first;

            this->mRHSOld.at(i).second = this->mRHSData.at(i).second;
            this->mRHSData.at(i).second = this->mRHSMatrix.at(i+start)*this->mSolution.at(i).second + ImExRK3::rhsN(step)*this->mRHSData.at(i).second;
         }
      } else
      {
         for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
         {
            tmp = this->mRHSData.at(i).first;
            this->mRHSData.at(i).first = this->mRHSMatrix.at(i+start)*this->mSolution.at(i).first + ImExRK3::rhsN(step)*this->mRHSData.at(i).first + ImExRK3::rhsNN(step)*this->mRHSOld.at(i).first;
            this->mRHSOld.at(i).first = tmp;

            tmp = this->mRHSData.at(i).second;
            this->mRHSData.at(i).second = this->mRHSMatrix.at(i+start)*this->mSolution.at(i).second + ImExRK3::rhsN(step)*this->mRHSData.at(i).second + ImExRK3::rhsNN(step)*this->mRHSOld.at(i).second;
            this->mRHSOld.at(i).second = tmp;
         }
      }
   }

   void  SparseDTimestepper::updateTimeMatrix(const MHDFloat lhsCoeff, const MHDFloat rhsCoeff, const int step)
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
            SparseMatrix::InnerIterator lhsIt(this->mLHSMatrix.at(start+i),lhsJ);
            SparseMatrix::InnerIterator rhsIt(this->mRHSMatrix.at(start+i),rhsJ);
            for (SparseMatrix::InnerIterator timeIt(this->mTMatrix.at(i),k); timeIt; ++timeIt)
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

   void SparseDTimestepper::initMatrices(const int n)
   {
      // Initialise base matrices
      SparseDLinearSolver::initMatrices(n);

      // Do not reinitialise if work already done by other field
      if(this->mRHSMatrix.size() == 0)
      {
         // Reserve space for the RHS matrices
         this->mRHSMatrix.reserve(n);

         // Initialise storage
         for(int i = 0; i < n; ++i)
         {
            // Create storage for LHS matrices
            this->mRHSMatrix.push_back(SparseMatrix());
         }
      }
   }

   SparseMatrix& SparseDTimestepper::rRHSMatrix(const int idx)
   {
      return this->mRHSMatrix.at(idx);
   }

   SparseMatrix& SparseDTimestepper::rTMatrix(const int idx)
   {
      return this->mTMatrix.at(idx);
   }

   void SparseDTimestepper::addStorage(const int rows, const int cols)
   {
      // Assert for non zero rows and columns
      assert(rows > 0);
      assert(cols > 0);

      SparseDLinearSolver::addStorage(rows,cols);

      // Add storage for old RHS data
      this->mRHSOld.push_back(std::make_pair(Matrix(rows,cols),Matrix(rows,cols)));
      this->mRHSOld.back().first.setZero();
      this->mRHSOld.back().second.setZero();

      // Add storage for the time matrix
      this->mTMatrix.push_back(SparseMatrix());
   }
}
}
