/** \file EquationDTimestepper.cpp
 *  \brief Implementation of a general timestepper structure
 */

// System includes
//

// External includes
//

// Class include
//
#include "Timesteppers/EquationDTimestepper.hpp"

// Project includes
//
#include "Timesteppers/ImExRK3.hpp"

namespace GeoMHDiSCC {

namespace Timestep {

   EquationDTimestepper::EquationDTimestepper(const int start)
      : SparseDLinearSolver(start)
   {
   }

   EquationDTimestepper::~EquationDTimestepper()
   {
   }

   void EquationDTimestepper::computeRHS(const int step)
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

   void  EquationDTimestepper::updateTimeMatrix(const MHDFloat lhsCoeff, const MHDFloat rhsCoeff, const int step)
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

   void EquationDTimestepper::initMatrices(const int n)
   {
      // Reserve space for the LHS matrices
      this->mLHSMatrix.reserve(n);

      // Reserve space for the RHS matrices
      this->mRHSMatrix.reserve(n);

      // Initialise storage
      for(int i = 0; i < n; ++i)
      {
         // Create storage for LHS matrices
         this->mLHSMatrix.push_back(SparseMatrix());

         // Create storage for LHS matrices
         this->mRHSMatrix.push_back(SparseMatrix());
      }
   }

   SparseMatrix& EquationDTimestepper::rRHSMatrix(const int idx)
   {
      return this->mRHSMatrix.at(idx);
   }

   SparseMatrix& EquationDTimestepper::rTMatrix(const int idx)
   {
      return this->mTMatrix.at(idx);
   }

   void EquationDTimestepper::addStorage(const int rows, const int cols)
   {
      // Assert for non zero rows and columns
      assert(rows > 0);
      assert(cols > 0);

      // Add storage for RHS data
      this->mRHSData.push_back(std::make_pair(Matrix(rows,cols),Matrix(rows,cols)));
      this->mRHSData.back().first.setZero();
      this->mRHSData.back().second.setZero();

      // Add storage for old RHS data
      this->mRHSOld.push_back(std::make_pair(Matrix(rows,cols),Matrix(rows,cols)));
      this->mRHSOld.back().first.setZero();
      this->mRHSOld.back().second.setZero();

      // Add storage for solution
      this->mSolution.push_back(std::make_pair(Matrix(rows,cols),Matrix(rows,cols)));
      this->mSolution.back().first.setZero();
      this->mSolution.back().second.setZero();

      // Add storage for the time matrix
      this->mTMatrix.push_back(SparseMatrix());
   }
}
}
