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

   EquationZTimestepper::EquationZTimestepper(const int nField)
      : EquationTimestepperBase(nField)
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
         for(size_t i = 0; i < this->mRHSData.size(); i++)
         {
            this->mRHSOld.at(i) = this->mRHSData.at(i);
            this->mRHSData.at(i) = this->mRHSMatrix.at(i+start)*this->mSolution.at(i) + ImExRK3::rhsN(step)*this->mRHSData.at(i);
         }
      } else
      {
         for(size_t i = 0; i < this->mRHSData.size(); i++)
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

      // Set m = 0 mode to zero
      this->mSolution.at(0).setZero();

      for(size_t i = 1; i < this->mRHSData.size(); i++)
      {
         this->mSolution.at(i) = this->mSolver.at(i+start)->solve(this->mRHSData.at(i));

         // Safety assert for successful solve
         assert(this->mSolver.at(i+start)->info() == Eigen::Success);
      }
   }

   void EquationZTimestepper::initSolver()
   {
      this->mSolver.reserve(this->mLHSMatrix.size());
      for(size_t i = 0; i < this->mLHSMatrix.size(); i++)
      {
         SharedPtrMacro<SparseSolverMacro<SparseMatrixZ> >  solver(new SparseSolverMacro<SparseMatrixZ>());

         this->mSolver.push_back(solver);
      }

      // Compute pattern and factorisation
      for(size_t i = 0; i < this->mLHSMatrix.size(); i++)
      {
         if(i % this->nSystem() != 0)
         {
            // Safety assert to make sur matrix is compressed
            assert(this->mLHSMatrix.at(i).isCompressed());

            this->mSolver.at(i)->compute(this->mLHSMatrix.at(i));

            // Safety assert for successful factorisation
            assert(this->mSolver.at(i)->info() == Eigen::Success);
         }
      }
   }

   void EquationZTimestepper::updateSolver()
   {
      // Compute factorisation
      for(size_t i = 0; i < this->mLHSMatrix.size(); i++)
      {
         if(i % this->nSystem() != 0)
         {
            // Safety assert to make sur matrix is compressed
            assert(this->mLHSMatrix.at(i).isCompressed());

            this->mSolver.at(i)->factorize(this->mLHSMatrix.at(i));

            // Safety assert for successful factorisation
            assert(this->mSolver.at(i)->info() == Eigen::Success);
         }
      }
   }

   void EquationZTimestepper::reserveMatrices(const int n)
   {
      // Reserve space for the LHS matrices
      this->mLHSMatrix.reserve(n);

      // Reserve space for the RHS matrices
      this->mRHSMatrix.reserve(n);
   }

   void EquationZTimestepper::addLHSMatrix(const DecoupledZSparse& lhs)
   {
      this->mLHSMatrix.push_back(lhs.first.cast<MHDComplex>() + MathConstants::cI*lhs.second);
   }

   void EquationZTimestepper::addRHSMatrix(const DecoupledZSparse& rhs)
   {
      this->mRHSMatrix.push_back(rhs.first.cast<MHDComplex>() + MathConstants::cI*rhs.second);
   }

   void EquationZTimestepper::setLHSMatrix(const int idx, const DecoupledZSparse& lhs)
   {
      this->mLHSMatrix.at(idx) = lhs.first.cast<MHDComplex>() + MathConstants::cI*lhs.second;
   }

   void EquationZTimestepper::setRHSMatrix(const int idx, const DecoupledZSparse& rhs)
   {
      this->mRHSMatrix.at(idx) = rhs.first.cast<MHDComplex>() + MathConstants::cI*rhs.second;
   }

   void EquationZTimestepper::completeLHSMatrix(const int idx, const DecoupledZSparse& lhs)
   {
      this->mLHSMatrix.at(idx) += lhs.first.cast<MHDComplex>() + MathConstants::cI*lhs.second;
   }

   void EquationZTimestepper::completeRHSMatrix(const int idx, const DecoupledZSparse& rhs)
   {
      this->mRHSMatrix.at(idx) += rhs.first.cast<MHDComplex>() + MathConstants::cI*rhs.second;
   }

   void EquationZTimestepper::addStorage(const int cols)
   {
      // get the required number of rows
      int rows = this->mLHSMatrix.back().cols();

      // Add storage for RHS data
      this->mRHSData.push_back(MatrixZ(rows,cols));
      this->mRHSData.back().setZero();

      // Add storage for old RHS
      this->mRHSOld.push_back(MatrixZ(rows,cols));
      this->mRHSOld.back().setZero();

      // Add storage for solution
      this->mSolution.push_back(MatrixZ(rows,cols));
      this->mSolution.back().setZero();
   }
}
}
