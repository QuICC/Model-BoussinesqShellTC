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

   EquationZTimestepper::EquationZTimestepper(const int nField, const int start)
      : EquationTimestepperBase(nField, start)
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

   void EquationZTimestepper::initSolver()
   {
      // Initialise the RHS matrices
      this->mRHSMatrix.reserve(this->mRHSTriplets.size());
      for(size_t i = 0; i < this->mRHSTriplets.size(); i++)
      {
         // Get matrix size
         int rows = this->mSize.at(i % this->nSystem()); 

         // Create RHS matrix from triplets
         SparseMatrixZ tmp(rows, rows);
         tmp.setFromTriplets(this->mRHSTriplets.at(i).begin(), this->mRHSTriplets.at(i).end());
         tmp.makeCompressed();
         
         // Add RHS matrix to storage
         this->mRHSMatrix.push_back(tmp);
      }

      // Initialise solver
      this->mSolver.reserve(this->mLHSTriplets.size());
      for(size_t i = 0; i < this->mLHSTriplets.size(); i++)
      {
         SharedPtrMacro<SparseSolverMacro<SparseMatrixZ> >  solver(new SparseSolverMacro<SparseMatrixZ>());

         this->mSolver.push_back(solver);
      }

      // Compute pattern and factorisation
      for(size_t i = 0; i < this->mLHSTriplets.size(); i++)
      {
         if(static_cast<int>(i) % this->nSystem() >= this->mZeroIdx)
         {
            // Get matrix size
            int rows = this->mSize.at(i % this->nSystem()); 

            // Create matrix from stored triplets
            SparseMatrixZ tmp(rows, rows);
            tmp.setFromTriplets(this->mLHSTriplets.at(i).begin(), this->mLHSTriplets.at(i).end());
            tmp.makeCompressed();

            // Safety assert to make sur matrix is compressed
            assert(tmp.isCompressed());

            this->mSolver.at(i)->compute(tmp);

            // Safety assert for successful factorisation
            assert(this->mSolver.at(i)->info() == Eigen::Success);
         }
      }
   }

   void EquationZTimestepper::updateSolver()
   {
      // Update the RHS matrices
      for(size_t i = 0; i < this->mRHSTriplets.size(); i++)
      {
         this->mRHSMatrix.at(i).setFromTriplets(this->mRHSTriplets.at(i).begin(), this->mRHSTriplets.at(i).end());
      }

      // Compute factorisation
      for(size_t i = 0; i < this->mLHSTriplets.size(); i++)
      {
         if(static_cast<int>(i) % this->nSystem() >= this->mZeroIdx)
         {
            // Get matrix size
            int rows = this->mSize.at(i % this->nSystem()); 

            // Create matrix from stored triplets
            SparseMatrixZ tmp(rows, rows);
            tmp.setFromTriplets(this->mLHSTriplets.at(i).begin(), this->mLHSTriplets.at(i).end());

            // Safety assert to make sur matrix is compressed
            assert(tmp.isCompressed());

            this->mSolver.at(i)->factorize(tmp);

            // Safety assert for successful factorisation
            assert(this->mSolver.at(i)->info() == Eigen::Success);
         }
      }
   }

   void EquationZTimestepper::reserveMatrices(const int n)
   {
      // Reserve space for the RHS matrices
      this->mRHSMatrix.reserve(n);

      // Reserve space for the LHS triplets
      this->mLHSTriplets.reserve(n);

      // Reserve space for the RHS triplets
      this->mRHSTriplets.reserve(n);

      // Initialise storage
      for(int i = 0; i < n; ++i)
      {
         this->mLHSTriplets.push_back(std::vector<TripletZ>());
         this->mRHSTriplets.push_back(std::vector<TripletZ>());
      }
   }

   std::vector<TripletZ>& EquationZTimestepper::rLHSTriplets(const int idx)
   {
      return this->mLHSTriplets.at(idx);
   }

   std::vector<TripletZ>& EquationZTimestepper::rRHSTriplets(const int idx)
   {
      return this->mRHSTriplets.at(idx);
   }

   void EquationZTimestepper::addStorage(const int rows, const int cols)
   {
      // Assert for non zero rows and columns
      assert(rows > 0);
      assert(cols > 0);

      // Set the matrix size
      this->mSize.push_back(rows);

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
