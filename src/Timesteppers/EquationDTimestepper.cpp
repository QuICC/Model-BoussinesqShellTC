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

   EquationDTimestepper::EquationDTimestepper(const int nField, const int start)
      : EquationTimestepperBase(nField, start)
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

   void EquationDTimestepper::solve(const int step)
   {
      int start = step*this->nSystem();

      for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
      {
         // Solve for the real component
         this->mSolution.at(i).first = this->mSolver.at(i+start)->solve(this->mRHSData.at(i).first);

         // Safety assert for successful solve
         assert(this->mSolver.at(i+start)->info() == Eigen::Success);

         // Solve for the imaginary component
         this->mSolution.at(i).second = this->mSolver.at(i+start)->solve(this->mRHSData.at(i).second);

         // Safety assert for successful solve
         assert(this->mSolver.at(i+start)->info() == Eigen::Success);
      }
   }

   void EquationDTimestepper::initSolver()
   {
      // Initialise the RHS matrices
      this->mRHSMatrix.reserve(this->mRHSTriplets.size());
      for(size_t i = 0; i < this->mRHSTriplets.size(); i++)
      {
         // Get matrix size
         int rows = this->mSize.at(i % this->nSystem()); 

         // Create RHS matrix from triplets
         SparseMatrix tmp(rows, rows);
         tmp.setFromTriplets(this->mRHSTriplets.at(i).begin(), this->mRHSTriplets.at(i).end());
         
         // Add RHS matrix to storage
         this->mRHSMatrix.push_back(tmp);
      }

      // Initialise the solver
      this->mSolver.reserve(this->mLHSTriplets.size());
      for(size_t i = 0; i < this->mLHSTriplets.size(); i++)
      {
         SharedPtrMacro<SparseSolverMacro<SparseMatrix> >  solver(new SparseSolverMacro<SparseMatrix>());

         this->mSolver.push_back(solver);
      }

      // Compute the pattern and the factorisations
      for(size_t i = 0; i < this->mLHSTriplets.size(); i++)
      {
         if(static_cast<int>(i) % this->nSystem() >= this->mZeroIdx)
         {
            // Get matrix size
            int rows = this->mSize.at(i % this->nSystem()); 

            // Create matrix from stored triplets
            SparseMatrix tmp(rows, rows);
            tmp.setFromTriplets(this->mLHSTriplets.at(i).begin(), this->mLHSTriplets.at(i).end());

            // Safety assert to make sur matrix is compressed
            assert(tmp.isCompressed());

            this->mSolver.at(i)->compute(tmp);

            // Safety assert for successful factorisation
            assert(this->mSolver.at(i)->info() == Eigen::Success);
         }
      }
   }

   void EquationDTimestepper::updateSolver()
   {
      // Update the RHS matrices
      for(size_t i = 0; i < this->mRHSTriplets.size(); i++)
      {
         this->mRHSMatrix.at(i).setFromTriplets(this->mRHSTriplets.at(i).begin(), this->mRHSTriplets.at(i).end());
      }

      // Compute the factorisations
      for(size_t i = 0; i < this->mLHSTriplets.size(); i++)
      {
         if(static_cast<int>(i) % this->nSystem() >= this->mZeroIdx)
         {
            // Get matrix size
            int rows = this->mSize.at(i % this->nSystem()); 

            // Create matrix from stored triplets
            SparseMatrix tmp(rows, rows);
            tmp.setFromTriplets(this->mLHSTriplets.at(i).begin(), this->mLHSTriplets.at(i).end());

            // Safety assert to make sur matrix is compressed
            assert(tmp.isCompressed());

            this->mSolver.at(i)->factorize(tmp);

            // Safety assert for successful factorisation
            assert(this->mSolver.at(i)->info() == Eigen::Success);
         }
      }
   }

   void EquationDTimestepper::reserveMatrices(const int n)
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
         this->mLHSTriplets.push_back(std::vector<Triplet>());
         this->mRHSTriplets.push_back(std::vector<Triplet>());
      }
   }

   std::vector<Triplet>& EquationDTimestepper::rLHSTriplets(const int idx)
   {
      return this->mLHSTriplets.at(idx);
   }

   std::vector<Triplet>& EquationDTimestepper::rRHSTriplets(const int idx)
   {
      return this->mRHSTriplets.at(idx);
   }

   void EquationDTimestepper::addStorage(const int rows, const int cols)
   {
      // Assert for non zero rows and columns
      assert(rows > 0);
      assert(cols > 0);

      // Set the matrix size
      this->mSize.push_back(rows);

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
   }
}
}
