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

   EquationDTimestepper::EquationDTimestepper(const int nField)
      : EquationTimestepperBase(nField)
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
         for(size_t i = 0; i < this->mRHSData.size(); i++)
         {
            this->mRHSOld.at(i).first = this->mRHSData.at(i).first;
            this->mRHSData.at(i).first = this->mRHSMatrix.at(i+start)*this->mSolution.at(i).first + ImExRK3::rhsN(step)*this->mRHSData.at(i).first;

            this->mRHSOld.at(i).second = this->mRHSData.at(i).second;
            this->mRHSData.at(i).second = this->mRHSMatrix.at(i+start)*this->mSolution.at(i).second + ImExRK3::rhsN(step)*this->mRHSData.at(i).second;
         }
      } else
      {
         for(size_t i = 0; i < this->mRHSData.size(); i++)
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

      for(size_t i = 0; i < this->mRHSData.size(); i++)
      {
         // Solve for the real component
         this->mSolution.at(i).first = this->mSolver.at(i+start)->solve(this->mRHSData.at(i).first);

         // Solve for the imaginary component
         this->mSolution.at(i).second = this->mSolver.at(i+start)->solve(this->mRHSData.at(i).second);
      }
   }

   void EquationDTimestepper::initSolver()
   {
      this->mSolver.reserve(this->mLHSMatrix.size());
      for(size_t i = 0; i < this->mLHSMatrix.size(); i++)
      {
         SharedPtrMacro<Eigen::SuperLU<SparseMatrix> >  solver(new Eigen::SuperLU<SparseMatrix>());

         this->mSolver.push_back(solver);
      }

      // Compute the pattern and the factorisations
      for(size_t i = 0; i < this->mLHSMatrix.size(); i++)
      {
         this->mSolver.at(i)->compute(this->mLHSMatrix.at(i));
      }
   }

   void EquationDTimestepper::updateSolver()
   {
      // Compute the factorisations
      for(size_t i = 0; i < this->mLHSMatrix.size(); i++)
      {
         this->mSolver.at(i)->factorize(this->mLHSMatrix.at(i));
      }
   }

   void EquationDTimestepper::addLHSMatrix(const DecoupledZSparse& lhs)
   {
      this->mLHSMatrix.push_back(lhs.first);
   }

   void EquationDTimestepper::addRHSMatrix(const DecoupledZSparse& rhs)
   {
      this->mRHSMatrix.push_back(rhs.first);
   }

   void EquationDTimestepper::setLHSMatrix(const int idx, const DecoupledZSparse& lhs)
   {
      this->mLHSMatrix.at(idx) = lhs.first;
   }

   void EquationDTimestepper::setRHSMatrix(const int idx, const DecoupledZSparse& rhs)
   {
      this->mRHSMatrix.at(idx) = rhs.first;
   }

   void EquationDTimestepper::completeLHSMatrix(const int idx, const DecoupledZSparse& lhs)
   {
      this->mLHSMatrix.at(idx) += lhs.first;
   }

   void EquationDTimestepper::completeRHSMatrix(const int idx, const DecoupledZSparse& rhs)
   {
      this->mRHSMatrix.at(idx) += rhs.first;
   }

   void EquationDTimestepper::addStorage(const int cols)
   {
      // get the required number of rows
      int rows = this->mLHSMatrix.back().cols();

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
