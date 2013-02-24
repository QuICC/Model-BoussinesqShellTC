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
         for(int i = 0; i < this->mRHSData.size(); i++)
         {
            this->mRHSOld.at(i) = this->mRHSData.at(i);
            this->mRHSData.at(i) = this->mRHSMatrix.at(i+start)*this->mSolution.at(i) + ImExRK3::rhsN(step)*this->mRHSData.at(i);
         }
      } else
      {
         for(int i = 0; i < this->mRHSData.size(); i++)
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

      // Set l = 0 mode to zero
      this->mSolution.at(0).setZero();

      // Create real solver for m=0
      // Eigen::SuperLU<SparseMatrix>  solverRe;
      // SparseMatrix   tmp = this->mLHSMatrix.at(0+start).real();
      // solverRe.compute(tmp);
      // Matrix rhs;
      //
      // rhs = this->mRHSData.at(0).real();
      // this->mSolution.at(0).real() = solverRe.solve(rhs);
      // rhs = this->mRHSData.at(0).imag();
      // this->mSolution.at(0).imag() = solverRe.solve(rhs);

      for(int i = 1; i < this->mRHSData.size(); i++)
      {
         this->mSolution.at(i) = this->mSolver.at(i+start)->solve(this->mRHSData.at(i));
         //std::cerr << (this->mSolution.at(i) - dealSol).array().abs().matrix().transpose() << std::endl;
         //std::cerr << this->mSolution.at(i).transpose() << std::endl;
      }
   }

   void EquationZTimestepper::initSolver()
   {
      this->mSolver.reserve(this->mLHSMatrix.size());
      for(int i = 0; i < this->mLHSMatrix.size(); i++)
      {
         SharedPtrMacro<Eigen::SuperLU<SparseMatrixZ> >  solver(new Eigen::SuperLU<SparseMatrixZ>());

         this->mSolver.push_back(solver);
      }

      for(int i = 0; i < this->mLHSMatrix.size(); i++)
      {
         if(i % this->nSystem() != 0)
         {
            std::cerr << " --------- i = "<< i << " ---------------" << std::endl;
            std::cerr << this->mLHSMatrix.at(i) << std::endl;
            this->mSolver.at(i)->compute(this->mLHSMatrix.at(i));
         }
      }
   }

   void EquationZTimestepper::addLHSMatrix(const DecoupledZSparse& lhs)
   {
      this->mLHSMatrix.push_back(lhs.first.cast<MHDComplex>() + MathConstants::cI*lhs.second);
   }

   void EquationZTimestepper::addRHSMatrix(const DecoupledZSparse& rhs)
   {
      this->mRHSMatrix.push_back(rhs.first.cast<MHDComplex>() + MathConstants::cI*rhs.second);
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
