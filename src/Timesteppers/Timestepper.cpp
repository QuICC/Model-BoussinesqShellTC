/** \file Timestepper.cpp
 *  \brief Implementation of a general timestepper structure
 */

// Debug includes
//
#include "Debug/DebuggerMacro.h"

// System includes
//

// External includes
//
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

// Class include
//
#include "Timesteppers/Timestepper.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"

namespace GeoMHDiSCC {

namespace Timestep {

   Timestepper::Timestepper()
      : SparseLinearCoordinatorBase(), mcMaxJump(1.602), mcUpWindow(1.05), mcMinDt(1e-8), mOldDt(this->mcMinDt), mDt(this->mcMinDt), mTime(0.0)
   {
   }

   Timestepper::~Timestepper()
   {
   }

   MHDFloat Timestepper::time() const
   {
      return this->mTime;
   }

   MHDFloat Timestepper::timestep() const
   {
      return this->mDt;
   }

   void Timestepper::update()
   {
      this->mTime += this->mDt;
   }

   void Timestepper::adaptTimestep(const MHDFloat cfl, const ScalarEquationRangeType& scalEq, const VectorEquationRangeType& vectEq)
   {
      // Flag to update timestep
      bool hasNewDt = false;

      // Check if CFL allows for a larger timestep
      if(cfl > this->mcUpWindow*this->mDt)
      {
         // Activate matrices update
         hasNewDt = true;

         // Set new timestep
         this->mOldDt = this->mDt;
         this->mDt = std::min(cfl, this->mcMaxJump*this->mDt);
      
      // Check if CFL is below minimal timestep or downard jump is large
      } else if(cfl < this->mcMinDt || cfl < this->mDt/this->mcMaxJump)
      {
         // Don't update matrices
         hasNewDt = false;
 
         // Signal simulation abort
         this->mOldDt = this->mDt;
         this->mDt = -cfl;
     
      // Check if CFL requires a lower timestep
      } else if(cfl < this->mDt)
      {
         // Activate matrices update
         hasNewDt = true;

         // Set new timestep
         this->mOldDt = this->mDt;
         this->mDt = cfl/this->mcUpWindow;

      // No need to change timestep
      } else
      {
         hasNewDt = false;
      }

//      this->mDt = this->mOldDt;
//      hasNewDt = false;

      //
      // Update the timestep matrices if necessary
      //
      if(hasNewDt)
      {
         DebuggerMacro_showValue("Updating timestep and matrices with new Dt = ", 0, this->mDt);

         DebuggerMacro_start("Update matrices", 0);
         // Update the time dependence in matrices
         this->updateMatrices();
         DebuggerMacro_stop("Update matrices t = ", 0);

         DebuggerMacro_start("Complex solver update", 0);
         // Update solvers from complex equation steppers
         std::vector<EquationZTimestepper>::iterator   zIt;
         for(zIt = this->mEqZStepper.begin(); zIt != this->mEqZStepper.end(); ++zIt)
         {
            zIt->updateSolver();
         }
         DebuggerMacro_stop("Complex solver update t = ", 0);

         DebuggerMacro_start("Real solver update", 0);
         // Update solvers from real equation steppers
         std::vector<EquationDTimestepper>::iterator   rIt;
         for(rIt = this->mEqDStepper.begin(); rIt != this->mEqDStepper.end(); ++rIt)
         {
            rIt->updateSolver();
         }
         DebuggerMacro_stop("Real solver update t = ", 0);

         // Reset the step index
         this->mStep = 0;
      }
   }

   void Timestepper::stepForward(const ScalarEquationRangeType& scalEq, const VectorEquationRangeType& vectEq)
   {
      // Update the equation input to the timestepper
      this->getInput(scalEq, vectEq);

      // Compute the RHS of the linear systems
      this->computeRHS();

      // Solve all the linear systems
      this->solve();
      
      // Transfer timestep output back to equations
      this->transferOutput(scalEq, vectEq);

      // Update the internal step counter, counting from 0 to steps - 1
      this->mStep = (this->mStep + 1) % ImExRK3::STEPS;
   }

   void Timestepper::init(const MHDFloat dt, const ScalarEquationRangeType& scalEq, const VectorEquationRangeType& vectEq)
   {
      // Set initial timestep
      this->mOldDt = dt;
      this->mDt = dt;
      DebuggerMacro_showValue("Creating timestepper with initial timestep Dt = ", 0, this->mDt);

      // Initialise solver
      SparseLinearCoordinatorBase::init(scalEq, vectEq);
   }

   void Timestepper::addSolverD(const int start)
   {
      SharedSparseDTimestepper spSolver(new SparseDTimestepper(start));

      this->mDSolvers.push_back(spSolver);
   }

   void Timestepper::addSolverZ(const int start)
   {
      SharedSparseZTimestepper spSolver(new SparseZTimestepper(start));

      this->mZSolvers.push_back(spSolver);
   }

   void Timestepper::updateMatrices()
   {
      // Loop over all substeps of timestepper
      for(int step = 0; step < ImExRK3::STEPS; ++step)
      {
         // Compute timestep correction coefficient for LHS matrix
         MHDFloat lhsCoeff = ImExRK3::lhsT(step)*(1.0/this->mOldDt - 1.0/this->mDt);

         // Compute timestep correction coefficient for RHS matrix
         MHDFloat rhsCoeff = ImExRK3::rhsT(step)*(1.0/this->mOldDt - 1.0/this->mDt);

         // Loop over all complex timesteppers
         std::vector<EquationZTimestepper>::iterator   zIt;
         for(zIt = this->mEqZStepper.begin(); zIt != this->mEqZStepper.end(); ++zIt)
         {
            zIt->updateTimeMatrix(lhsCoeff, rhsCoeff, step);
         }

         // Loop over all real timesteppers
         std::vector<EquationDTimestepper>::iterator   rIt;
         for(rIt = this->mEqDStepper.begin(); rIt != this->mEqDStepper.end(); ++rIt)
         {
            rIt->updateTimeMatrix(lhsCoeff, rhsCoeff, step);
         }
      }
   }

   void Timestepper::computeRHS()
   {
      // Compute RHS component for complex linear systems
      std::vector<EquationZTimestepper>::iterator   zIt;
      for(zIt = this->mEqZStepper.begin(); zIt != this->mEqZStepper.end(); ++zIt)
      {
         // Compute linear solve RHS
         zIt->computeRHS(this->mStep);
      }

      std::vector<EquationDTimestepper>::iterator   dIt;
      for(dIt = this->mEqDStepper.begin(); dIt != this->mEqDStepper.end(); ++dIt)
      {
         // Compute linear solve RHS
         dIt->computeRHS(this->mStep);
      }
   }

   void Timestepper::buildTimeMatrix(SparseMatrix& timeMatrix, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      if(timeMatrix.size() == 0)
      {
         timeMatrix.resize(spEq->couplingInfo(comp).systemN(idx), spEq->couplingInfo(comp).systemN(idx));
      }

      timeMatrix += spEq->operatorRow(Equations::IEquation::TIMEROW, comp, idx).first;
   }

   void Timestepper::buildTimeMatrix(SparseMatrixZ& timeMatrix, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      if(timeMatrix.size() == 0)
      {
         timeMatrix.resize(spEq->couplingInfo(comp).systemN(idx), spEq->couplingInfo(comp).systemN(idx));
      }

      DecoupledZSparse tRow = spEq->operatorRow(Equations::IEquation::TIMEROW, comp, idx);

      timeMatrix += tRow.first.cast<MHDComplex>() + MathConstants::cI*tRow.second;
   }

   void Timestepper::buildSolverMatrix(SolverD_iterator solDIt, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const bool isLhs)
   {
      // Operator coefficients
      MHDFloat timeCoeff;
      MHDFloat linearCoeff;

      if(isLhs)
      {
         // Set time coefficients for LHS Matrix
         timeCoeff = ImExRK3::lhsT(this->mStep)*1.0/this->mDt;

         // Set linear coefficients for LHS Matrix
         linearCoeff = ImExRK3::lhsL(this->mStep);
      } else
      {
         // Set time coefficients for RHS Matrix
         timeCoeff = ImExRK3::rhsT(this->mStep)*1.0/this->mDt;

         // Set linear coefficients for RHS Matrix
         linearCoeff = -ImExRK3::rhsL(this->mStep);
      }

      if(solverMatrix.size() == 0)
      {
         solverMatrix.resize(spEq->couplingInfo(comp).systemN(idx), spEq->couplingInfo(comp).systemN(idx));
      }

      // Add boundary row for LHS operator
      if(isLhs)
      {
         solverMatrix += spEq->operatorRow(Equations::IEquation::BOUNDARYROW, comp, idx).first;
      }

      solverMatrix += linearCoeff*spEq->operatorRow(Equations::IEquation::LINEARROW, comp, idx).first - timeCoeff*spEq->operatorRow(Equations::IEquation::TIMEROW, comp, idx).first;
   }

   void Timestepper::buildSolverMatrix(SolverZ_iterator solZIt, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      // Operator coefficients
      MHDFloat timeCoeff;
      MHDFloat linearCoeff;

      if(isLhs)
      {
         // Set time coefficients for LHS Matrix
         timeCoeff = ImExRK3::lhsT(this->mStep)*1.0/this->mDt;

         // Set linear coefficients for LHS Matrix
         linearCoeff = ImExRK3::lhsL(this->mStep);
      } else
      {
         // Set time coefficients for RHS Matrix
         timeCoeff = ImExRK3::rhsT(this->mStep)*1.0/this->mDt;

         // Set linear coefficients for RHS Matrix
         linearCoeff = -ImExRK3::rhsL(this->mStep);
      }

      if(solverMatrix.size() == 0)
      {
         solverMatrix.resize(spEq->couplingInfo(comp).systemN(idx), spEq->couplingInfo(comp).systemN(idx));
      }

      // Add boundary row for LHS operator
      if(isLhs)
      {
         DecoupledZSparse bcRow = spEq->operatorRow(Equations::IEquation::BOUNDARYROW, comp, idx);
         solverMatrix += bcRow.first.cast<MHDComplex>() + MathConstants::cI*bcRow.second;
      }

      DecoupledZSparse linRow = spEq->operatorRow(Equations::IEquation::LINEARROW, comp, idx);
      DecoupledZSparse tRow = spEq->operatorRow(Equations::IEquation::TIMEROW, comp, idx);
      solverMatrix += linearCoeff*linRow.first.cast<MHDComplex>() + MathConstants::cI*linearCoeff*linRow.second - timeCoeff*tRow.first.cast<MHDComplex>() - MathConstants::cI*timeCoeff*tRow.second;
   }

}
}
