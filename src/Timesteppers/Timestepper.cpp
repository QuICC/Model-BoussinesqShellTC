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
      this->mNSteps = ImExRK3::STEPS;
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

   void Timestepper::adaptTimestep(const MHDFloat cfl, const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
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
         SolverZ_iterator   solZIt;
         for(solZIt = this->mZSolvers.begin(); solZIt != this->mZSolvers.end(); ++solZIt)
         {
            (*solZIt)->updateSolver();
         }
         DebuggerMacro_stop("Complex solver update t = ", 0);

         DebuggerMacro_start("Real solver update", 0);
         // Update solvers from real equation steppers
         SolverD_iterator   solDIt;
         for(solDIt = this->mDSolvers.begin(); solDIt != this->mDSolvers.end(); ++solDIt)
         {
            (*solDIt)->updateSolver();
         }
         DebuggerMacro_stop("Real solver update t = ", 0);

         // Reset the step index
         this->mStep = 0;
      }
   }

   void Timestepper::stepForward(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
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
      this->mStep = (this->mStep + 1) % this->mNStep;
   }

   void Timestepper::init(const MHDFloat dt, const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
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
      for(int step = 0; step < this->mNStep; ++step)
      {
         // Compute timestep correction coefficient for LHS matrix
         MHDFloat lhsCoeff = ImExRK3::lhsT(step)*(1.0/this->mOldDt - 1.0/this->mDt);

         // Compute timestep correction coefficient for RHS matrix
         MHDFloat rhsCoeff = ImExRK3::rhsT(step)*(1.0/this->mOldDt - 1.0/this->mDt);

         // Loop over all complex timesteppers
         SolverZ_iterator   solZIt;
         for(solZIt = this->mZSolvers.begin(); solZIt != this->mZSolvers.end(); ++solZIt)
         {
            (*solZIt)->updateTimeMatrix(lhsCoeff, rhsCoeff, step);
         }

         // Loop over all real timesteppers
         SolverD_iterator   solDIt;
         for(solDIt = this->mDSolvers.begin(); solDIt != this->mDSolvers.end(); ++solZIt)
         {
            (*solDIt)->updateTimeMatrix(lhsCoeff, rhsCoeff, step);
         }
      }
   }

   void Timestepper::computeRHS()
   {
      // Compute RHS component for complex linear systems
      SolverZ_iterator   solZIt;
      for(solZIt = this->mZSolvers.begin(); solZIt != this->mZSolvers.end(); ++solZIt)
      {
         // Compute linear solve RHS
         (*solZIt)->computeRHS(this->mStep);
      }

      std::vector<EquationDTimestepper>::iterator   dIt;
      for(solDIt = this->mDSolvers.begin(); solDIt != this->mDSolvers.end(); ++solDIt)
      {
         // Compute linear solve RHS
         (*solDIt)->computeRHS(this->mStep);
      }
   }

   void Timestepper::buildSolverMatrix(SharedSparseDLinearSolver spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      // Operator coefficients
      MHDFloat lhsTCoeff;
      MHDFloat rhsTCoeff;
      MHDFloat lhsLCoeff;
      MHDFloat rhsLCoeff;

      // Set time coefficients for LHS Matrix
      lhsTCoeff = ImExRK3::lhsT(this->mStep)*1.0/this->mDt;

      // Set linear coefficients for LHS Matrix
      lhsLCoeff = ImExRK3::lhsL(this->mStep);

      // Set time coefficients for RHS Matrix
      rhsTCoeff = ImExRK3::rhsT(this->mStep)*1.0/this->mDt;

      // Set linear coefficients for RHS Matrix
      rhsLCoeff = -ImExRK3::rhsL(this->mStep);

      // Resize LHS matrix if necessary
      if(spSolver->rLHSMatrix(matIdx).size() == 0)
      {
         spSolver->rLHSMatrix(matIdx).resize(spEq->couplingInfo(comp).systemN(idx), spEq->couplingInfo(comp).systemN(idx));
      }

      // Resize RHS matrix if necessary
      if(spSolver->rRHSMatrix(matIdx).size() == 0)
      {
         spSolver->rRHSMatrix(matIdx).resize(spEq->couplingInfo(comp).systemN(idx), spEq->couplingInfo(comp).systemN(idx));
      }

      DecoupledZSparse linRow = spEq->operatorRow(Equations::IEquation::LINEARROW, comp, idx);
      DecoupledZSparse tRow = spEq->operatorRow(Equations::IEquation::TIMEROW, comp, idx);

      // Compute LHS matrix
      spSolver->rLHSMatrix(matIdx) += spEq->operatorRow(Equations::IEquation::BOUNDARYROW, comp, idx).first;

      spSolver->rLHSMatrix(matIdx) += lhsLCoeff*linRow.first - lhsTCoeff*tRow.first;

      // Compute RHS matrix
      spSolver->rRHSMatrix(matIdx) += rhsLCoeff*linRow.first - rhsCoeff*tRow.first;

      // Set time matrix for timestep updates
      if(matIdx == idx)
      {
         // Resize time matrix if necessary
         if(spSolver->rTMatrix(i).size() == 0)
         {
            spSolver->rTMatrix(idx).resize(spEq->couplingInfo(comp).systemN(idx), spEq->couplingInfo(comp).systemN(idx));
         }

         // Set time matrix
         spSolver->rTMatrix(i) += tRow.first;
      }
   }

   void Timestepper::buildSolverMatrix(SharedSparseZLinearSolver spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      // Operator coefficients
      MHDFloat lhsTCoeff;
      MHDFloat rhsTCoeff;
      MHDFloat lhsLCoeff;
      MHDFloat rhsLCoeff;

      // Set time coefficients for LHS Matrix
      lhsTCoeff = ImExRK3::lhsT(this->mStep)*1.0/this->mDt;

      // Set linear coefficients for LHS Matrix
      lhsLCoeff = ImExRK3::lhsL(this->mStep);

      // Set time coefficients for RHS Matrix
      rhsTCoeff = ImExRK3::rhsT(this->mStep)*1.0/this->mDt;

      // Set linear coefficients for RHS Matrix
      rhsLCoeff = -ImExRK3::rhsL(this->mStep);

      // Resize LHS matrix if necessary
      if(spSolver->rLHSMatrix(matIdx).size() == 0)
      {
         spSolver->rLHSMatrix(matIdx).resize(spEq->couplingInfo(comp).systemN(idx), spEq->couplingInfo(comp).systemN(idx));
      }

      // Resize RHS matrix if necessary
      if(spSolver->rRHSMatrix(matIdx).size() == 0)
      {
         spSolver->rRHSMatrix(matIdx).resize(spEq->couplingInfo(comp).systemN(idx), spEq->couplingInfo(comp).systemN(idx));
      }

      // Add boundary row for LHS operator
      DecoupledZSparse bcRow = spEq->operatorRow(Equations::IEquation::BOUNDARYROW, comp, idx);
      spSolver->rLHSMatrix(matIdx) += bcRow.first.cast<MHDComplex>() + MathConstants::cI*bcRow.second;

      DecoupledZSparse linRow = spEq->operatorRow(Equations::IEquation::LINEARROW, comp, idx);
      DecoupledZSparse tRow = spEq->operatorRow(Equations::IEquation::TIMEROW, comp, idx);


      // Set LHS matrix
      spSolver->rLHSMatrix(matIdx) += lhsLCoeff*linRow.first.cast<MHDComplex>() + MathConstants::cI*lhsLCoeff*linRow.second - lhsTCoeff*tRow.first.cast<MHDComplex>() - MathConstants::cI*lhsTCoeff*tRow.second;

      // Set RHS matrix
      spSolver->rRHSMatrix(matIdx) += rhsLCoeff*linRow.first.cast<MHDComplex>() + MathConstants::cI*rhsLCoeff*linRow.second - rhsTCoeff*tRow.first.cast<MHDComplex>() - MathConstants::cI*rhsTCoeff*tRow.second;

      // Set time matrix for timestep updates
      if(matIdx == idx)
      {
         // Resize time matrix if necessary
         if(spSolver->rTMatrix.at(i).size() == 0)
         {
            spSolver->rTMatrix.at(idx).resize(spEq->couplingInfo(comp).systemN(idx), spEq->couplingInfo(comp).systemN(idx));
         }

         // Set time matrix
         spSolver->rTMatrix.at(i) += tRow.first.cast<MHDComplex>() + MathConstants::cI*tRow.second;
      }
   }

}
}
