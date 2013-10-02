/** 
 * @file TimestepCoordinator.cpp
 * @brief Implementation of a general timestep coordinator structure
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Debug includes
//
#include "Debug/DebuggerMacro.h"
#include "Profiler/ProfilerMacro.h"

// System includes
//

// External includes
//
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

// Class include
//
#include "Timesteppers/TimestepCoordinator.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"
#include "TypeSelectors/TimeSchemeSelector.hpp"

namespace GeoMHDiSCC {

namespace Timestep {

   TimestepCoordinator::TimestepCoordinator()
      : SparseLinearCoordinatorBase(), mcMaxJump(1.602), mcUpWindow(1.05), mcMinDt(1e-8), mOldDt(this->mcMinDt), mDt(this->mcMinDt), mTime(0.0)
   {
      this->mNStep = TimeSchemeType::STEPS;
   }

   TimestepCoordinator::~TimestepCoordinator()
   {
   }

   MHDFloat TimestepCoordinator::time() const
   {
      return this->mTime;
   }

   MHDFloat TimestepCoordinator::timestep() const
   {
      return this->mDt;
   }

   void TimestepCoordinator::update()
   {
      this->mTime += this->mDt;
   }

   void TimestepCoordinator::adaptTimestep(const MHDFloat cfl, const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
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
         SolverRZ_iterator   solRZIt;
         for(solRZIt = this->mRZSolvers.begin(); solRZIt != this->mRZSolvers.end(); ++solRZIt)
         {
            (*solRZIt)->updateSolver();
         }
         DebuggerMacro_stop("Real solver update t = ", 0);

         // Reset the step index
         this->mStep = 0;
      }
   }

   void TimestepCoordinator::stepForward(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq, const ScalarVariable_map& scalVar, const VectorVariable_map& vectVar)
   {

      DetailedProfilerMacro_start(ProfilerMacro::TSTEPIN);
      // Update the equation input to the timestepper
      this->getInput(scalEq, vectEq, scalVar, vectVar);
      DetailedProfilerMacro_stop(ProfilerMacro::TSTEPIN);

      DetailedProfilerMacro_start(ProfilerMacro::TSTEPRHS);
      // Compute the RHS of the linear systems
      this->computeRHS();
      DetailedProfilerMacro_stop(ProfilerMacro::TSTEPRHS);

      DetailedProfilerMacro_start(ProfilerMacro::TSTEPSOLVE);
      // Solve all the linear systems
      this->solveSystems();
      DetailedProfilerMacro_stop(ProfilerMacro::TSTEPSOLVE);
      
      DetailedProfilerMacro_start(ProfilerMacro::TSTEPOUT);
      // Transfer timestep output back to equations
      this->transferOutput(scalEq, vectEq);
      DetailedProfilerMacro_stop(ProfilerMacro::TSTEPOUT);

      // Update the internal step counter, counting from 0 to steps - 1
      this->mStep = (this->mStep + 1) % this->mNStep;
   }

   void TimestepCoordinator::init(const MHDFloat dt, const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
   {
      // Set initial timestep
      this->mOldDt = dt;
      this->mDt = dt;
      DebuggerMacro_showValue("Creating timestepper with initial timestep Dt = ", 0, this->mDt);

      // Initialise solver
      SparseLinearCoordinatorBase::init(scalEq, vectEq);
   }

   void TimestepCoordinator::addSolverRZ(const int start)
   {
      SharedSparseRZTimestepper spSolver(new SparseRZTimestepper(start));

      this->mRZSolvers.push_back(spSolver);
   }

   void TimestepCoordinator::addSolverZ(const int start)
   {
      SharedSparseZTimestepper spSolver(new SparseZTimestepper(start));

      this->mZSolvers.push_back(spSolver);
   }

   void TimestepCoordinator::updateMatrices()
   {
      // Loop over all substeps of timestepper
      for(int step = 0; step < this->mNStep; ++step)
      {
         // Compute timestep correction coefficient for LHS matrix
         MHDFloat lhsCoeff = TimeSchemeType::lhsT(step)*(1.0/this->mOldDt - 1.0/this->mDt);

         // Compute timestep correction coefficient for RHS matrix
         MHDFloat rhsCoeff = TimeSchemeType::rhsT(step)*(1.0/this->mOldDt - 1.0/this->mDt);

         // Loop over all complex timesteppers
         SolverZ_iterator   solZIt;
         for(solZIt = this->mZSolvers.begin(); solZIt != this->mZSolvers.end(); ++solZIt)
         {
            std::tr1::static_pointer_cast<SparseZTimestepper>(*solZIt)->updateTimeMatrix(lhsCoeff, rhsCoeff, step);
         }

         // Loop over all real timesteppers
         SolverRZ_iterator   solRZIt;
         for(solRZIt = this->mRZSolvers.begin(); solRZIt != this->mRZSolvers.end(); ++solZIt)
         {
            std::tr1::static_pointer_cast<SparseRZTimestepper>(*solRZIt)->updateTimeMatrix(lhsCoeff, rhsCoeff, step);
         }
      }
   }

   void TimestepCoordinator::computeRHS()
   {
      // Compute RHS component for complex linear systems
      SolverZ_iterator   solZIt;
      for(solZIt = this->mZSolvers.begin(); solZIt != this->mZSolvers.end(); ++solZIt)
      {
         // Compute linear solve RHS
         std::tr1::static_pointer_cast<SparseZTimestepper>(*solZIt)->computeRHS(this->mStep);
      }

      SolverRZ_iterator   solRZIt;
      for(solRZIt = this->mRZSolvers.begin(); solRZIt != this->mRZSolvers.end(); ++solRZIt)
      {
         // Compute linear solve RHS
         std::tr1::static_pointer_cast<SparseRZTimestepper>(*solRZIt)->computeRHS(this->mStep);
      }
   }

   void TimestepCoordinator::buildSolverMatrix(Solver::SharedSparseRZLinearSolver spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      // Operator coefficients
      MHDFloat lhsTCoeff;
      MHDFloat rhsTCoeff;
      MHDFloat lhsLCoeff;
      MHDFloat rhsLCoeff;

      // Set time coefficients for LHS Matrix
      lhsTCoeff = TimeSchemeType::lhsT(this->mStep)*1.0/this->mDt;

      // Set linear coefficients for LHS Matrix
      lhsLCoeff = TimeSchemeType::lhsL(this->mStep);

      // Set time coefficients for RHS Matrix
      rhsTCoeff = TimeSchemeType::rhsT(this->mStep)*1.0/this->mDt;

      // Set linear coefficients for RHS Matrix
      rhsLCoeff = -TimeSchemeType::rhsL(this->mStep);

      // Resize LHS matrix if necessary
      if(spSolver->rLHSMatrix(matIdx).size() == 0)
      {
         spSolver->rLHSMatrix(matIdx).resize(spEq->couplingInfo(comp).systemN(idx), spEq->couplingInfo(comp).systemN(idx));
      }

      // Resize RHS matrix if necessary
      if(std::tr1::static_pointer_cast<SparseRZTimestepper>(spSolver)->rRHSMatrix(matIdx).size() == 0)
      {
         std::tr1::static_pointer_cast<SparseRZTimestepper>(spSolver)->rRHSMatrix(matIdx).resize(spEq->couplingInfo(comp).systemN(idx), spEq->couplingInfo(comp).systemN(idx));
      }

      DecoupledZSparse linRow = spEq->operatorRow(Equations::IEquation::LINEARROW, comp, idx);
      DecoupledZSparse tRow = spEq->operatorRow(Equations::IEquation::TIMEROW, comp, idx);

      // Compute LHS matrix
      spSolver->rLHSMatrix(matIdx) += spEq->operatorRow(Equations::IEquation::BOUNDARYROW, comp, idx).real();

      spSolver->rLHSMatrix(matIdx) += lhsLCoeff*linRow.real() - lhsTCoeff*tRow.real();

      // Compute RHS matrix
      std::tr1::static_pointer_cast<SparseRZTimestepper>(spSolver)->rRHSMatrix(matIdx) += rhsLCoeff*linRow.real() - rhsTCoeff*tRow.real();

      // Set time matrix for timestep updates
      if(matIdx == idx)
      {
         // Resize time matrix if necessary
         if(std::tr1::static_pointer_cast<SparseRZTimestepper>(spSolver)->rTMatrix(idx).size() == 0)
         {
            std::tr1::static_pointer_cast<SparseRZTimestepper>(spSolver)->rTMatrix(idx).resize(spEq->couplingInfo(comp).systemN(idx), spEq->couplingInfo(comp).systemN(idx));
         }

         // Set time matrix
         std::tr1::static_pointer_cast<SparseRZTimestepper>(spSolver)->rTMatrix(idx) += tRow.real();
      }
   }

   void TimestepCoordinator::buildSolverMatrix(Solver::SharedSparseZLinearSolver spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      // Operator coefficients
      MHDFloat lhsTCoeff;
      MHDFloat rhsTCoeff;
      MHDFloat lhsLCoeff;
      MHDFloat rhsLCoeff;

      // Set time coefficients for LHS Matrix
      lhsTCoeff = TimeSchemeType::lhsT(this->mStep)*1.0/this->mDt;

      // Set linear coefficients for LHS Matrix
      lhsLCoeff = TimeSchemeType::lhsL(this->mStep);

      // Set time coefficients for RHS Matrix
      rhsTCoeff = TimeSchemeType::rhsT(this->mStep)*1.0/this->mDt;

      // Set linear coefficients for RHS Matrix
      rhsLCoeff = -TimeSchemeType::rhsL(this->mStep);

      // Resize LHS matrix if necessary
      if(spSolver->rLHSMatrix(matIdx).size() == 0)
      {
         spSolver->rLHSMatrix(matIdx).resize(spEq->couplingInfo(comp).systemN(idx), spEq->couplingInfo(comp).systemN(idx));
      }

      // Resize RHS matrix if necessary
      if(std::tr1::static_pointer_cast<SparseZTimestepper>(spSolver)->rRHSMatrix(matIdx).size() == 0)
      {
         std::tr1::static_pointer_cast<SparseZTimestepper>(spSolver)->rRHSMatrix(matIdx).resize(spEq->couplingInfo(comp).systemN(idx), spEq->couplingInfo(comp).systemN(idx));
      }

      // Add boundary row for LHS operator
      DecoupledZSparse bcRow = spEq->operatorRow(Equations::IEquation::BOUNDARYROW, comp, idx);
      spSolver->rLHSMatrix(matIdx) += bcRow.real().cast<MHDComplex>() + MathConstants::cI*bcRow.imag();

      DecoupledZSparse linRow = spEq->operatorRow(Equations::IEquation::LINEARROW, comp, idx);
      DecoupledZSparse tRow = spEq->operatorRow(Equations::IEquation::TIMEROW, comp, idx);


      // Set LHS matrix
      spSolver->rLHSMatrix(matIdx) += lhsLCoeff*linRow.real().cast<MHDComplex>() + MathConstants::cI*lhsLCoeff*linRow.imag() - lhsTCoeff*tRow.real().cast<MHDComplex>() - MathConstants::cI*lhsTCoeff*tRow.imag();

      // Set RHS matrix
      std::tr1::static_pointer_cast<SparseZTimestepper>(spSolver)->rRHSMatrix(matIdx) += rhsLCoeff*linRow.real().cast<MHDComplex>() + MathConstants::cI*rhsLCoeff*linRow.imag() - rhsTCoeff*tRow.real().cast<MHDComplex>() - MathConstants::cI*rhsTCoeff*tRow.imag();

      // Set time matrix for timestep updates
      if(matIdx == idx)
      {
         // Resize time matrix if necessary
         if(std::tr1::static_pointer_cast<SparseZTimestepper>(spSolver)->rTMatrix(idx).size() == 0)
         {
            std::tr1::static_pointer_cast<SparseZTimestepper>(spSolver)->rTMatrix(idx).resize(spEq->couplingInfo(comp).systemN(idx), spEq->couplingInfo(comp).systemN(idx));
         }

         // Set time matrix
         std::tr1::static_pointer_cast<SparseZTimestepper>(spSolver)->rTMatrix(idx) += tRow.real().cast<MHDComplex>() + MathConstants::cI*tRow.imag();
      }
   }

}
}
