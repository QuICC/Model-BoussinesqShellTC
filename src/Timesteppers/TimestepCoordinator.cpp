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
#include "TypeSelectors/TimeSchemeSelector.hpp"

namespace GeoMHDiSCC {

namespace Timestep {

   TimestepCoordinator::TimestepCoordinator()
      : SparseLinearCoordinatorBase<SparseTimestepper>(), mcMaxJump(1.602), mcUpWindow(1.05), mcMinDt(1e-8), mOldDt(this->mcMinDt), mDt(this->mcMinDt), mTime(0.0)
   {
      this->mNStep = IntegratorSelector::STEPS;
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

         DebuggerMacro_start("Complex operator, complex field solver update", 0);
         // Update solvers from complex operator, complex field steppers
         SolverZZ_iterator   solZIt;
         for(solZIt = this->mZZSolvers.begin(); solZIt != this->mZZSolvers.end(); ++solZIt)
         {
            (*solZIt)->updateSolver();
         }
         DebuggerMacro_stop("Complex operator, complex field solver update t = ", 0);

         DebuggerMacro_start("Real operator, real field solver update", 0);
         // Update solvers from real operator, real field steppers
         SolverRR_iterator   solRRIt;
         for(solRRIt = this->mRRSolvers.begin(); solRRIt != this->mRRSolvers.end(); ++solRRIt)
         {
            (*solRRIt)->updateSolver();
         }
         DebuggerMacro_stop("Real operator, real field solver update t = ", 0);

         DebuggerMacro_start("Real operator, complex field solver update", 0);
         // Update solvers from real operator, complex field steppers
         SolverRZ_iterator   solRZIt;
         for(solRZIt = this->mRZSolvers.begin(); solRZIt != this->mRZSolvers.end(); ++solRZIt)
         {
            (*solRZIt)->updateSolver();
         }
         DebuggerMacro_stop("Real operator, complex field solver update t = ", 0);

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

   void TimestepCoordinator::updateMatrices()
   {
      // Loop over all substeps of timestepper
      for(int step = 0; step < this->mNStep; ++step)
      {
         // Compute timestep correction coefficient for LHS matrix
         MHDFloat lhsCoeff = IntegratorSelector::lhsT(step)*(1.0/this->mOldDt - 1.0/this->mDt);

         // Compute timestep correction coefficient for RHS matrix
         MHDFloat rhsCoeff = IntegratorSelector::rhsT(step)*(1.0/this->mOldDt - 1.0/this->mDt);

         // Loop over all complex operator, complex field timesteppers
         SolverZZ_iterator   solZIt;
         for(solZIt = this->mZZSolvers.begin(); solZIt != this->mZZSolvers.end(); ++solZIt)
         {
            (*solZIt)->updateTimeMatrix(lhsCoeff, rhsCoeff, step);
         }

         // Loop over all real operator, real field timesteppers
         SolverRR_iterator   solRRIt;
         for(solRRIt = this->mRRSolvers.begin(); solRRIt != this->mRRSolvers.end(); ++solRRIt)
         {
            (*solRRIt)->updateTimeMatrix(lhsCoeff, rhsCoeff, step);
         }

         // Loop over all real operator, complex field timesteppers
         SolverRZ_iterator   solRZIt;
         for(solRZIt = this->mRZSolvers.begin(); solRZIt != this->mRZSolvers.end(); ++solRZIt)
         {
            (*solRZIt)->updateTimeMatrix(lhsCoeff, rhsCoeff, step);
         }
      }
   }

   void TimestepCoordinator::computeRHS()
   {
      // Compute RHS component for complex operator, complex field linear systems
      SolverZZ_iterator   solZIt;
      for(solZIt = this->mZZSolvers.begin(); solZIt != this->mZZSolvers.end(); ++solZIt)
      {
         // Compute linear solve RHS
         (*solZIt)->computeRHS(this->mStep);
      }

      // Compute RHS component for real operator, real field linear systems
      SolverRR_iterator   solRRIt;
      for(solRRIt = this->mRRSolvers.begin(); solRRIt != this->mRRSolvers.end(); ++solRRIt)
      {
         // Compute linear solve RHS
         (*solRRIt)->computeRHS(this->mStep);
      }

      // Compute RHS component for real operator, complex field linear systems
      SolverRZ_iterator   solRZIt;
      for(solRZIt = this->mRZSolvers.begin(); solRZIt != this->mRZSolvers.end(); ++solRZIt)
      {
         // Compute linear solve RHS
         (*solRZIt)->computeRHS(this->mStep);
      }
   }

   void TimestepCoordinator::buildSolverMatrix(TimestepCoordinator::SharedRRSolverType spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      this->buildSolverMatrixWrapper(spSolver, matIdx, spEq, comp, idx);
   }

   void TimestepCoordinator::buildSolverMatrix(TimestepCoordinator::SharedRZSolverType spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      this->buildSolverMatrixWrapper(spSolver, matIdx, spEq, comp, idx);
   }

   void TimestepCoordinator::buildSolverMatrix(TimestepCoordinator::SharedZZSolverType spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      this->buildSolverMatrixWrapper(spSolver, matIdx, spEq, comp, idx);
   }

}
}
