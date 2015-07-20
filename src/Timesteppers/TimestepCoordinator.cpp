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
      : Solver::SparseLinearCoordinatorBase<SparseTimestepper>(), mcMaxJump(1.602), mcUpWindow(1.05), mcMinDt(1e-11), mcMaxDt(1e-1), mOldDt(this->mcMinDt), mDt(this->mcMinDt), mTime(0.0), mCnstSteps(0.0), mStepTime(0.0)
   {
      this->mNStep = IntegratorSelector::STEPS;

      // Create CFL writer
      IoAscii::SharedCflWriter   spCflWriter = IoAscii::SharedCflWriter(new IoAscii::CflWriter());
      this->mspIo = spCflWriter;
      this->mspIo->init();
   }

   TimestepCoordinator::~TimestepCoordinator()
   {
      this->mspIo->finalize();
   }

   void TimestepCoordinator::tuneAdaptive(const MHDFloat time)
   {
      this->mStepTime = time;
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

         DebuggerMacro_start("Complex operator update", 0);
         // Update solvers from complex operator, complex field steppers
         Solver::updateSolvers<SparseTimestepper, Solver::SparseCoordinatorBase<SparseTimestepper>::ComplexSolver_iterator>(*this);
         DebuggerMacro_stop("Complex operator solver update t = ", 0);

         DebuggerMacro_start("Real operator solver update", 0);
         // Update solvers from real operator, complex field steppers
         Solver::updateSolvers<SparseTimestepper, Solver::SparseCoordinatorBase<SparseTimestepper>::RealSolver_iterator>(*this);
         DebuggerMacro_stop("Real operator solver update t = ", 0);

         // Reset the step index
         this->mStep = 0;
      } else
      {
         this->mCnstSteps += 1.0;
      }

      // Update CFL writer
      this->mspIo->setSimTime(this->mTime, this->mDt, this->mCnstSteps);
      this->mspIo->write();

      if(hasNewDt)
      {
         this->mCnstSteps = 0.0;
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

      // Update the internal step counter
      this->updateStep();

      // Clear the solver RHS
      this->clearSolvers();
   }

   void TimestepCoordinator::init(const MHDFloat time, const MHDFloat dt, const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
   {
      // Set initial time
      this->mTime = time;

      // Set initial timestep
      this->mOldDt = dt;
      this->mDt = dt;
      DebuggerMacro_showValue("Creating timestepper with initial timestep Dt = ", 0, this->mDt);

      // Initialise solver
      Solver::SparseLinearCoordinatorBase<SparseTimestepper>::init(scalEq, vectEq);
   }

   void TimestepCoordinator::updateMatrices()
   {
      // Loop over all substeps of timestepper
      for(int step = 0; step < this->mNStep; ++step)
      {
         // Compute timestep correction coefficient for LHS matrix
         MHDFloat lhsCoeff = IntegratorSelector::lhsT(step)*(1.0/this->mOldDt - 1.0/this->mDt);

         // Compute timestep correction coefficient for RHS matrix at t_n
         MHDFloat rhsCoeff = IntegratorSelector::rhsT(0, step)*(1.0/this->mOldDt - 1.0/this->mDt);

         // Compute timestep correction coefficient for RHS matrix at t_(n-i), i > 0
         std::vector<MHDFloat> oldRhsCoeff;
         for(int i = 0; i < IntegratorSelector::FIELD_MEMORY; ++i)
         {
            oldRhsCoeff.push_back(IntegratorSelector::rhsT(i+1, step)*(1.0/this->mOldDt - 1.0/this->mDt));
         }

         // Loop over all complex operator, complex field timesteppers
         Solver::updateTimeMatrixSolvers<SparseTimestepper, Solver::SparseCoordinatorBase<SparseTimestepper>::ComplexSolver_iterator>(*this, lhsCoeff, rhsCoeff, oldRhsCoeff, step);

         // Loop over all real operator, complex field timesteppers
         Solver::updateTimeMatrixSolvers<SparseTimestepper, Solver::SparseCoordinatorBase<SparseTimestepper>::RealSolver_iterator>(*this, lhsCoeff, rhsCoeff, oldRhsCoeff, step);
      }
   }

   void TimestepCoordinator::computeRHS()
   {
      // Compute RHS component for complex operator, complex field linear systems
      Solver::computeRHSSolvers<SparseTimestepper, Solver::SparseCoordinatorBase<SparseTimestepper>::ComplexSolver_iterator>(*this, this->mStep);

      // Compute RHS component for real operator, complex field linear systems
      Solver::computeRHSSolvers<SparseTimestepper, Solver::SparseCoordinatorBase<SparseTimestepper>::RealSolver_iterator>(*this, this->mStep);
   }

   void TimestepCoordinator::computeTimeCoeffs(MHDFloat& lhsL, MHDFloat& lhsT, MHDFloat& rhsL, MHDFloat& rhsT, std::vector<MHDFloat>& oldRhsL, std::vector<MHDFloat>& oldRhsT)
   {
      // Set time coefficients for LHS Matrix
      lhsT = IntegratorSelector::lhsT(this->mStep)*1.0/this->mDt;

      // Set linear coefficients for LHS Matrix
      lhsL = IntegratorSelector::lhsL(this->mStep);

      // Set time coefficients for RHS Matrix
      rhsT = IntegratorSelector::rhsT(0, this->mStep)*1.0/this->mDt;

      // Set linear coefficients for RHS Matrix
      rhsL = IntegratorSelector::rhsL(0, this->mStep);

      // Compute timestep correction coefficient for RHS matrix at t_(n-i), i > 0
      for(int i = 0; i < IntegratorSelector::FIELD_MEMORY; ++i)
      {
         oldRhsT.push_back(IntegratorSelector::rhsT(i+1, this->mStep)*1.0/this->mDt);
         oldRhsL.push_back(IntegratorSelector::rhsL(i+1, this->mStep));
      }
   }

   void TimestepCoordinator::buildSolverMatrix(TimestepCoordinator::SharedRealSolverType spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      // Operator coefficients
      MHDFloat lhsLCoeff, lhsTCoeff, rhsLCoeff, rhsTCoeff;
      std::vector<MHDFloat>   oldRhsLCoeff;
      std::vector<MHDFloat>   oldRhsTCoeff;

      this->computeTimeCoeffs(lhsLCoeff, lhsTCoeff, rhsLCoeff, rhsTCoeff, oldRhsLCoeff, oldRhsTCoeff);
      
      buildSolverMatrixWrapper(spSolver, matIdx, spEq, comp, idx, lhsLCoeff, lhsTCoeff, rhsLCoeff, rhsTCoeff, oldRhsLCoeff, oldRhsTCoeff);
   }

   void TimestepCoordinator::buildSolverMatrix(TimestepCoordinator::SharedComplexSolverType spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      // Operator coefficients
      MHDFloat lhsLCoeff, lhsTCoeff, rhsLCoeff, rhsTCoeff;
      std::vector<MHDFloat>   oldRhsLCoeff;
      std::vector<MHDFloat>   oldRhsTCoeff;

      this->computeTimeCoeffs(lhsLCoeff, lhsTCoeff, rhsLCoeff, rhsTCoeff, oldRhsLCoeff, oldRhsTCoeff);

      buildSolverMatrixWrapper(spSolver, matIdx, spEq, comp, idx, lhsLCoeff, lhsTCoeff, rhsLCoeff, rhsTCoeff, oldRhsLCoeff, oldRhsTCoeff);
   }

}
}
