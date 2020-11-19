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
#include "IoTools/Formatter.hpp"

namespace QuICC {

namespace Timestep {

   TimestepCoordinator::TimestepCoordinator()
      : Solver::SparseLinearCoordinatorBase<TimeSchemeTypeSelector>(), mcMinCnst(2), mcMaxJump(1.602), mcUpWindow(1.05), mcMinDt(1e-11), mcMaxDt(1e-1), mMaxError(-1.0), mOldDt(this->mcMinDt), mDt(2,1), mTime(0.0), mRefTime(0.0), mCnstSteps(0.0), mStepTime(0.0)
   {
      this->mDt(0,0) = this->mcMinDt;
      this->mDt(1,0) = -100.0;

      // Initialize timestepper
      TimeSchemeSelector::init();

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
      return this->mDt(0,0);
   }

   void TimestepCoordinator::update()
   {
      this->mTime = this->mRefTime + this->timestep();
      this->mRefTime = this->mTime;
   }

   void TimestepCoordinator::adaptTimestep(const Matrix& cfl, const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
   {
      // Store old timestep
      this->mOldDt = this->timestep();

      // Update CFL information
      this->mDt.block(0, 1, this->mDt.rows(), cfl.cols()-1) = cfl.rightCols(cfl.cols()-1);

      // New compute CFL
      MHDFloat compCfl = cfl(0,0);

      // Check if CFL allows for a larger timestep
      MHDFloat newCflDt = 0.0;
      if(compCfl > this->mcUpWindow*this->timestep())
      {
         if(this->mCnstSteps >= this->mcMinCnst)
         {
            // Set new timestep
            newCflDt = std::min(compCfl, this->mcMaxJump*this->timestep());
         } else
         {
            // Reuse same timestep
            newCflDt = this->timestep();
         }
      
      // Check if CFL is below minimal timestep or downard jump is large
      } else if(compCfl < this->mcMinDt || compCfl < this->timestep()/this->mcMaxJump)
      {
         // Signal simulation abort
         newCflDt = -compCfl;
     
      // Check if CFL requires a lower timestep
      } else if(compCfl < this->timestep()*(2.0-this->mcUpWindow))
      {
         // Set new timestep
         newCflDt = compCfl;

      } else
      {
         newCflDt = this->timestep();
      }
      
      // Gather error across processes
      #ifdef QUICC_MPI
      if(this->mError > 0.0)
      {
         MPI_Allreduce(MPI_IN_PLACE, &this->mError, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      }
      #endif //QUICC_MPI
      
      // No error control and no CFL condition
      MHDFloat newErrorDt = 0.0;

      // Use what ever condition is used by CFL
      if(this->mError < 0)
      {
         newErrorDt = -1.0;

      // Error is too large, reduce timestep
      } else if(this->mError > this->mMaxError)
      {
         newErrorDt = this->timestep()*std::pow(this->mMaxError/this->mError,1./TimeSchemeSelector::ORDER)/this->mcUpWindow;

      // Error is small, increase timestep
      } else if(this->mError < this->mMaxError/(this->mcMaxJump*0.9) && this->mCnstSteps >= this->mcMinCnst)
      {
         newErrorDt = std::min(this->timestep()*std::pow(this->mMaxError/this->mError,1./TimeSchemeSelector::ORDER), this->timestep()*this->mcMaxJump);

      // Timestep should not be increased
      } else
      {
         newErrorDt = this->timestep();
      }

      // Update error details
      if(this->mMaxError > 0.0)
      {
         this->mDt(0, this->mDt.cols()-1) = newErrorDt;
         this->mDt(1, this->mDt.cols()-1) = this->mError;
      }

      // CFL condition requested abort!
      if(newCflDt < 0.0)
      {
         this->mDt(0,0) = newCflDt;
         this->mDt(1,0) = cfl(1,0);

      // Get minimum between both conditions
      } else if(newCflDt > 0.0 && newErrorDt > 0.0)
      {
         if(newCflDt < newErrorDt)
         {
            this->mDt(0,0) = newCflDt;
            this->mDt(1,0) = cfl(1,0);
         } else
         {
            this->mDt(0,0) = newErrorDt;
            this->mDt(1,0) = -200.0;
         }

      // Use CFL condition
      } else if(newCflDt > 0.0)
      {
         if(this->timestep() != newCflDt)
         {
            this->mDt(0,0) = newCflDt;
            this->mDt(1,0) = cfl(1,0);
         }

      // Use error condition
      } else if(newErrorDt > 0.0)
      {
         this->mDt(0,0) = newErrorDt;
         this->mDt(1,0) = -200.0;
      }

      //
      // Update the timestep matrices if necessary
      //
      if(this->timestep() != this->mOldDt && this->timestep() > 0.0)
      {
         DebuggerMacro_showValue("Updating timestep and matrices with new Dt = ", 0, this->timestep());

         DebuggerMacro_start("Update matrices", 0);
         // Update the time dependence in matrices
         this->updateMatrices();
         DebuggerMacro_stop("Update matrices t = ", 0);

         DebuggerMacro_start("Complex operator update", 0);
         // Update solvers from complex operator, complex field steppers
         Solver::updateSolvers<TimeSchemeTypeSelector, Solver::SparseCoordinatorBase<TimeSchemeTypeSelector>::ComplexSolver_iterator>(*this);
         DebuggerMacro_stop("Complex operator solver update t = ", 0);

         DebuggerMacro_start("Real operator solver update", 0);
         // Update solvers from real operator, complex field steppers
         Solver::updateSolvers<TimeSchemeTypeSelector, Solver::SparseCoordinatorBase<TimeSchemeTypeSelector>::RealSolver_iterator>(*this);
         DebuggerMacro_stop("Real operator solver update t = ", 0);
      } else
      {
         this->mCnstSteps += 1.0;
      }

      // Update CFL writer
      this->mspIo->setSimTime(this->mTime, this->mDt, this->mCnstSteps);
      this->mspIo->write();

      if(this->timestep() != this->mOldDt && this->timestep() > 0.0)
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

      DetailedProfilerMacro_start(ProfilerMacro::TSTEPSOLVE);
      // Solve all the linear systems
      this->solveSystems();
      DetailedProfilerMacro_stop(ProfilerMacro::TSTEPSOLVE);
      
      DetailedProfilerMacro_start(ProfilerMacro::TSTEPOUT);
      // Transfer timestep output back to equations
      this->transferOutput(scalEq, vectEq);
      DetailedProfilerMacro_stop(ProfilerMacro::TSTEPOUT);

      // Clear the solver RHS
      this->clearSolvers();

      // Update current time
      this->mTime = this->mRefTime + this->stepFraction()*this->timestep();
   }

   void TimestepCoordinator::init(const MHDFloat time, const Matrix& cfl, const MHDFloat maxError, const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
   {
      // Set initial time
      this->mTime = time;
      this->mRefTime = this->mTime;

      // Use embedded scheme to compute error
      if(maxError > 0.0)
      {
         TimeSchemeSelector::useEmbedded();
         this->mMaxError = maxError;
      }

      // Set initial timestep
      this->mOldDt = cfl(0,0);

      // Update CFL details
      if(this->mMaxError > 0.0)
      {
         this->mDt.resize(cfl.rows(), cfl.cols()+1);
         this->mDt.leftCols(cfl.cols()) = cfl;
         this->mDt.rightCols(1)(0) = 0.0;
         this->mDt.rightCols(1)(1) = -200.0;
      } else
      {
         this->mDt = cfl;
      }

      DebuggerMacro_showValue("Creating timestepper with initial timestep Dt = ", 0, this->timestep());

      // Initialise solver
      Solver::SparseLinearCoordinatorBase<TimeSchemeTypeSelector>::init(scalEq, vectEq);
   }

   void TimestepCoordinator::updateMatrices()
   {
      // Loop over all complex operator, complex field timesteppers
      Solver::updateTimeMatrixSolvers<TimeSchemeTypeSelector, Solver::SparseCoordinatorBase<TimeSchemeTypeSelector>::ComplexSolver_iterator>(*this, this->timestep());

      // Loop over all real operator, complex field timesteppers
      Solver::updateTimeMatrixSolvers<TimeSchemeTypeSelector, Solver::SparseCoordinatorBase<TimeSchemeTypeSelector>::RealSolver_iterator>(*this, this->timestep());
   }

   void TimestepCoordinator::buildSolverMatrix(TimestepCoordinator::SharedRealSolverType spSolver, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      buildTimestepMatrixWrapper(spSolver, spEq, comp, this->timestep(), idx);
   }

   void TimestepCoordinator::buildSolverMatrix(TimestepCoordinator::SharedComplexSolverType spSolver, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      buildTimestepMatrixWrapper(spSolver, spEq, comp, this->timestep(), idx);
   }

   void TimestepCoordinator::printInfo(std::ostream& stream)
   {
      // Create nice looking ouput header
      IoTools::Formatter::printNewline(stream);
      IoTools::Formatter::printLine(stream, '-');
      IoTools::Formatter::printCentered(stream, "Timestepper information", '*');
      IoTools::Formatter::printLine(stream, '-');

      std::stringstream oss;
      int base = 20;

      // Timestep scheme
      oss << "Timestepper: " << TimeSchemeSelector::NAME << " (" << TimeSchemeSelector::ORDER << ")";
      IoTools::Formatter::printCentered(stream, oss.str(), ' ', base);
      oss.str("");

      // General linear solver
      oss << "General solver: ";
      #if defined QUICC_SPLINALG_MUMPS
         oss << "MUMPS";
      #elif defined QUICC_SPLINALG_UMFPACK
         oss << "UmfPack";
      #elif defined QUICC_SPLINALG_SPARSELU
         oss << "SparseLU";
      #else
         oss << "(unknown)";
      #endif //defined QUICC_SPLINALG_MUMPS

      IoTools::Formatter::printCentered(stream, oss.str(), ' ', base);
      oss.str("");

      // Triangular linear solver
      oss << "Triangular solver: ";
      #if defined QUICC_SPTRILINALG_SPARSELU
         oss << "SparseLU";
      #elif defined QUICC_SPTRILINALG_MUMPS
         oss << "MUMPS";
      #elif defined QUICC_SPTRILINALG_UMFPACK
         oss << "UmfPack";
      #else
         oss << "(unknown)";
      #endif //defined QUICC_SPTRILINALG_SPARSELU

      IoTools::Formatter::printCentered(stream, oss.str(), ' ', base);
      oss.str("");

      // SPD linear solver
      oss << "SPD solver: ";
      #if defined QUICC_SPSPDLINALG_SIMPLICIALLDLT
         oss << "SimplicialLDLT";
      #elif defined QUICC_SPSPDLINALG_SIMPLICIALLLT
         oss << "SimplicialLLT";
      #elif defined QUICC_SPSPDLINALG_MUMPS
         oss << "MUMPS";
      #elif defined QUICC_SPSPDLINALG_UMFPACK
         oss << "UmfPack";
      #elif defined QUICC_SPSPDLINALG_SPARSELU
         oss << "SparseLU";
      #else
         oss << "(unknown)";
      #endif //defined QUICC_SPSPDLINALG_SIMPLICIALLDLT

      IoTools::Formatter::printCentered(stream, oss.str(), ' ', base);
      oss.str("");

      IoTools::Formatter::printLine(stream, '*');
      IoTools::Formatter::printNewline(stream);
   }

}
}
