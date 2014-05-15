/** 
 * @file TimestepRCoordinatorBase.cpp
 * @brief Implementation of a real field timestep coordinator base
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

// Class include
//
#include "Timesteppers/TimestepRCoordinatorBase.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Timestep {

   TimestepCoordinatorBase<false>::TimestepCoordinatorBase()
      : SparseLinearCoordinatorBase<SparseTimestepper,false>()
   {
   }

   TimestepCoordinatorBase<false>::~TimestepCoordinatorBase()
   {
   }

   void TimestepCoordinatorBase<false>::updateSolvers()
   {
      DebuggerMacro_start("Real operator, real field solver update", 0);
      // Update solvers from real operator, real field steppers
      SolverRR_iterator   solRRIt;
      for(solRRIt = this->mRRSolvers.begin(); solRRIt != this->mRRSolvers.end(); ++solRRIt)
      {
         (*solRRIt)->updateSolver();
      }
      DebuggerMacro_stop("Real operator, real field solver update t = ", 0);
   }

   void TimestepCoordinatorBase<false>::updateTimeMatrices(const MHDFloat lhsCoeff, const MHDFloat rhsCoeff, const int step)
   {
      // Loop over all real operator, real field timesteppers
      SolverRR_iterator   solRRIt;
      for(solRRIt = this->mRRSolvers.begin(); solRRIt != this->mRRSolvers.end(); ++solRRIt)
      {
         (*solRRIt)->updateTimeMatrix(lhsCoeff, rhsCoeff, step);
      }
   }

   void TimestepCoordinatorBase<false>::computeRHS()
   {
      // Compute RHS component for real operator, real field linear systems
      SolverRR_iterator   solRRIt;
      for(solRRIt = this->mRRSolvers.begin(); solRRIt != this->mRRSolvers.end(); ++solRRIt)
      {
         // Compute linear solve RHS
         (*solRRIt)->computeRHS(this->mStep);
      }
   }

   void TimestepCoordinatorBase<false>::buildSolverMatrix(TimestepCoordinatorBase<false>::SharedRRSolverType spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      this->buildSolverMatrixWrapper(spSolver, matIdx, spEq, comp, idx);
   }

}
}
