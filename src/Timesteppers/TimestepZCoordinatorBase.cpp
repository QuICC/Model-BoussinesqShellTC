/** 
 * @file TimestepZCoordinatorBase.cpp
 * @brief Implementation of a complex field timestep coordinator structure
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
#include "Timesteppers/TimestepZCoordinatorBase.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Timestep {

   TimestepCoordinatorBase<true>::TimestepCoordinatorBase()
      : SparseLinearCoordinatorBase<SparseTimestepper,true>()
   {
   }

   TimestepCoordinatorBase<true>::~TimestepCoordinatorBase()
   {
   }

   void TimestepCoordinatorBase<true>::adaptSolvers()
   {
      DebuggerMacro_start("Complex operator, complex field solver update", 0);
      // Update solvers from complex operator, complex field steppers
      SolverZZ_iterator   solZIt;
      for(solZIt = this->mZZSolvers.begin(); solZIt != this->mZZSolvers.end(); ++solZIt)
      {
         (*solZIt)->updateSolver();
      }
      DebuggerMacro_stop("Complex operator, complex field solver update t = ", 0);

      DebuggerMacro_start("Real operator, complex field solver update", 0);
      // Update solvers from real operator, complex field steppers
      SolverRZ_iterator   solRZIt;
      for(solRZIt = this->mRZSolvers.begin(); solRZIt != this->mRZSolvers.end(); ++solRZIt)
      {
         (*solRZIt)->updateSolver();
      }
      DebuggerMacro_stop("Real operator, complex field solver update t = ", 0);
   }

   void TimestepCoordinatorBase<true>::updateTimeMatrices(const MHDFloat lhsCoeff, const MHDFloat rhsCoeff, const int step)
   {
      // Loop over all complex operator, complex field timesteppers
      SolverZZ_iterator   solZIt;
      for(solZIt = this->mZZSolvers.begin(); solZIt != this->mZZSolvers.end(); ++solZIt)
      {
         (*solZIt)->updateTimeMatrix(lhsCoeff, rhsCoeff, step);
      }

      // Loop over all real operator, complex field timesteppers
      SolverRZ_iterator   solRZIt;
      for(solRZIt = this->mRZSolvers.begin(); solRZIt != this->mRZSolvers.end(); ++solRZIt)
      {
         (*solRZIt)->updateTimeMatrix(lhsCoeff, rhsCoeff, step);
      }
   }

   void TimestepCoordinatorBase<true>::computeRHS()
   {
      // Compute RHS component for complex operator, complex field linear systems
      SolverZZ_iterator   solZIt;
      for(solZIt = this->mZZSolvers.begin(); solZIt != this->mZZSolvers.end(); ++solZIt)
      {
         // Compute linear solve RHS
         (*solZIt)->computeRHS(this->mStep);
      }

      // Compute RHS component for real operator, complex field linear systems
      SolverRZ_iterator   solRZIt;
      for(solRZIt = this->mRZSolvers.begin(); solRZIt != this->mRZSolvers.end(); ++solRZIt)
      {
         // Compute linear solve RHS
         (*solRZIt)->computeRHS(this->mStep);
      }
   }

   void TimestepCoordinatorBase<true>::buildSolverMatrix(TimestepCoordinatorBase<true>::SharedRZSolverType spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      this->buildSolverMatrixWrapper(spSolver, matIdx, spEq, comp, idx);
   }

   void TimestepCoordinatorBase<true>::buildSolverMatrix(TimestepCoordinatorBase<true>::SharedZZSolverType spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      this->buildSolverMatrixWrapper(spSolver, matIdx, spEq, comp, idx);
   }

}
}
