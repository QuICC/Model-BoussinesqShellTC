/** 
 * @file SparseLinearCoordinator.cpp
 * @brief Implementation of a general linear solver structure
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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
#include "SparseSolvers/SparseLinearCoordinator.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"

namespace GeoMHDiSCC {

namespace Solver {

   SparseLinearCoordinator::SparseLinearCoordinator()
      : SparseLinearCoordinatorBase<SparseLinearSolver>()
   {
   }

   SparseLinearCoordinator::~SparseLinearCoordinator()
   {
   }

   void SparseLinearCoordinator::solve(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq, const ScalarVariable_map& scalVar, const VectorVariable_map& vectVar)
   {
      // Update the equation input to the timestepper
      this->getInput(scalEq, vectEq, scalVar, vectVar);

      // Solve all the linear systems
      this->solveSystems();
      
      // Transfer timestep output back to equations
      this->transferOutput(scalEq, vectEq);

      // Clear the solver RHS
      this->clearSolvers();
   }

   void SparseLinearCoordinator::buildSolverMatrix(SparseLinearCoordinator::SharedRealSolverType spSolver, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      buildLinearSolverMatrixWrapper(spSolver, spEq, comp, idx);
   }

   void SparseLinearCoordinator::buildSolverMatrix(SparseLinearCoordinator::SharedComplexSolverType spSolver, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      buildLinearSolverMatrixWrapper(spSolver, spEq, comp, idx);
   }

}
}
