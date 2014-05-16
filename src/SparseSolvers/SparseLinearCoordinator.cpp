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

      // Update the internal step counter, counting from 0 to steps - 1
      this->mStep = (this->mStep + 1) % this->mNStep;
   }

   void SparseLinearCoordinator::buildSolverMatrix(SparseLinearCoordinator::SharedRRSolverType spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      this->buildSolverMatrixWrapper(spSolver, matIdx, spEq, comp, idx);
   }

   void SparseLinearCoordinator::buildSolverMatrix(SparseLinearCoordinator::SharedRZSolverType spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      this->buildSolverMatrixWrapper(spSolver, matIdx, spEq, comp, idx);
   }

   void SparseLinearCoordinator::buildSolverMatrix(SparseLinearCoordinator::SharedZZSolverType spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      this->buildSolverMatrixWrapper(spSolver, matIdx, spEq, comp, idx);
   }

}
}
