/** 
 * @file SparseLinearCoordinatorImpl.cpp
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
#include "SparseSolvers/SparseLinearCoordinatorImpl.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"

namespace GeoMHDiSCC {

namespace Solver {

   SparseLinearCoordinatorImpl::SparseLinearCoordinatorImpl()
      : SparseLinearCoordinatorBase<SparseLinearSolver,Datatypes::ScalarSelector<Dimensions::Transform::TRA1D>::FIELD_IS_COMPLEX>()
   {
   }

   SparseLinearCoordinatorImpl::~SparseLinearCoordinatorImpl()
   {
   }

   void SparseLinearCoordinatorImpl::solve(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq, const ScalarVariable_map& scalVar, const VectorVariable_map& vectVar)
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

}
}
