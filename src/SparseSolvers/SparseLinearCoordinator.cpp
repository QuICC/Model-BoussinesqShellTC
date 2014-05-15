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

   SparseLinearCoordinator<false>::SparseLinearCoordinator()
      : SparseLinearCoordinatorImpl()
   {
   }

   SparseLinearCoordinator<false>::~SparseLinearCoordinator()
   {
   }

   void SparseLinearCoordinator<false>::buildSolverMatrix(SparseLinearCoordinator<false>::SharedRRSolverType spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      this->buildSolverMatrixWrapper(spSolver, matIdx, spEq, comp, idx);
   }




   SparseLinearCoordinator<true>::SparseLinearCoordinator()
      : SparseLinearCoordinatorImpl()
   {
   }

   SparseLinearCoordinator<true>::~SparseLinearCoordinator()
   {
   }

   void SparseLinearCoordinator<true>::buildSolverMatrix(SparseLinearCoordinator<true>::SharedRZSolverType spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      this->buildSolverMatrixWrapper(spSolver, matIdx, spEq, comp, idx);
   }

   void SparseLinearCoordinator<true>::buildSolverMatrix(SparseLinearCoordinator<true>::SharedZZSolverType spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      this->buildSolverMatrixWrapper(spSolver, matIdx, spEq, comp, idx);
   }

}
}
