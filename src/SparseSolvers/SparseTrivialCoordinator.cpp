/** \file SparseTrivialCoordinator.cpp
 *  \brief Implementation of the for a general trivial solver structure
 */

// Debug includes
//
#include "Debug/DebuggerMacro.h"

// System includes
//

// External includes
//
#include <Eigen/Sparse>

// Class include
//
#include "SparseSolvers/SparseTrivialCoordinator.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"

namespace GeoMHDiSCC {

namespace Solver {

   SparseTrivialCoordinator::SparseTrivialCoordinator()
      : SparseCoordinatorBase<SharedSparseZTrivialSolver, SharedSparseDTrivialSolver>()
   {
   }

   SparseTrivialCoordinator::~SparseTrivialCoordinator()
   {
   }

   void SparseTrivialCoordinator::solve(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
   {
      // Update the equation input to the timestepper
      this->getInput(scalEq, vectEq);
      
      // Transfer timestep output back to equations
      this->transferOutput(scalEq, vectEq);

      // Update the internal step counter, counting from 0 to steps - 1
      this->mStep = (this->mStep + 1) % this->mNStep;
   }

   void SparseTrivialCoordinator::init(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
   {
      //
      // Create real/complex solvers
      //

      DebuggerMacro_start("Create trivial solvers", 0);
      // Loop over all scalar equations
      std::vector<Equations::SharedIScalarEquation>::const_iterator scalEqIt;
      for(scalEqIt = scalEq.first; scalEqIt < scalEq.second; scalEqIt++)
      {
         // Get type information for the solvers
         this->createSolver((*scalEqIt), FieldComponents::Spectral::SCALAR);
      }

      // Loop over all vector equations
      std::vector<Equations::SharedIVectorEquation>::const_iterator vectEqIt;
      for(vectEqIt = vectEq.first; vectEqIt < vectEq.second; vectEqIt++)
      {
         // Get type information for the solvers for the first component
         this->createSolver((*vectEqIt), FieldComponents::Spectral::ONE);

         // Get type information for the solvers for the first component
         this->createSolver((*vectEqIt), FieldComponents::Spectral::TWO);
      }
      DebuggerMacro_stop("Create trivial solvers t = ", 0);

      //
      // Create the timestep matrices
      //

      DebuggerMacro_start("Create solver storage", 0);
      // Loop over all substeps of timestepper
      for(this->mStep = 0; this->mStep < this->mNStep; this->mStep++)
      {
         // Loop over all scalar equations
         for(scalEqIt = scalEq.first; scalEqIt < scalEq.second; scalEqIt++)
         {
            // Create (coupled) matrices
            this->createStorage((*scalEqIt), FieldComponents::Spectral::SCALAR);
         }

         // Loop over all vector equations
         for(vectEqIt = vectEq.first; vectEqIt < vectEq.second; vectEqIt++)
         {
            // Create (coupled) matrices
            this->createStorage((*vectEqIt), FieldComponents::Spectral::ONE);

            // Create (coupled) matrices
            this->createStorage((*vectEqIt), FieldComponents::Spectral::TWO);
         }
      }
      DebuggerMacro_stop("Create solver storage t = ", 0);

      //
      // Initialise the solvers initial state
      //

      // Reset the step index
      this->mStep = 0;

      // Initialise with initial state
      this->initSolution(scalEq, vectEq);
   }

   void SparseTrivialCoordinator::addSolverD(const int start)
   {
      SharedSparseDTrivialSolver spSolver(new SparseDTrivialSolver(start));

      this->mDSolvers.push_back(spSolver);
   }

   void SparseTrivialCoordinator::addSolverZ(const int start)
   {
      SharedSparseZTrivialSolver spSolver(new SparseZTrivialSolver(start));

      this->mZSolvers.push_back(spSolver);
   }

}
}
