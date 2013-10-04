/** 
 * @file SparseTrivialCoordinator.cpp
 * @brief Implementation of the for a general trivial solver structure
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

// Class include
//
#include "SparseSolvers/SparseTrivialCoordinator.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"

namespace GeoMHDiSCC {

namespace Solver {

   SparseTrivialCoordinator::SparseTrivialCoordinator()
      : SparseCoordinatorBase<SparseZTrivialSolver, SparseRZTrivialSolver>()
   {
   }

   SparseTrivialCoordinator::~SparseTrivialCoordinator()
   {
   }

   void SparseTrivialCoordinator::solve(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq, const ScalarVariable_map& scalVar, const VectorVariable_map& vectVar)
   {
      // Update the equation input to the timestepper
      this->getInput(scalEq, vectEq, scalVar, vectVar);
      
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

      DebuggerMacro_start("Trivial: create solvers", 2);
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
         // Get type information for the solvers
         Equations::IVectorEquation::SpectralComponent_iterator compIt;
         Equations::IVectorEquation::SpectralComponent_range  compRange = (*vectEqIt)->spectralRange();
         for(compIt = compRange.first; compIt != compRange.second; ++compIt)
         {
            this->createSolver((*vectEqIt), *compIt);
         }
      }
      DebuggerMacro_stop("Trivial: create solvers t = ", 2);

      //
      // Initialise the solver storage
      //

      DebuggerMacro_start("Trivial: create storage", 2);
      // Loop over all substeps of timestepper
      for(this->mStep = 0; this->mStep < this->mNStep; this->mStep++)
      {
         // Loop over all scalar equations
         for(scalEqIt = scalEq.first; scalEqIt < scalEq.second; scalEqIt++)
         {
            // Create storage 
            this->createStorage((*scalEqIt), FieldComponents::Spectral::SCALAR);
         }

         // Loop over all vector equations
         for(vectEqIt = vectEq.first; vectEqIt < vectEq.second; vectEqIt++)
         {
            // Create storage 
            Equations::IVectorEquation::SpectralComponent_iterator compIt;
            Equations::IVectorEquation::SpectralComponent_range  compRange = (*vectEqIt)->spectralRange();
            for(compIt = compRange.first; compIt != compRange.second; ++compIt)
            {
               this->createStorage((*vectEqIt), *compIt);
            }
         }
      }
      DebuggerMacro_stop("Trivial: create storage t = ", 2);

      //
      // Initialise the solvers initial state
      //

      // Reset the step index
      this->mStep = 0;

      // Initialise with initial state
      this->initSolution(scalEq, vectEq);
   }

}
}
