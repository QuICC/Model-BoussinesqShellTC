/** 
 * @file SparseLinearCoordinatorBase.cpp
 * @brief Implementation of the for a general linear solver structure
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
#include "SparseSolvers/SparseLinearCoordinatorBase.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"

namespace GeoMHDiSCC {

namespace Solver {

   SparseLinearCoordinatorBase::SparseLinearCoordinatorBase()
      : SparseCoordinatorBase<SharedSparseZLinearSolver,SharedSparseDLinearSolver>()
   {
   }

   SparseLinearCoordinatorBase::~SparseLinearCoordinatorBase()
   {
   }

   void SparseLinearCoordinatorBase::init(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
   {
      //
      // Create real/complex solvers
      //

      DebuggerMacro_start("Linear: create solvers", 2);
      // Loop over all scalar equations
      std::vector<Equations::SharedIScalarEquation>::const_iterator scalEqIt;
      for(scalEqIt = scalEq.first; scalEqIt < scalEq.second; scalEqIt++)
      {
         DebuggerMacro_msg("---> scalar solver", 2);

         // Get type information for the linear solvers
         this->createSolver((*scalEqIt), FieldComponents::Spectral::SCALAR);
      }

      // Loop over all vector equations
      std::vector<Equations::SharedIVectorEquation>::const_iterator vectEqIt;
      for(vectEqIt = vectEq.first; vectEqIt < vectEq.second; vectEqIt++)
      {
         DebuggerMacro_msg("---> vector solver", 2);

         // Get type information for the linear solvers 
         Equations::IVectorEquation::SpectralComponent_iterator compIt;
         Equations::IVectorEquation::SpectralComponent_range  compRange = (*vectEqIt)->spectralRange();
         for(compIt = compRange.first; compIt != compRange.second; ++compIt)
         {
            this->createSolver((*vectEqIt), *compIt);
         }
      }
      DebuggerMacro_stop("Linear: create solvers t = ", 2);

      //
      // Initialise the solver storage
      //

      DebuggerMacro_start("Linear: create storage", 2);
      // Loop over all substeps of timestepper
      for(this->mStep = 0; this->mStep < this->mNStep; this->mStep++)
      {
         // Loop over all scalar equations
         for(scalEqIt = scalEq.first; scalEqIt < scalEq.second; scalEqIt++)
         {
            DebuggerMacro_msg("---> scalar storage", 2);

            // Create storage
            this->createStorage((*scalEqIt), FieldComponents::Spectral::SCALAR);
         }

         // Loop over all vector equations
         for(vectEqIt = vectEq.first; vectEqIt < vectEq.second; vectEqIt++)
         {
            DebuggerMacro_msg("---> vector storage", 2);

            // Create storage
            Equations::IVectorEquation::SpectralComponent_iterator compIt;
            Equations::IVectorEquation::SpectralComponent_range  compRange = (*vectEqIt)->spectralRange();
            for(compIt = compRange.first; compIt != compRange.second; ++compIt)
            {
               this->createStorage((*vectEqIt), *compIt);
            }
         }
      }
      DebuggerMacro_stop("Linear: create storage t = ", 2);

      //
      // Create the timestep matrices
      //

      DebuggerMacro_start("Linear: create operators", 2);
      // Loop over all substeps of timestepper
      for(this->mStep = 0; this->mStep < this->mNStep; this->mStep++)
      {
         // Loop over all scalar equations
         for(scalEqIt = scalEq.first; scalEqIt < scalEq.second; scalEqIt++)
         {
            DebuggerMacro_msg("---> scalar operators", 2);

            // Create (coupled) matrices
            this->createMatrices((*scalEqIt), FieldComponents::Spectral::SCALAR);
         }

         // Loop over all vector equations
         for(vectEqIt = vectEq.first; vectEqIt < vectEq.second; vectEqIt++)
         {
            DebuggerMacro_msg("---> vector operators", 2);

            // Create (coupled) matrices
            Equations::IVectorEquation::SpectralComponent_iterator compIt;
            Equations::IVectorEquation::SpectralComponent_range  compRange = (*vectEqIt)->spectralRange();
            for(compIt = compRange.first; compIt != compRange.second; ++compIt)
            {
               this->createMatrices((*vectEqIt), *compIt);
            }
         }
      }
      DebuggerMacro_stop("Linear: create operators t = ", 2);

      //
      // Initialise the solvers and the initial state
      //

      DebuggerMacro_start("Linear: complex init", 2);
      // Initialise solvers from complex equation steppers
      SolverZ_iterator   solZIt;
      for(solZIt = this->mZSolvers.begin(); solZIt != this->mZSolvers.end(); ++solZIt)
      {
         DebuggerMacro_msg("---> complex solver", 2);

         (*solZIt)->initSolver();
      }
      DebuggerMacro_stop("Linear: complex init t = ", 2);

      DebuggerMacro_start("Linear: real init", 2);
      // Initialise solvers from real equation steppers
      SolverD_iterator   solDIt;
      for(solDIt = this->mDSolvers.begin(); solDIt != this->mDSolvers.end(); ++solDIt)
      {
         DebuggerMacro_msg("---> real solver", 2);

         (*solDIt)->initSolver();
      }
      DebuggerMacro_stop("Linear: real init t = ", 2);

      // Reset the step index
      this->mStep = 0;

      // Initialise with initial state
      this->initSolution(scalEq, vectEq);
   }

   void SparseLinearCoordinatorBase::createMatrices(Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp)
   {
      // ID of the current field
      SpectralFieldId myId = std::make_pair(spEq->name(),comp);

      // Index of the current field
      int myIdx = spEq->couplingInfo(myId.second).solverIndex();

      // Complex matrices in linear solve
      if(spEq->couplingInfo(myId.second).isComplex())
      {
         // Create iterator to current complex solver
         SolverZ_iterator solZIt = this->mZSolvers.begin();
         std::advance(solZIt, myIdx);

         // Build solver matrices
         this->buildSolverMatrices(spEq, myId, solZIt);

      // Real matrices in linear solve
      } else
      {
         // Create iterator to current real solver
         SolverD_iterator solDIt = this->mDSolvers.begin();
         std::advance(solDIt, myIdx);

         // Build solver matrices
         this->buildSolverMatrices(spEq, myId, solDIt);
      }
   }

   void SparseLinearCoordinatorBase::solveSystems()
   {
      // Solve complex linear systems
      SolverZ_iterator   solZIt;
      for(solZIt = this->mZSolvers.begin(); solZIt != this->mZSolvers.end(); ++solZIt)
      {
         // Compute linear solve RHS
         (*solZIt)->solve(this->mStep);
      }

      // Solve real linear systems
      SolverD_iterator   solDIt;
      for(solDIt = this->mDSolvers.begin(); solDIt != this->mDSolvers.end(); ++solDIt)
      {
         // Compute linear solve RHS
         (*solDIt)->solve(this->mStep);
      }
   }

}
}
