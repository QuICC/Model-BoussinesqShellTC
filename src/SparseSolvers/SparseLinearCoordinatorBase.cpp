/** \file SparseLinearCoordinatorBase.cpp
 *  \brief Implementation of the for a general linear solver structure
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
      : mStep(0)
   {
   }

   SparseLinearCoordinatorBase::~SparseLinearCoordinatorBase()
   {
   }

   bool SparseLinearCoordinatorBase::finishedStep() const
   {
      return this->mStep == 0;
   }

   void SparseLinearCoordinatorBase::init(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
   {
      //
      // Create real/complex solvers
      //

      DebuggerMacro_start("Create linear solvers", 0);
      // Loop over all scalar equations
      std::vector<Equations::SharedIScalarEquation>::const_iterator scalEqIt;
      for(scalEqIt = scalEq.first; scalEqIt < scalEq.second; scalEqIt++)
      {
         // Get type information for the linear solvers
         this->createsolver((*scalEqIt), FieldComponents::Spectral::SCALAR);
      }

      // Loop over all vector equations
      std::vector<Equations::SharedIVectorEquation>::const_iterator vectEqIt;
      for(vectEqIt = vectEq.first; vectEqIt < vectEq.second; vectEqIt++)
      {
         // Get type information for the linear solvers for the first component
         this->createSolver((*vectEqIt), FieldComponents::Spectral::ONE);

         // Get type information for the linear solvers for the first component
         this->createSolver((*vectEqIt), FieldComponents::Spectral::TWO);
      }
      DebuggerMacro_stop("Create linear solvers t = ", 0);

      //
      // Create the timestep matrices
      //

      DebuggerMacro_start("Create linear operators", 0);
      // Loop over all substeps of timestepper
AARRRRRGHHHHHH
      for(this->mStep = 0; this->mStep < ImExRK3::STEPS; this->mStep++)
      {
         // Loop over all scalar equations
         for(scalEqIt = scalEq.first; scalEqIt < scalEq.second; scalEqIt++)
         {
            // Create (coupled) matrices
            this->createMatrices((*scalEqIt), FieldComponents::Spectral::SCALAR);
         }

         // Loop over all vector equations
         for(vectEqIt = vectEq.first; vectEqIt < vectEq.second; vectEqIt++)
         {
            // Create (coupled) matrices
            this->createMatrices((*vectEqIt), FieldComponents::Spectral::ONE);

            // Create (coupled) matrices
            this->createMatrices((*vectEqIt), FieldComponents::Spectral::TWO);
         }
      }
      DebuggerMacro_stop("Create linear operators t = ", 0);

      //
      // Initialise the solvers and the initial state
      //

      DebuggerMacro_start("Complex solver init", 0);
      // Initialise solvers from complex equation steppers
      SolverZ_iterator   solZIt;
      for(solZIt = this->mZSolvers.begin(); solZIt != this->mZSolvers.end(); ++solZIt)
      {
         solZIt->initSolver();
      }
      DebuggerMacro_stop("Complex solver init t = ", 0);

      DebuggerMacro_start("Real solver init", 0);
      // Initialise solvers from real equation steppers
      SolverD_iterator   solDIt;
      for(solDIt = this->mDSolvers.begin(); solDIt != this->mDSolvers.end(); ++solDIt)
      {
         solDIt->initSolver();
      }
      DebuggerMacro_stop("Real solver init t = ", 0);

      // Reset the step index
      this->mStep = 0;

      // Initialise with initial state
      this->initSolution(scalEq, vectEq);
   }

   void SparseLinearCoordinatorBase::createSolver(Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp)
   {
      // Equation is part of a complex system
      if(spEq->couplingInfo(comp).isComplex())
      {
         // Add linear solver if system index does not yet exist
         if(spEq->couplingInfo(comp).solverIndex() > static_cast<int>(this->mZSolvers.size()) - 1)
         {
            this->addSolverZ(spEq->couplingInfo(comp).fieldStart());
         }

      // Equation is part of a real system
      } else
      {
         // Add equation stepper if system index does not yet exist
         if(spEq->couplingInfo(comp).solverIndex() > static_cast<int>(this->mDSolvers.size()) - 1)
         {
            this->addSolverD(spEq->couplingInfo(comp).fieldStart());
         }
      }
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
         this->buildSolverMatrix(spEq, myId, solZIt);

      // Real matrices in linear solve
      } else
      {
         // Create iterator to current real solver
         SolverD_iterator solDIt = this->mDSolvers.begin();
         std::advance(solDIt, myIdx);

         // Build solver matrices
         this->buildSolverMatrix(spEq, myId, solDIt);
      }
   }

   void SparseLinearCoordinatorBase::initSolution(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
   {
      // Storage for information and identity
      SpectralFieldId myId;

      // Loop over all scalar equations
      std::vector<Equations::SharedIScalarEquation>::const_iterator scalEqIt;
      for(scalEqIt = scalEq.first; scalEqIt < scalEq.second; scalEqIt++)
      {
         // Get field identity
         myId = std::make_pair((*scalEqIt)->name(), FieldComponents::Spectral::SCALAR);

         // Get index of solver
         int myIdx = (*scalEqIt)->couplingInfo(myId.second).solverIndex();

         // Linear solve matrices are complex
         if((*scalEqIt)->couplingInfo(myId.second).isComplex())
         {
            // Create iterator to current complex solver
            SolverZ_iterator solZIt = this->mZSolvers.begin();
            std::advance(solZIt, myIdx);

            // Get solver input
            for(int i = 0; i < solZIt->nSystem(); i++)
            {
               Equations::copyUnknown(*(*scalEqIt), myId.second, solZIt->rSolution(i), i, solZIt->startRow(myId,i));
            }

         // Linear solve matrices are real
         } else
         {
            // Create iterator to current real solver
            SolverD_iterator solDIt = this->mDSolvers.begin();
            std::advance(solDIt, myIdx);

            // Get solver input
            for(int i = 0; i < solDIt->nSystem(); i++)
            {
               Equations::copyUnknown(*(*scalEqIt), myId.second, solDIt->rSolution(i), i, solDIt->startRow(myId,i));
            }
         }
      }

      // Loop over all vector equations
      std::vector<Equations::SharedIVectorEquation>::const_iterator vectEqIt;
      for(vectEqIt = vectEq.first; vectEqIt < vectEq.second; vectEqIt++)
      {
         // Get field identity
         myId = std::make_pair((*vectEqIt)->name(), FieldComponents::Spectral::ONE);

         // Get index of solver
         int myIdx = (*vectEqIt)->couplingInfo(myId.second).solverIndex();

         // Linear solve matrices are complex
         if((*vectEqIt)->couplingInfo(myId.second).isComplex())
         {
            // Create iterator to current complex solver
            SolverZ_iterator solZIt = this->mZSolvers.begin();
            std::advance(solZIt, myIdx);

            // Get solver input for toroidal component
            for(int i = 0; i < solZIt->nSystem(); i++)
            {
               Equations::copyUnknown(*(*vectEqIt), myId.second, solZIt->rSolution(i), i, solZIt->startRow(myId,i));
            }

         // Linear solve matrices are real
         } else
         {
            // Create iterator to current real solver
            SolverD_iterator solDIt = this->mDSolvers.begin();
            std::advance(solDIt, myIdx);

            // Get solver input for toroidal component
            for(int i = 0; i < solDIt->nSystem(); i++)
            {
               Equations::copyUnknown(*(*vectEqIt), myId.second, solDIt->rSolution(i), i, solDIt->startRow(myId,i));
            }
         }

         // Get field identity
         myId = std::make_pair((*vectEqIt)->name(), FieldComponents::Spectral::TWO);

         // Get index of solver
         myIdx = (*vectEqIt)->couplingInfo(myId.second).solverIndex();

         // Linear solve matrices are complex
         if((*vectEqIt)->couplingInfo(myId.second).isComplex())
         {
            // Create iterator to current complex solver
            SolverZ_iterator solZIt = this->mZSolvers.begin();
            std::advance(solZIt, myIdx);

            // Get solver input for poloidal component
            for(int i = 0; i < solZIt->nSystem(); i++)
            {
               Equations::copyUnknown(*(*vectEqIt), myId.second, solZIt->rSolution(i), i, solZIt->startRow(myId,i));
            }

         // Linear solve matrices are real
         } else
         {
            // Create iterator to current real solver
            SolverD_iterator solDIt = this->mDSolvers.begin();
            std::advance(solDIt, myIdx);

            // Get solver input for poloidal component
            for(int i = 0; i < solDIt->nSystem(); i++)
            {
               Equations::copyUnknown(*(*vectEqIt), myId.second, solDIt->rSolution(i), i, solDIt->startRow(myId,i));
            }
         }
      }
   }

   void SparseLinearCoordinatorBase::getInput(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
   {
      // Storage for information and identity
      SpectralFieldId myId;

      // Loop over all scalar equations
      std::vector<Equations::SharedIScalarEquation>::const_iterator scalEqIt;
      for(scalEqIt = scalEq.first; scalEqIt < scalEq.second; scalEqIt++)
      {
         // Get field identity
         myId = std::make_pair((*scalEqIt)->name(), FieldComponents::Spectral::SCALAR);

         // Get index of solver
         int myIdx = (*scalEqIt)->couplingInfo(myId.second).solverIndex();

         // Linear solve matrices are complex
         if((*scalEqIt)->couplingInfo(myId.second).isComplex())
         {
            // Create iterator to current complex solver
            SolverZ_iterator solZIt = this->mZSolvers.begin();
            std::advance(solZIt, myIdx);

            // Get solver input
            this->getSolverInput(scalEqIt, myId, solZIt);

         // Linear solve matrices are real
         } else
         {
            // Create iterator to current real solver
            SolverD_iterator solDIt = this->mDSolvers.begin();
            std::advance(solDIt, myIdx);

            // Get solver input
            this->getSolverInput(scalEqIt, myId, solDIt);
         }
      }

      // Loop over all vector equations
      std::vector<Equations::SharedIVectorEquation>::const_iterator vectEqIt;
      for(vectEqIt = vectEq.first; vectEqIt < vectEq.second; vectEqIt++)
      {
         // Get field identity for first component
         myId = std::make_pair((*vectEqIt)->name(), FieldComponents::Spectral::ONE);

         // Get index of solver
         int myIdx = (*vectEqIt)->couplingInfo(myId.second).solverIndex();

         // Linear solve matrices are complex
         if((*vectEqIt)->couplingInfo(myId.second).isComplex())
         {
            // Create iterator to current complex solver
            SolverZ_iterator solZIt = this->mZSolvers.begin();
            std::advance(solZIt, myIdx);

            // Get solver input
            this->getSolverInput(vectEqIt, myId, solZIt);

         // Linear solve matrices are real
         } else
         {
            // Create iterator to current real solver
            SolverD_iterator solDIt = this->mDSolvers.begin();
            std::advance(solDIt, myIdx);

            // Get solver input
            this->getSolverInput(vectEqIt, myId, solDIt);
         }

         // Get field identity for second component
         myId = std::make_pair((*vectEqIt)->name(), FieldComponents::Spectral::TWO);

         // Get index of solver
         myIdx = (*vectEqIt)->couplingInfo(myId.second).solverIndex();

         // Linear solve matrices are complex
         if((*vectEqIt)->couplingInfo(myId.second).isComplex())
         {
            // Create iterator to current complex solver
            SolverZ_iterator solZIt = this->mZSolvers.begin();
            std::advance(solZIt, myIdx);

            // Get solver input
            this->getSolverInput(vectEqIt, myId, solZIt);

         // Linear solve matrices are real
         } else
         {
            // Create iterator to current real solver
            SolverD_iterator solDIt = this->mDSolvers.begin();
            std::advance(solDIt, myIdx);

            // Get solver input
            this->getSolverInput(vectEqIt, myId, solDIt);
         }
      }
   }

   void SparseLinearCoordinatorBase::solveSystems()
   {
      // Solve complex linear systems
      SolverZ_iterator   solZIt;
      for(soZIt = this->mZSolvers.begin(); solZIt != this->mZSolvers.end(); ++solZIt)
      {
         // Compute linear solve RHS
         solZIt->solve(this->mStep);
      }

      // Solve real linear systems
      SolverD_iterator   solDIt;
      for(solDIt = this->mDSolvers.begin(); solDIt != this->mDSolvers.end(); ++solDIt)
      {
         // Compute linear solve RHS
         solDIt->solve(this->mStep);
      }
   }

   void SparseLinearCoordinatorBase::transferOutput(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
   {
      // Storage for identity
      SpectralFieldId myId;

      // Loop over all scalar equations
      std::vector<Equations::SharedIScalarEquation>::const_iterator scalEqIt;
      for(scalEqIt = scalEq.first; scalEqIt < scalEq.second; scalEqIt++)
      {
         // Get field identity
         myId = std::make_pair((*scalEqIt)->name(), FieldComponents::Spectral::SCALAR);

         // Get index of solver
         int myIdx = (*scalEqIt)->couplingInfo(myId.second).solverIndex();

         // Linear solve matrices are complex
         if((*scalEqIt)->couplingInfo(myId.second).isComplex())
         {
            // Create iterator to current complex solver
            SolverZ_iterator solZIt = this->mZSolvers.begin();
            std::advance(solZIt, myIdx);

            // Get solver output
            for(int i = 0; i < solZIt->nSystem(); i++)
            {
               (*scalEqIt)->storeSolution(myId.second, solZIt->solution(i), i, solZIt->startRow(myId,i));
            }

         // Linear solve matrices are real
         } else
         {
            // Create iterator to current real solver
            SolverD_iterator solDIt = this->mDSolvers.begin();
            std::advance(solDIt, myIdx);

            // Get solver output
            for(int i = 0; i < solDIt->nSystem(); i++)
            {
               (*scalEqIt)->storeSolution(myId.second, solDIt->solution(i), i, solDIt->startRow(myId,i));
            }
         }
      }

      // Loop over all vector equations
      std::vector<Equations::SharedIVectorEquation>::const_iterator vectEqIt;
      for(vectEqIt = vectEq.first; vectEqIt < vectEq.second; vectEqIt++)
      {
         // Get field identity
         myId = std::make_pair((*vectEqIt)->name(), FieldComponents::Spectral::ONE);

         // Get index of solver
         int myIdx = (*vectEqIt)->couplingInfo(myId.second).solverIndex();

         // Linear solve matrices are complex
         if((*vectEqIt)->couplingInfo(myId.second).isComplex())
         {
            // Create iterator to current complex solver
            SolverZ_iterator solZIt = this->mZSolvers.begin();
            std::advance(solZIt, myIdx);

            // Get solver output for first component
            for(int i = 0; i < solZIt->nSystem(); i++)
            {
               (*vectEqIt)->storeSolution(myId.second, solZIt->solution(i), i, solZIt->startRow(myId,i));
            }

         // Linear solve matrices are real
         } else
         {
            // Create iterator to current real solver
            SolverD_iterator solDIt = this->mDSolvers.begin();
            std::advance(solDIt, myIdx);

            // Get solver output for first component
            for(int i = 0; i < solDIt->nSystem(); i++)
            {
               (*vectEqIt)->storeSolution(myId.second, solDIt->solution(i), i, solDIt->startRow(myId,i));
            }
         }

         // Get field identity
         myId = std::make_pair((*vectEqIt)->name(), FieldComponents::Spectral::TWO);

         // Get index of solver
         myIdx = (*vectEqIt)->couplingInfo(myId.second).solverIndex();

         // Linear solve matrices are complex
         if((*vectEqIt)->couplingInfo(myId.second).isComplex())
         {
            // Create iterator to current complex solver
            SolverZ_iterator solZIt = this->mZSolvers.begin();
            std::advance(solZIt, myIdx);

            // Get solver input for second component
            for(int i = 0; i < solZIt->nSystem(); i++)
            {
               (*vectEqIt)->storeSolution(myId.second, solZIt->solution(i), i, solZIt->startRow(myId,i));
            }

         // Linear solve matrices are real
         } else
         {
            // Create iterator to current real solver
            SolverD_iterator solDIt = this->mDSolvers.begin();
            std::advance(solDIt, myIdx);

            // Get solver output for second component
            for(int i = 0; i < solDIt->nSystem(); i++)
            {
               (*vectEqIt)->storeSolution(myId.second, solDIt->solution(i), i, solDIt->startRow(myId,i));
            }
         }
      }
   }

}
}
