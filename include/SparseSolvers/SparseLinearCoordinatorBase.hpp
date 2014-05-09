/** 
 * @file SparseLinearCoordinatorBase.hpp
 * @brief Implementation of the base for a general sparse linear solver coordinator
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SPARSELINEARCOORDINATORBASE_HPP
#define SPARSELINEARCOORDINATORBASE_HPP

// Debug includes
//
#include "Debug/DebuggerMacro.h"

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

// Project includes
//
#include "Base/MathConstants.hpp"
#include "SparseSolvers/SparseCoordinatorBase.hpp"
#include "Equations/IScalarEquation.hpp"
#include "Equations/IVectorEquation.hpp"

namespace GeoMHDiSCC {

namespace Solver {

   namespace internal {

      void addOperators(SparseMatrix& mat, const MHDFloat c, const DecoupledZSparse& decMat);

      void addOperators(SparseMatrixZ& mat, const MHDFloat c, const DecoupledZSparse& decMat);
   }

   /**
    * @brief Implementation of the base for a general sparse linear solver coordinator
    */
   template <typename TZSolver, typename TRSolver> class SparseLinearCoordinatorBase: public SparseCoordinatorBase<TZSolver,TRSolver>
   {
      public:
         /**
          * @brief Constructor
          */
         SparseLinearCoordinatorBase();

         /**
          * @brief Destructor
          */
         virtual ~SparseLinearCoordinatorBase();

         /**
          * @brief Initialise solver coordinator
          *
          * @param scalEq  Shared scalar equations
          * @param vectEq  Shared vector equations
          */
         void init(const typename SparseLinearCoordinatorBase<TZSolver,TRSolver>::ScalarEquation_range& scalEq, const typename SparseLinearCoordinatorBase<TZSolver,TRSolver>::VectorEquation_range& vectEq);
         
      protected:
         /**
          * @brief Build the solver matrices independently of solver type
          */
         template <typename TSolverIt> void buildSolverMatrices(Equations::SharedIEquation spEq, const SpectralFieldId id, const TSolverIt solveIt);

         /**
          * @brief Compute (coupled) matrices
          */
         void createMatrices(Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp);

         /**
          * @brief Solve all the linear systems
          */
         void solveSystems();

         /**
          * @brief Build the real solver matrix
          *
          * @param spSolver   Shared sparse real solver
          * @param matIdx     Index of the solver matrix
          * @param spEq    Shared pointer to equation
          * @param comp    Field component
          * @param idx     Matrix index
          */
         virtual void buildSolverMatrix(typename SparseLinearCoordinatorBase<TZSolver,TRSolver>::SharedRSolverType spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx) = 0;

         /**
          * @brief Build the complex solver matrix
          *
          * @param spSolver   Shared sparse real solver
          * @param matIdx     Index of the solver matrix
          * @param spEq    Shared pointer to equation
          * @param comp    Field component
          * @param idx     Matrix index
          */
         virtual void buildSolverMatrix(typename SparseLinearCoordinatorBase::SharedZSolverType spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx) = 0;

      private:
   };

   template <typename TZSolver, typename TRSolver> template <typename TSolverIt> void SparseLinearCoordinatorBase<TZSolver,TRSolver>::buildSolverMatrices(Equations::SharedIEquation spEq, const SpectralFieldId id, const TSolverIt solveIt)
   {
      // Number of linear systems
      int nSystems = spEq->couplingInfo(id.second).nSystems();

      // Loop over all substeps of solver
      for(this->mStep = 0; this->mStep < this->mNStep; this->mStep++)
      {
         // start index for matrices
         int start = this->mStep*nSystems;

         // Reserve storage for matrices and initialise vectors
         (*solveIt)->initMatrices(this->mNStep*nSystems);

         // Build the solver matrices
         for(int i = 0; i < nSystems; i++)
         {
            // Build LHS solver matrix
            this->buildSolverMatrix(*solveIt, start+i, spEq, id.second, i);
         }
      }
   }

   template <typename TZSolver, typename TRSolver> SparseLinearCoordinatorBase<TZSolver,TRSolver>::SparseLinearCoordinatorBase()
      : SparseCoordinatorBase<TZSolver,TRSolver>()
   {
   }

   template <typename TZSolver, typename TRSolver> SparseLinearCoordinatorBase<TZSolver,TRSolver>::~SparseLinearCoordinatorBase()
   {
   }

   template <typename TZSolver, typename TRSolver> void SparseLinearCoordinatorBase<TZSolver,TRSolver>::init(const typename SparseLinearCoordinatorBase<TZSolver,TRSolver>::ScalarEquation_range& scalEq, const typename SparseLinearCoordinatorBase<TZSolver,TRSolver>::VectorEquation_range& vectEq)
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
      // Loop over all substeps of solver
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
      // Create the problem matrices
      //

      DebuggerMacro_start("Linear: create operators", 2);
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
      DebuggerMacro_stop("Linear: create operators t = ", 2);

      //
      // Initialise the solvers and the initial state
      //

      DebuggerMacro_start("Linear: complex init", 2);
      // Initialise solvers from complex equation steppers
      typename SparseLinearCoordinatorBase<TZSolver,TRSolver>::SolverZ_iterator   solZIt;
      for(solZIt = this->mZSolvers.begin(); solZIt != this->mZSolvers.end(); ++solZIt)
      {
         DebuggerMacro_msg("---> complex solver", 2);

         (*solZIt)->initSolver();
      }
      DebuggerMacro_stop("Linear: complex init t = ", 2);

      DebuggerMacro_start("Linear: real init", 2);
      // Initialise solvers from real equation steppers
      typename SparseLinearCoordinatorBase<TZSolver,TRSolver>::SolverR_iterator   solRIt;
      for(solRIt = this->mRSolvers.begin(); solRIt != this->mRSolvers.end(); ++solRIt)
      {
         DebuggerMacro_msg("---> real solver", 2);

         (*solRIt)->initSolver();
      }
      DebuggerMacro_stop("Linear: real init t = ", 2);

      // Reset the step index
      this->mStep = 0;

      // Initialise with initial state
      this->initSolution(scalEq, vectEq);
   }

   template <typename TZSolver, typename TRSolver> void SparseLinearCoordinatorBase<TZSolver,TRSolver>::createMatrices(Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp)
   {
      // ID of the current field
      SpectralFieldId myId = std::make_pair(spEq->name(),comp);

      // Get solver index
      int myIdx = spEq->couplingInfo(myId.second).solverIndex();

      // Complex matrices in linear solve
      if(spEq->couplingInfo(myId.second).isComplex())
      {
         // Create iterator to current complex solver
         typename SparseLinearCoordinatorBase<TZSolver,TRSolver>::SolverZ_iterator solZIt = this->mZSolvers.begin();
         std::advance(solZIt, myIdx);

         // Build solver matrices
         if(! (*solZIt)->isInitialized())
         {
            this->buildSolverMatrices(spEq, myId, solZIt);
         }

      // Real matrices in linear solve
      } else
      {
         // Create iterator to current real solver
         typename SparseLinearCoordinatorBase<TZSolver,TRSolver>::SolverR_iterator solRIt = this->mRSolvers.begin();
         std::advance(solRIt, myIdx);

         // Build solver matrices
         if(! (*solRIt)->isInitialized())
         {
            this->buildSolverMatrices(spEq, myId, solRIt);
         }
      }
   }

   template <typename TZSolver, typename TRSolver> void SparseLinearCoordinatorBase<TZSolver,TRSolver>::solveSystems()
   {
      // Solve complex linear systems
      typename SparseLinearCoordinatorBase<TZSolver,TRSolver>::SolverZ_iterator   solZIt;
      for(solZIt = this->mZSolvers.begin(); solZIt != this->mZSolvers.end(); ++solZIt)
      {
         // Compute linear solve RHS
         (*solZIt)->solve(this->mStep);
      }

      // Solve real linear systems
      typename SparseLinearCoordinatorBase<TZSolver,TRSolver>::SolverR_iterator   solRIt;
      for(solRIt = this->mRSolvers.begin(); solRIt != this->mRSolvers.end(); ++solRIt)
      {
         // Compute linear solve RHS
         (*solRIt)->solve(this->mStep);
      }
   }

   namespace internal {

      inline void addOperators(SparseMatrix& mat, const MHDFloat c, const DecoupledZSparse& decMat)
      {
         assert(decMat.real().size() > 0);
         assert(decMat.imag().size() == 0 || decMat.imag().nonZeros() == 0);

         if(c != 1.0)
         {
            mat += c*decMat.real();
         } else
         {
            mat += decMat.real();
         }
      }

      inline void addOperators(SparseMatrixZ& mat, const MHDFloat c, const DecoupledZSparse& decMat)
      {
         assert(decMat.real().size() > 0);
         assert(decMat.imag().size() > 0);
         assert(decMat.real().size() == decMat.imag().size());

         if(c != 1.0)
         {
            mat += c*decMat.real().cast<MHDComplex>() + c*Math::cI*decMat.imag();
         } else
         {
            mat += decMat.real().cast<MHDComplex>() + Math::cI*decMat.imag();
         }
      }
   }
}
}

#endif // SPARSELINEARCOORDINATORBASE_HPP
