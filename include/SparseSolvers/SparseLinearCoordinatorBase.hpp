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
#include "TypeSelectors/ScalarSelector.hpp"
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
   template <template <class,class> class TSolver> class SparseLinearCoordinatorBase: public SparseCoordinatorBase<TSolver>
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
         void init(const typename SparseLinearCoordinatorBase<TSolver>::ScalarEquation_range& scalEq, const typename SparseLinearCoordinatorBase<TSolver>::VectorEquation_range& vectEq);

         /**
          * @brief Build the solver matrices independently of solver type
          *
          * \mhdBug Should not be public
          */
         template <typename TSolverIt> void buildSolverMatrices(Equations::SharedIEquation spEq, const SpectralFieldId id, const TSolverIt solveIt);
         
      protected:

         /**
          * @brief Compute (coupled) matrices
          */
         void createMatrices(Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp);

         /**
          * @brief Solve all the linear systems
          */
         void solveSystems();

         /**
          * @brief Build the real operator
          *
          * @param spSolver   Shared sparse real solver
          * @param matIdx     Index of the solver matrix
          * @param spEq    Shared pointer to equation
          * @param comp    Field component
          * @param idx     Matrix index
          */
         virtual void buildSolverMatrix(typename SparseLinearCoordinatorBase<TSolver>::SharedRealSolverType spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx) = 0;

         /**
          * @brief Build the complex operator
          *
          * @param spSolver   Shared sparse real solver
          * @param matIdx     Index of the solver matrix
          * @param spEq    Shared pointer to equation
          * @param comp    Field component
          * @param idx     Matrix index
          */
         virtual void buildSolverMatrix(typename SparseLinearCoordinatorBase::SharedComplexSolverType spSolver, const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx) = 0;

      private:
   };

   template <template <class,class> class TSolver, typename TSolverIt> void buildSolverMatricesSolver(SparseLinearCoordinatorBase<TSolver>& coord, Equations::SharedIEquation spEq, const int idx, SpectralFieldId id);

   //
   //
   //

   template <template <class,class> class TSolver> template <typename TSolverIt> void SparseLinearCoordinatorBase<TSolver>::buildSolverMatrices(Equations::SharedIEquation spEq, const SpectralFieldId id, const TSolverIt solveIt)
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

   template <template <class,class> class TSolver> SparseLinearCoordinatorBase<TSolver>::SparseLinearCoordinatorBase()
      : SparseCoordinatorBase<TSolver>()
   {
   }

   template <template <class,class> class TSolver> SparseLinearCoordinatorBase<TSolver>::~SparseLinearCoordinatorBase()
   {
   }

   template <template <class,class> class TSolver> void SparseLinearCoordinatorBase<TSolver>::init(const typename SparseLinearCoordinatorBase<TSolver>::ScalarEquation_range& scalEq, const typename SparseLinearCoordinatorBase<TSolver>::VectorEquation_range& vectEq)
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

      // Loop over all substeps of solver
      DebuggerMacro_start("Linear: create storage", 2);
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

      // Initialise solvers from complex equation steppers
      DebuggerMacro_start("Linear: complex operator init", 2);
      initSolvers<TSolver,typename SparseCoordinatorBase<TSolver>::ComplexSolver_iterator>(*this);
      DebuggerMacro_stop("Linear: complex operator init t = ", 2);

      // Initialise solvers from real equation steppers
      DebuggerMacro_start("Linear: real operator init", 2);
      initSolvers<TSolver,typename SparseCoordinatorBase<TSolver>::RealSolver_iterator>(*this);
      DebuggerMacro_stop("Linear: real operator init t = ", 2);

      // Reset the step index
      this->mStep = 0;

      // Initialise with initial state
      this->initSolution(scalEq, vectEq);
   }

   template <template <class,class> class TSolver> void SparseLinearCoordinatorBase<TSolver>::createMatrices(Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp)
   {
      // ID of the current field
      SpectralFieldId myId = std::make_pair(spEq->name(),comp);

      // Get solver index
      int myIdx = spEq->couplingInfo(myId.second).solverIndex();

      // System operator is complex
      if(spEq->couplingInfo(myId.second).isComplex())
      {
         buildSolverMatricesSolver<TSolver,typename SparseCoordinatorBase<TSolver>::ComplexSolver_iterator>(*this, spEq, myIdx, myId);

      // System operator is real
      } else
      {
         buildSolverMatricesSolver<TSolver,typename SparseCoordinatorBase<TSolver>::RealSolver_iterator>(*this, spEq, myIdx, myId);
      }
   }

   template <template <class,class> class TSolver> void SparseLinearCoordinatorBase<TSolver>::solveSystems()
   {
      // Solve complex operator, complex field linear systems
      solveSolvers<TSolver,typename SparseCoordinatorBase<TSolver>::ComplexSolver_iterator>(*this, this->mStep);

      // Solve real operator, complex field linear systems
      solveSolvers<TSolver,typename SparseCoordinatorBase<TSolver>::RealSolver_iterator>(*this, this->mStep);
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

   //
   //
   //

   template <template <class,class> class TSolver, typename TSolverIt> void buildSolverMatricesSolver(SparseLinearCoordinatorBase<TSolver>& coord, Equations::SharedIEquation spEq, const int idx, SpectralFieldId id)
   {
      // Create iterator to current solver
      TSolverIt solIt;
      coord.setIterator(solIt);
      std::advance(solIt, idx);

      // Build solver matrices
      if(! (*solIt)->isInitialized())
      {
         coord.buildSolverMatrices(spEq, id, solIt);
      }
   }
}
}

#endif // SPARSELINEARCOORDINATORBASE_HPP
