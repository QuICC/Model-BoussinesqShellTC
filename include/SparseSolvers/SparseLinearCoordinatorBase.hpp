/** \file SparseLinearCoordinatorBase.hpp
 *  \brief Implementation of the base for a general sparse linear solver coordinator
 */

#ifndef SPARSELINEARCOORDINATORBASE_HPP
#define SPARSELINEARCOORDINATORBASE_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "SparseLinearSolver/SparseDLinearSolver.hpp"
#include "SparseLinearSolver/SparseZLinearSolver.hpp"
#include "Equations/IScalarEquation.hpp"
#include "Equations/IVectorEquation.hpp"

namespace GeoMHDiSCC {

namespace Solver {

   /**
    * @brief Implementation of the base for a general sparse linear solver coordinatore
    */
   class SparseLinearCoordinatorBase
   {
      public:
         /// Typedef for an iterator to a complex linear solver
         typedef std::vector<SharedSparseZLinearSolver>::iterator   SolverZ_iterator;

         /// Typedef for an iterator to a real linear solver
         typedef std::vector<SharedSparseDLinearSolver>::iterator   SolverD_iterator;

         /// Typedef for a shared scalar equation iterator
         typedef std::vector<Equations::SharedIScalarEquation>::iterator   ScalarEquation_iterator;

         /// Typedef for a shared vector equation iterator
         typedef std::vector<Equations::SharedIVectorEquation>::iterator   VectorEquation_iterator;

         /// Typedef for a shared scalar equation range
         typedef std::pair<ScalarEquation_iterator, ScalarEquation_iterator>  ScalarEquation_range;

         /// Typedef for a shared vector equation range
         typedef std::pair<VectorEquation_iterator, VectorEquation_iterator>  VectorEquation_range;

         /**
          * @brief Constructor
          */
         SparseLinearCoordinatorBase();

         /**
          * @brief Destructor
          */
         ~SparseLinearCoordinatorBase();

         /**
          * @brief Finished computation of solver step?
          */
         bool finishedStep() const;

         /**
          * @brief Initialise solver coordinator
          *
          * @param scalEq  Shared scalar equations
          * @param vectEq  Shared vector equations
          */
         void init(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq);
         
      protected:
         /**
          * @brief Create a real linear solver
          */
         virtual void addSolverD(const int start) = 0;

         /**
          * @brief Create a complex linear solver
          */
         virtual void addSolverZ(const int start) = 0;

         /**
          * @brief Get the solver input independently of solver type
          */
         template <typename TEquationIt, typename TSolverIt> void getSolverInput(const TEquationIt eqIt, const SpectralFieldId id, const TSolverIt solveIt);

         /**
          * @brief Build the solver matrices independently of solver type
          */
         template <typename TSolverIt> void buildSolverMatrices(Equations::SharedIEquation spEq, const SpectralFieldId id, const TSolverIt solveIt);

         /**
          * @brief Create the correct linear solver
          */
         void createSolver(Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp);

         /**
          * @brief Compute (coupled) matrices
          */
         void createMatrices(Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp);

         /**
          * @brief Initialise the solution
          *
          * @param scalEq Scalar equations
          * @param vectEq Vector equations
          */
         void initSolution(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq);

         /**
          * @brief Update equation input to timestepper
          *
          * @param scalEq Scalar equations
          * @param vectEq Vector equations
          */
         void getInput(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq);

         /**
          * @brief Solve all the linear systems
          */
         void solveSystems();

         /**
          * @brief Update equation unkowns with timestepper output 
          *
          * @param scalEq Shared scalar equations
          * @param vectEq Shared vector equations
          */
         void transferOutput(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq);

         /**
          * @brief Build the solver matrix
          *
          * @param solDIt  Real solver iterator
          * @param spEq    Shared pointer to equation
          * @param comp    Field component
          * @param idx     Matrix index
          * @param isLhs   Flag to update LHS and RHS time dependent matrix
          */
         virtual void buildSolverMatrix(const int matIdx, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx) = 0;

         /**
          * @brief Build the solver matrix
          *
          * @param solZIt  Complex solver iterator
          * @param spEq    Shared pointer to equation
          * @param comp    Field component
          * @param idx     Matrix index
          * @param isLhs   Flag to update LHS and RHS time dependent matrix
          */
         virtual void buildSolverMatrix(SparseMatrixZ& solverMatrix, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx) = 0;

         /**
          * @brief Current solver step
          */
         int   mStep;

         /**
          * @brief Vector of (coupled) real linear solvers
          */
         std::vector<SharedSparseDLinearSolver> mDSolvers;

         /**
          * @brief Vector of (coupled) complex linear solvers
          */
         std::vector<SharedSparseZLinearSolver> mZSolvers;

      private:
   };

   template <typename TSolverIt> void SparseLinearCoordinatorBase::buildSolverMatrices(Equations::SharedIEquation spEq, const SpectralFieldId id, const TSolverIt solveIt)
   {
      // Number of linear systems
      int nSystems = spEq->couplingInfo(id.second).nSystems();

      // start index for matrices
      int start = this->mStep*nSystems;

      // Start row for storage information
      ArrayI startRow(nSystems);

      // Initialise the linear solver
      if(solveIt->nSystem() == 0)
      {
         // Reserve storage for matrice and initialise vectors
         solveIt->initMatrices(ImExRK3::STEPS*nSystems);

         // Initialise field storage and information
         for(int i = 0; i < nSystems; i++)
         {
            // Create data storage
            solveIt->addStorage(spEq->couplingInfo(id.second).systemN(i), spEq->couplingInfo(id.second).rhsCols(i));
         }
      }

      // Build the solver matrices
      for(int i = 0; i < nSystems; i++)
      {
         // Build LHS solver matrix
         this->buildSolverMatrix(solveIt, spEq, id.second, i);

         // Store the start row
         startRow(i) = spEq->couplingInfo(id.second).fieldIndex()*spEq->couplingInfo(id.second).blockN(i);
      }

      // Store storage information
      solveIt->addInformation(id,startRow);
   }

   template <typename TEquationIt, typename TSolverIt> void SparseLinearCoordinatorBase::getSolverInput(const TEquationIt eqIt, const SpectralFieldId id, const TSolverIt solveIt)
   {
      // Get timestep input
      for(int i = 0; i < solveIt->nSystem(); i++)
      {
         // Copy field values into solver input
         Equations::copyUnknown(*(*eqIt), id.second, solveIt->rRHSData(i), i, solveIt->startRow(id,i));

         // Apply quasi-inverse to nonlinear terms
         Equations::applyQuasiInverse(*(*eqIt), id.second, solveIt->rRHSData(i), i, solveIt->startRow(id,i));

         // Loop over all complex solvers
         for(SolverZ_iterator zIt = this->mZSolvers.begin(); zIt != this->mZSolvers.end(); ++zIt)
         {
            // Loop over all fields
            SparseZLinearSolver::FieldId_range   fRange = zIt->fieldRange();
            for(SparseZLinearSolver::FieldId_iterator  fIt = fRange.first; fIt != fRange.second; ++fIt)
            {
               Equations::addExplicitLinear(*(*eqIt), id.second, solveIt->rRHSData(i), solveIt->startRow(id,i), *fIt, zIt->rRHSData(i), zIt->startRow(*fIt,i), i);
            }
         }

         // Loop over all real solvers
         for(SolverD_iterator dIt = this->mDSolvers.begin(); dIt != this->mDSolvers.end(); ++dIt)
         {
            // Loop over all fields
            SparseDLinearSolver::FieldId_range   fRange = dIt->fieldRange();
            for(SparseDLinearSolver::FieldId_iterator  fIt = fRange.first; fIt != fRange.second; ++fIt)
            {
               Equations::addExplicitLinear(*(*eqIt), id.second, solveIt->rRHSData(i), solveIt->startRow(id,i), *fIt, dIt->rRHSData(i), dIt->startRow(*fIt,i), i);
            }
         }
      }
   }
}
}

#endif // SPARSELINEARCOORDINATORBASE_HPP
