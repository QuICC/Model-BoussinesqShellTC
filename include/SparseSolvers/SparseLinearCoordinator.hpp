/** \file SparseLinearCoordinator.hpp
 *  \brief Implementation of a general sparse linear solver coordinator
 */

#ifndef SPARSELINEARCOORDINATOR_HPP
#define SPARSELINEARCOORDINATOR_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "SparseSolvers/SparseLinearCoordinatorBase.hpp"
#include "SparseLinearSolver/SparseDLinearSolver.hpp"
#include "SparseLinearSolver/SparseZLinearSolver.hpp"
#include "Equations/IScalarEquation.hpp"
#include "Equations/IVectorEquation.hpp"

namespace GeoMHDiSCC {

namespace Solver {

   /**
    * @brief Implementation of general sparse linear solver coordinatore
    */
   class SparseLinearCoordinator: public SparseLinearCoordinatorBase
   {
      public:
         /**
          * @brief Constructor
          */
         SparseLinearCoordinator();

         /**
          * @brief Destructor
          */
         ~SparseLinearCoordinator();

         /**
          * @brief Solve the equations
          *
          * @param scalEq Shared scalar equations
          * @param vectEq Shared vector equations
          */
         void solve(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq);
      protected:
         /**
          * @brief Create a real linear solver
          */
         virtual void addSolverD(const int start);

         /**
          * @brief Create a complex linear solver
          */
         virtual void addSolverZ(const int start);

      private:
         /**
          * @brief Build the solver matrix
          *
          * @param solverMatrix  Storage for solver matrix
          * @param spEq       Shared pointer to equation
          * @param comp       Field component
          * @param idx        Matrix index
          * @param isLhs      Flag to update LHS and RHS time dependent matrix
          */
         void buildSolverMatrix(SparseMatrix& solverMatrix, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const bool isLhs);

         /**
          * @brief Build the solver matrix
          *
          * @param solverMatrix   Storage for the solver matrix
          * @param spEq       Shared pointer to equation
          * @param comp       Field component
          * @param idx        Matrix index
          * @param isLhs      Flag to update LHS and RHS time dependent matrix
          */
         void buildSolverMatrix(SparseMatrixZ& solverMatrix, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const bool isLhs);
   };

   template <typename TSolverIt> void SparseLinearCoordinator::buildSolverMatrix(Equations::SharedIEquation spEq, const SpectralFieldId id, const TSolverIt solveIt)
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
         this->buildSolverMatrix(solveIt->rLHSMatrix(start+i), spEq, id.second, i);

         // Store the start row
         startRow(i) = spEq->couplingInfo(id.second).fieldIndex()*spEq->couplingInfo(id.second).blockN(i);
      }

      // Store storage information
      solveIt->addInformation(id,startRow);
   }
}
}

#endif // SPARSELINEARCOORDINATOR_HPP
