/** \file Timestepper.hpp
 *  \brief Implementation of a general timestepper structure
 */

#ifndef TIMESTEPPER_HPP
#define TIMESTEPPER_HPP

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
#include "Timesteppers/EquationDTimestepper.hpp"
#include "Timesteppers/EquationZTimestepper.hpp"
#include "Equations/IScalarEquation.hpp"
#include "Equations/IVectorEquation.hpp"
#include "Timesteppers/ImExRK3.hpp"

namespace GeoMHDiSCC {

namespace Timestep {

   /**
    * @brief Implementation of general timestepper structure
    */
   class Timestepper: public SparseLinearCoordinatorBase
   {
      public:
         /**
          * @brief Constructor
          */
         Timestepper();

         /**
          * @brief Destructor
          */
         ~Timestepper();

         /**
          * @brief Initialise timestepper
          *
          * @param dt      Initial timestep value
          * @param scalEq  Shared scalar equations
          * @param vectEq  Shared vector equations
          */
         void init(const MHDFloat dt, const ScalarEquationRangeType& scalEq, const VectorEquationRangeType& vectEq);

         /**
          * @brief Adapt the timestep used
          *
          * @param cfl     CFL condition
          * @param scalEq  Shared scalar equations
          * @param vectEq  Shared vector equations
          */
         void adaptTimestep(const MHDFloat cfl, const ScalarEquationRangeType& scalEq, const VectorEquationRangeType& vectEq);

         /**
          * @brief Update control status
          */
         void update();

         /**
          * @brief Compute (partial) forward step
          *
          * @param scalEq Shared scalar equations
          * @param vectEq Shared vector equations
          */
         void stepForward(const ScalarEquationRangeType& scalEq, const VectorEquationRangeType& vectEq);

         /**
          * @brief Get current simulation time
          */
         MHDFloat time() const;

         /**
          * @brief Get current simulation timestep
          */
         MHDFloat timestep() const;
         
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
          * @brief Build the solver matrices independenlty of solver type
          */
         template <typename TSolverIt> void buildSolverMatrix(Equations::SharedIEquation spEq, const SpectralFieldId id, const TSolverIt solveIt);

         /**
          * @brief Update time dependence
          */
         void updateMatrices();

         /**
          * @brief Compute the RHS of all linear systems
          */
         void computeRHS();

         /**
          * @brief Build the time matrix
          *
          * @param timeMatrix Storage for time matrix
          * @param spEq       Shared pointer to equation
          * @param comp       Field component
          * @param idx        Matrix index
          */
         void buildTimeMatrix(SparseMatrix& timeMatrix, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx);

         /**
          * @brief Build the time matrix
          *
          * @param timeMatrix Storage for the time matrix
          * @param spEq       Shared pointer to equation
          * @param comp       Field component
          * @param idx        Matrix index
          */
         void buildTimeMatrix(SparseMatrixZ& timeMatrix, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx);

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

         /**
          * @brief Maximum timestep jump per step (See Soederlind)
          */
         const MHDFloat mcMaxJump;

         /**
          * @brief No update window for timestep increase
          */
         const MHDFloat mcUpWindow;

         /**
          * @brief Minimal timestep allowed before simulation abort
          */
         const MHDFloat mcMinDt;

         /**
          * @brief Previous timestep length
          */
         MHDFloat mOldDt;

         /**
          * @brief Timestep length
          */
         MHDFloat mDt;

         /**
          * @brief Current time
          */
         MHDFloat mTime;
   };

   template <typename TSolverIt> void Timestepper::buildSolverMatrix(Equations::SharedIEquation spEq, const SpectralFieldId id, const TSolverIt solveIt)
   {
      // Number of linear systems
      int nSystems = spEq->couplingInfo(id.second).nSystems();

      // start index for matrices
      int start = this->mStep*nSystems;

      // Start row for storage information
      ArrayI startRow(nSystems);

      // Initialise the equation stepper
      if(solveIt->nSystem() == 0)
      {
         // Reserve storage for matrice and initialise vectors
         solveIt->initMatrices(ImExRK3::STEPS*nSystems);

         // Initialise field storage and information
         for(int i = 0; i < nSystems; i++)
         {
            // Create RHS and solution data storage
            solveIt->addStorage(spEq->couplingInfo(id.second).systemN(i), spEq->couplingInfo(id.second).rhsCols(i));

            // Build time matrix
            this->buildTimeMatrix(solveIt->rTMatrix(i), spEq, id.second, i);
         }
      }

      // Build the solver matrices
      for(int i = 0; i < nSystems; i++)
      {
         // Build LHS solver matrix
         this->buildSolverMatrix(solveIt->rLHSMatrix(start+i), spEq, id.second, i, true);

         // Build RHS solver matrix
         this->buildSolverMatrix(solveIt->rRHSMatrix(start+i), spEq, id.second, i, false);

         // Store the start row
         startRow(i) = spEq->couplingInfo(id.second).fieldIndex()*spEq->couplingInfo(id.second).blockN(i);
      }

      // Store storage information
      solveIt->addInformation(id,startRow);
   }
}
}

#endif // TIMESTEPPER_HPP
