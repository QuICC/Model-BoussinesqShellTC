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
   class Timestepper
   {
      public:
         /// Typedef for an iterator to a complex equation timstepper
         typedef std::vector<EquationZTimestepper>::iterator   eqz_iterator;

         /// Typedef for an iterator to a real equation timstepper
         typedef std::vector<EquationDTimestepper>::iterator   eqd_iterator;

         /// Typedef for a shared scalar equation iterator
         typedef std::vector<Equations::SharedIScalarEquation>::iterator   ScalarEquationIteratorType;

         /// Typedef for a shared vector equation iterator
         typedef std::vector<Equations::SharedIVectorEquation>::iterator   VectorEquationIteratorType;

         /// Typedef for a shared scalar equation range
         typedef std::pair<ScalarEquationIteratorType, ScalarEquationIteratorType>  ScalarEquationRangeType;

         /// Typedef for a shared vector equation range
         typedef std::pair<VectorEquationIteratorType, VectorEquationIteratorType>  VectorEquationRangeType;

         /**
          * @brief Constructor
          */
         Timestepper();

         /**
          * @brief Destructor
          */
         ~Timestepper();

         /**
          * @brief Finished computation of full timestep?
          */
         bool finishedStep() const;

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

      private:
         /**
          * @brief Get the solver input independenlty of solver type
          */
         template <typename TEquationIt, typename TSolverIt> void getSolverInput(const TEquationIt eqIt, const SpectralFieldId id, const TSolverIt solveIt);

         /**
          * @brief Build the solver matrices independenlty of solver type
          */
         template <typename TSolverIt> void buildSolverMatrix(Equations::SharedIEquation spEq, const SpectralFieldId id, const TSolverIt solveIt);

         /**
          * @brief Create the correct equation steppers
          */
         void createEqStepper(Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp);

         /**
          * @brief Compute (coupled) matrices
          */
         void createMatrices(Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp);

         /**
          * @brief Update time dependence
          */
         void updateMatrices();

         /**
          * @brief Initialise the solution
          *
          * @param scalEq Scalar equations
          * @param vectEq Vector equations
          */
         void initSolution(const ScalarEquationRangeType& scalEq, const VectorEquationRangeType& vectEq);

         /**
          * @brief Update equation input to timestepper
          *
          * @param scalEq Scalar equations
          * @param vectEq Vector equations
          */
         void getInput(const ScalarEquationRangeType& scalEq, const VectorEquationRangeType& vectEq);

         /**
          * @brief Compute the RHS of all linear systems
          */
         void computeRHS();

         /**
          * @brief Solve all the linear systems
          */
         void solve();

         /**
          * @brief Update equation unkowns with timestepper output 
          *
          * @param scalEq Shared scalar equations
          * @param vectEq Shared vector equations
          */
         void transferOutput(const ScalarEquationRangeType& scalEq, const VectorEquationRangeType& vectEq);

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
          * @brief Current timestepper step (local)
          */
         int   mStep;

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

         /**
          * @brief Vector of (coupled) real equation timesteppers
          */
         std::vector<EquationDTimestepper> mEqDStepper;

         /**
          * @brief Vector of (coupled) complex equation timesteppers
          */
         std::vector<EquationZTimestepper> mEqZStepper;
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

   template <typename TEquationIt, typename TSolverIt> void Timestepper::getSolverInput(const TEquationIt eqIt, const SpectralFieldId id, const TSolverIt solveIt)
   {
      // Get timestep input
      for(int i = 0; i < solveIt->nSystem(); i++)
      {
         // Copy field values into timestep input
         Equations::copyUnknown(*(*eqIt), id.second, solveIt->rRHSData(i), i, solveIt->startRow(id,i));

         // Apply quasi-inverse to nonlinear terms
         Equations::applyQuasiInverse(*(*eqIt), id.second, solveIt->rRHSData(i), i, solveIt->startRow(id,i));

         // Loop over all complex solvers
         for(std::vector<EquationZTimestepper>::iterator zIt = this->mEqZStepper.begin(); zIt != this->mEqZStepper.end(); ++zIt)
         {
            // Loop over all fields
            EquationZTimestepper::field_iterator_range   fRange = zIt->fieldRange();
            for(EquationZTimestepper::field_iterator  fIt = fRange.first; fIt != fRange.second; ++fIt)
            {
               Equations::addExplicitLinear(*(*eqIt), id.second, solveIt->rRHSData(i), solveIt->startRow(id,i), *fIt, zIt->rRHSData(i), zIt->startRow(*fIt,i), i);
            }
         }

         // Loop over all real solvers
         for(std::vector<EquationDTimestepper>::iterator dIt = this->mEqDStepper.begin(); dIt != this->mEqDStepper.end(); ++dIt)
         {
            // Loop over all fields
            EquationDTimestepper::field_iterator_range   fRange = dIt->fieldRange();
            for(EquationDTimestepper::field_iterator  fIt = fRange.first; fIt != fRange.second; ++fIt)
            {
               Equations::addExplicitLinear(*(*eqIt), id.second, solveIt->rRHSData(i), solveIt->startRow(id,i), *fIt, dIt->rRHSData(i), dIt->startRow(*fIt,i), i);
            }
         }
      }
   }
}
}

#endif // TIMESTEPPER_HPP
