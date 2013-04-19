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

namespace GeoMHDiSCC {

namespace Timestep {

   /**
    * @brief Implementation of general timestepper structure
    */
   class Timestepper
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
         void init(const MHDFloat dt, const std::vector<Equations::SharedIScalarEquation>& scalEq, const std::vector<Equations::SharedIVectorEquation>& vectEq);

         /**
          * @brief Adapt the timestep used
          *
          * @param cfl     CFL condition
          * @param scalEq  Shared scalar equations
          * @param vectEq  Shared vector equations
          */
         void adaptTimestep(const MHDFloat cfl, const std::vector<Equations::SharedIScalarEquation>& scalEq, const std::vector<Equations::SharedIVectorEquation>& vectEq);

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
         void stepForward(const std::vector<Equations::SharedIScalarEquation>& scalEq, const std::vector<Equations::SharedIVectorEquation>& vectEq);

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
          * @brief Create the correct equation steppers
          */
         void createEqStepper(Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp);

         /**
          * @brief Compute (coupled) matrices
          */
         void createMatrices(Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp);

         /**
          * @brief Update (coupled) matrices
          */
         void updateMatrices(Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp);

         /**
          * @brief Initialise the solution
          *
          * @param scalEq Scalar equations
          * @param vectEq Vector equations
          */
         void initSolution(const std::vector<Equations::SharedIScalarEquation>& scalEq, const std::vector<Equations::SharedIVectorEquation>& vectEq);

         /**
          * @brief Update equation input to timestepper
          *
          * @param scalEq Scalar equations
          * @param vectEq Vector equations
          */
         void getInput(const std::vector<Equations::SharedIScalarEquation>& scalEq, const std::vector<Equations::SharedIVectorEquation>& vectEq);

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
         void transferOutput(const std::vector<Equations::SharedIScalarEquation>& scalEq, const std::vector<Equations::SharedIVectorEquation>& vectEq);

         /**
          * @brief Build the solver matrix
          *
          * @param solverMatrix  Storage for solver matrix
          * @param spEq       Shared pointer to equation
          * @param comp       Field component
          * @param idx        Matrix index
          * @param isLhs      Flag to update LHS and RHS time dependent matrix
          */
         void buildSolverMatrix(SparseMatrix& solverMatrix, Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const bool isLhs);

         /**
          * @brief Build the solver matrix
          *
          * @param solverMatrix   Storage for the solver matrix
          * @param spEq       Shared pointer to equation
          * @param comp       Field component
          * @param idx        Matrix index
          * @param isLhs      Flag to update LHS and RHS time dependent matrix
          */
         void buildSolverMatrix(SparseMatrixZ& solverMatrix, Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const bool isLhs);

         /**
          * @brief Update the LHS matrix triplets
          *
          * @param oldTime    Storage for the triplets
          * @param spEq       Shared pointer to equation
          * @param comp       Field component
          * @param idx        Matrix index
          * @param isLhs      Flag to update LHS and RHS time dependent matrix
          */
         void updateTimeMatrix(SparseMatrix& oldTime, Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const bool isLhs);

         /**
          * @brief Update the LHS matrix triplets
          *
          * @param oldTime    Storage for the triplets
          * @param spEq       Shared pointer to equation
          * @param comp       Field component
          * @param idx        Matrix index
          * @param isLhs      Flag to update LHS and RHS time dependent matrix
          */
         void updateTimeMatrix(SparseMatrixZ& oldTime, Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const bool isLhs);

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
}
}

#endif // TIMESTEPPER_HPP
