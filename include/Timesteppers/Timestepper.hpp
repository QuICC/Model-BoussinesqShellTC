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
#include "Timesteppers/TimestepCoupling.hpp"
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
          * @brief Compute the type of the equation stepper (real or complex)
          */
         void getEqStepperType(Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp);

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
          * @brief Build the LHS matrix triplets
          *
          * @param triplets   Storage for the triplets
          * @param spEq       Shared pointer to equation
          * @param comp       Field component
          * @param idx        Matrix index
          * @param nC         Number of coupled fields
          * @param cRow       Row index of coupled field
          */
         int buildLHSMatrix(std::vector<Triplet>& triplets, Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const int nC = 0, const int cRow = 0);

         /**
          * @brief Build the LHS matrix triplets
          *
          * @param triplets   Storage for the triplets
          * @param spEq       Shared pointer to equation
          * @param comp       Field component
          * @param idx        Matrix index
          * @param nC         Number of coupled fields
          * @param cRow       Row index of coupled field
          */
         int buildLHSMatrix(std::vector<TripletZ>& triplets, Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const int nC = 0, const int cRow = 0);

         /**
          * @brief Build the RHS matrix triplets
          *
          * @param triplets   Storage for the triplets
          * @param spEq       Shared pointer to equation
          * @param comp       Field component
          * @param idx        Matrix index
          * @param nC         Number of coupled fields
          * @param cRow       Row index of coupled field
          */
         void buildRHSMatrix(std::vector<Triplet>& triplets, Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const int nC = 0, const int cRow = 0);

         /**
          * @brief Build the RHS matrix triplets
          *
          * @param triplets   Storage for the triplets
          * @param spEq       Shared pointer to equation
          * @param comp       Field component
          * @param idx        Matrix index
          * @param nC         Number of coupled fields
          * @param cRow       Row index of coupled field
          */
         void buildRHSMatrix(std::vector<TripletZ>& triplets, Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const int nC = 0, const int cRow = 0);

         /**
          * @brief Update the LHS matrix triplets
          *
          * @param triplets   Storage for the triplets
          * @param spEq       Shared pointer to equation
          * @param comp       Field component
          * @param idx        Matrix index
          * @param nC         Number of coupled fields
          * @param cRow       Row index of coupled field
          */
         void updateLHSMatrix(std::vector<Triplet>& triplets, Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const int nC = 0, const int cRow = 0);

         /**
          * @brief Update the LHS matrix triplets
          *
          * @param triplets   Storage for the triplets
          * @param spEq       Shared pointer to equation
          * @param comp       Field component
          * @param idx        Matrix index
          * @param nC         Number of coupled fields
          * @param cRow       Row index of coupled field
          */
         void updateLHSMatrix(std::vector<TripletZ>& triplets, Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const int nC = 0, const int cRow = 0);

         /**
          * @brief Update the RHS matrix triplets
          *
          * @param triplets   Storage for the triplets
          * @param spEq       Shared pointer to equation
          * @param comp       Field component
          * @param idx        Matrix index
          * @param nC         Number of coupled fields
          * @param cRow       Row index of coupled field
          */
         void updateRHSMatrix(std::vector<Triplet>& triplets, Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const int nC = 0, const int cRow = 0);

         /**
          * @brief Update the RHS matrix triplets
          *
          * @param triplets   Storage for the triplets
          * @param spEq       Shared pointer to equation
          * @param comp       Field component
          * @param idx        Matrix index
          * @param nC         Number of coupled fields
          * @param cRow       Row index of coupled field
          */
         void updateRHSMatrix(std::vector<TripletZ>& triplets, Equations::SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const int nC = 0, const int cRow = 0);

         /**
          * @brief Add real triplets corresponding to sparse real matrix
          */
         void addTriplets(std::vector<Triplet>& triplets, const SparseMatrix& mat, const int rowShift, const int colShift, const MHDFloat c);

         /**
          * @brief Add complex triplets corresponding to sparse complex matrix
          */
         void addTriplets(std::vector<TripletZ>& triplets, const SparseMatrixZ& mat, const int rowShift, const int colShift, const MHDComplex c);

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
          * @brief Timestep equation coupling
          */
         TimestepCoupling mTimeCoupling;

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
