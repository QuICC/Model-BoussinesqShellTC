/** \file Timestepper.hpp
 *  \brief Implementation of a general timestepper structure
 *
 *  \mhdBug Needs test
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

   /**
    * @brief Implementation of general timestepper structure
    */
   class Timestepper
   {
      public:
         /// Typedef for a field ID
         typedef std::pair<PhysicalNames::Id, FieldComponents::Spectral::Id>  FieldIdType;

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
          * @param scalEq Scalar equations
          * @param vectEq Vector equations
          */
         void init(const std::vector<SharedIScalarEquation>& scalEq, const std::vector<SharedIVectorEquation>& vectEq);

         /**
          * @brief Update control status
          */
         void update();

         /**
          * @brief Compute (partial) forward step
          *
          * @param scalEq Scalar equations
          * @param vectEq Vector equations
          */
         void stepForward(const std::vector<SharedIScalarEquation>& scalEq, const std::vector<SharedIVectorEquation>& vectEq);
         
      protected:

      private:
         /**
          * @brief Compute the type of the equation stepper (real or complex)
          */
         void getEqStepperType(SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp, std::map<FieldIdType, int>& coupled, std::vector<bool>& typeInfo);

         /**
          * @brief Create the correct equation steppers
          */
         void createEqStepper(SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp, std::map<FieldIdType, std::pair<bool,int> >& coupled, const std::vector<bool>& typeInfo);

         /**
          * @brief Compute (coupled) matrices
          */
         void createMatrices(SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp);

         /**
          * @brief Initialise the solution
          *
          * @param scalEq Scalar equations
          * @param vectEq Vector equations
          */
         void initSolution(const std::vector<SharedIScalarEquation>& scalEq, const std::vector<SharedIVectorEquation>& vectEq);

         /**
          * @brief Update equation input to timestepper
          *
          * @param scalEq Scalar equations
          * @param vectEq Vector equations
          */
         void getInput(const std::vector<SharedIScalarEquation>& scalEq, const std::vector<SharedIVectorEquation>& vectEq);

         /**
          * @brief Compute the RHS of all linear systems
          */
         void computeRHS();

         /**
          * @brief Solve the all the linear systems
          */
         void solve();

         /**
          * @brief Update equation unkowns with timestepper output 
          *
          * @param scalEq Scalar equations
          * @param vectEq Vector equations
          */
         void transferOutput(const std::vector<SharedIScalarEquation>& scalEq, const std::vector<SharedIVectorEquation>& vectEq);

         /**
          * @brief Build the LHS matrix
          *
          * @param spEq Shared pointer to equation
          * @param comp Field component
          * @param idx  Matrix index
          * @param nC   Number of coupled fields
          * @param cRow Row index of coupled field
          */
         DecoupledZSparse buildLHSMatrix(SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const int nC = 0, const int cRow = 0);

         /**
          * @brief Build the RHS matrix
          *
          * @param spEq Shared pointer to equation
          * @param comp Field component
          * @param idx  Matrix index
          * @param nC   Number of coupled fields
          * @param cRow Row index of coupled field
          */
         DecoupledZSparse buildRHSMatrix(SharedIEvolutionEquation spEq, FieldComponents::Spectral::Id comp, const int idx, const int nC = 0, const int cRow = 0);

         /**
          * @brief Current timestepper step (local)
          */
         int   mStep;

         /**
          * @brief Timestep length
          */
         MHDFloat mDt;

         /**
          * @brief Current time
          */
         MHDFloat mTime;

         /**
          * @brief Map between field ID and corresponding storage
          */
         std::map<FieldIdType,std::pair<bool,int> > mEqInformation;

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

#endif // TIMESTEPPER_HPP
