/** \file SpatialSimulationBase.hpp
 *  \brief Building block for the implementation of a spatial simulation
 */

#ifndef SPATIALSIMULATIONBASE_HPP
#define SPATIALSIMULATIONBASE_HPP

// Configuration includes
//
#include "Simulation/PrepMacros/VariableTypedefsMacro.h"
#include "Simulation/PrepMacros/TransformTypedefsMacro.h"
#include "Simulation/PrepMacros/EquationParametersMacro.h"
#include "Simulation/PrepMacros/TransformGrouperMacro.h"

// System includes
//
#include <map>

// External includes
//

// Project includes
//
#include "Simulation/Enums/PhysicalNames.hpp"
#include "Simulation/System/SystemBase.hpp"
#include "Simulation/Equations/ScalarEquation.hpp"
#include "Simulation/Equations/VectorEquation.hpp"

namespace EPMPhoenix {

   /**
    * \brief Building block for the implementation of a spatial simulation
    */
   class SpatialSimulationBase: public SystemBase
   {
      public:
         /**
          * @brief Simple empty destructor
          */
         virtual ~SpatialSimulationBase() {};

      protected:
         /**
          * @brief Constructor
          */
         SpatialSimulationBase();

         /**
          * @brief Initialise the spatial simulation base
          */
         void initSpatial(const int nCpu, const ArrayI& dim);

         /**
          * @brief Initialise the equations
          */
         virtual void initEquations() = 0;

         /**
          * @brief Setup the equations
          */
         void setupEquations();

         /**
          * @brief Initialise the variables 
          *
          * The required information is gathered from the set of equations
          */
         void initVariables();

         /**
          * @brief Compute the nonlinear terms
          */
         void computeNonlinear();

         /**
          * @brief Get/Set the transform coordinator
          */
         Code::TransformCoordinator&  transformCoordinator();

         /**
          * @brief Storage for scalar equations
          */
         std::vector<SharedScalarEquation> mScalarEquations;

         /**
          * @brief Storage for vector equations
          */
         std::vector<SharedVectorEquation> mVectorEquations;

         /**
          * @brief Map between name and pointer for the scalar variables
          */
         std::map<PhysicalNames::Id, Code::SharedScalarVariable>  mScalarVariables;

         /**
          * @brief Map between name and pointer for the vector variables
          */
         std::map<PhysicalNames::Id, Code::SharedVectorVariable>  mVectorVariables;

         /**
          * @brief Storage for a shared forward transform grouper
          */
         SharedForwardGrouperBase   mspFwdGrouper;

         /**
          * @brief Storage for a shared backward transform grouper
          */
         SharedBackwardGrouperBase   mspBwdGrouper;

      private:
         /**
          * @brief Transform coordinator
          */
         Code::TransformCoordinator mTransformCoordinator;
   };

   inline Code::TransformCoordinator& SpatialSimulationBase::transformCoordinator()
   {
      return this->mTransformCoordinator;
   }
}

#endif // SPATIALSIMULATIONBASE_HPP
