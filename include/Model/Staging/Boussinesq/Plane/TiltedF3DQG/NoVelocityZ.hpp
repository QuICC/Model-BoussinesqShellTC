/**
 * @file BoussinesqTiltedFPlane3DQGNoVelocityZ.hpp
 * @brief Implementation of the non orthogonal vertical velocity computation for the Boussinesq tilted F-plane 3DQG model 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BOUSSINESQTILTEDFPLANE3DQGNOVELOCITYZ_HPP
#define BOUSSINESQTILTEDFPLANE3DQGNOVELOCITYZ_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "Equations/IScalarEquation.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Implementation of the non orthogonal vertical velocity computation for the Boussinesq tilted F-plane 3DQG model
    */
   class BoussinesqTiltedFPlane3DQGNoVelocityZ: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          * @param Solver timing
          */
         BoussinesqTiltedFPlane3DQGNoVelocityZ(SharedEquationParameters spEqParams, const SolveTiming::Id time);

         /**
          * @brief Simple empty destructor
          */
         virtual ~BoussinesqTiltedFPlane3DQGNoVelocityZ();

         /**
          * @brief Compute the nonlinear interaction term
          *
          * @param rNLComp Nonlinear term component
          * @param id      ID of the component (allows for a more general implementation)
          */
         virtual void computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const;
         
      protected:
         /**
          * @brief Set variable requirements
          */
         virtual void setRequirements();

         /**
          * @brief Set the equation coupling information
          */
         virtual void setCoupling();

      private:
   };

}
}

#endif // BOUSSINESQTILTEDFPLANE3DQGNOVELOCITYZ_HPP
