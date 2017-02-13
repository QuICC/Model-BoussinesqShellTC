/**
 * @file BoussinesqFPlane3DQGVelocityZ.hpp
 * @brief Implementation of the upright vertical velocity equation for the Boussinesq F-plane 3DQG model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BOUSSINESQFPLANE3DQGVELOCITYZ_HPP
#define BOUSSINESQFPLANE3DQGVELOCITYZ_HPP

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
    * @brief Implementation of the upright vertical velocity equation for the Boussinesq F-plane 3DQG model
    */
   class BoussinesqFPlane3DQGVelocityZ: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         BoussinesqFPlane3DQGVelocityZ(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~BoussinesqFPlane3DQGVelocityZ();
         
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

#endif // BOUSSINESQFPLANE3DQGVELOCITYZ_HPP
