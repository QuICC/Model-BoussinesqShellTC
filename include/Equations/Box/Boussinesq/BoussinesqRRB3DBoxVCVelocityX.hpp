/**
 * @file BoussinesqRRB3DBoxVCVelocityX.hpp
 * @brief Implementation of the momentum equation for the X component for rotating Rayleigh-Benard convection in a 3D box (velocity-continuity formulation)
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BOUSSINESQRRB3DBOXVCVELOCITYX_HPP
#define BOUSSINESQRRB3DBOXVCVELOCITYX_HPP

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

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * @brief Implementation of the momentum equation for the X component for rotating Rayleigh-Benard convection in a 3D box (velocity-continuity formulation)
    */
   class BoussinesqRRB3DBoxVCVelocityX: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams    Shared equation parameters
          */
         BoussinesqRRB3DBoxVCVelocityX(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~BoussinesqRRB3DBoxVCVelocityX();

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

#endif // BOUSSINESQRRB3DBOXVCVELOCITYX_HPP
