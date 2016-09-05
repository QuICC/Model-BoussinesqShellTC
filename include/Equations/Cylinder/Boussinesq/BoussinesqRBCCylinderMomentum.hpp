/**
 * @file BoussinesqRBCCylinderMomentum.hpp
 * @brief Implementation of the vector momentum equation for Rayleigh-Benard convection in a cylinder (toroidal-poloidal formulation)
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BOUSSINESQRBCCYLINDERMOMENTUM_HPP
#define BOUSSINESQRBCCYLINDERMOMENTUM_HPP

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
#include "Equations/IVectorEquation.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * @brief Implementation of the vector momentum equation for Rayleigh-Benard convection in a cylinder (toroidal-poloidal formulation)
    */
   class BoussinesqRBCCylinderMomentum: public IVectorEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams    Shared equation parameters
          */
         BoussinesqRBCCylinderMomentum(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~BoussinesqRBCCylinderMomentum();

         /**
          * @brief Compute the nonlinear interaction term
          *
          * @param rNLComp Nonlinear term component
          * @param compId  ID of the component (allows for a more general implementation)
          */
         virtual void computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const;

      protected:
         /**
          * @brief Set variable requirements
          */
         virtual void setRequirements();

         /**
          * @brief Set the equation coupling information
          */
         virtual void setCoupling();

         /**
          * @brief Set the nonlinear integration components
          */
         virtual void setNLComponents();

      private:
   };

}
}

#endif // BOUSSINESQRBCCYLINDERMOMENTUM_HPP
