/**
 * @file BoussinesqRBAnnulusVCMomentum.hpp
 * @brief Implementation of the vector momentum equation for Rayleigh-Benard convection in a cylindrical annulus (velocity-continuity formulation)
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BOUSSINESQRBANNULUSVCMOMENTUM_HPP
#define BOUSSINESQRBANNULUSVCMOMENTUM_HPP

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
    * @brief Implementation of the vector momentum equation for Rayleigh-Benard convection in a cylindrical annulus (velocity-continuity formulation)
    */
   class BoussinesqRBAnnulusVCMomentum: public IVectorEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams    Shared equation parameters
          */
         BoussinesqRBAnnulusVCMomentum(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~BoussinesqRBAnnulusVCMomentum();

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

#endif // BOUSSINESQRBANNULUSVCMOMENTUM_HPP
