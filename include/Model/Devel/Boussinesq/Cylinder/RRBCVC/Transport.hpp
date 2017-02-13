/**
 * @file BoussinesqRRBCCylinderVCTransport.hpp
 * @brief Implementation of the transport equation for rotating Rayleigh-Benard convection in a cylinder (velocity-continuity formulation)
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BOUSSINESQRRBCCYLINDERVCTRANSPORT_HPP
#define BOUSSINESQRRBCCYLINDERVCTRANSPORT_HPP

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
    * @brief Implementation of the transport equation for rotating Rayleigh-Benard convection in a cylinder (velocity-continuity formulation)
    */
   class BoussinesqRRBCCylinderVCTransport: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams    Shared equation parameters
          */
         BoussinesqRRBCCylinderVCTransport(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~BoussinesqRRBCCylinderVCTransport();

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

#endif // BOUSSINESQRRBCCYLINDERVCTRANSPORT_HPP
