/**
 * @file BoussinesqPrecessionRTCDynamoSphereTransport.hpp
 * @brief Implementation of the transport equation for the Boussinesq combined precession and thermal convection dynamo in a sphere
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BOUSSINESQPRECESSIONRTCDYNAMOSPHERETRANSPORT_HPP
#define BOUSSINESQPRECESSIONRTCDYNAMOSPHERETRANSPORT_HPP

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
    * @brief Implementation of the transport equation for the Boussinesq combined precession and thermal convection dynamo in a sphere 
    */
   class BoussinesqPrecessionRTCDynamoSphereTransport: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         BoussinesqPrecessionRTCDynamoSphereTransport(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~BoussinesqPrecessionRTCDynamoSphereTransport();

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

#endif // BOUSSINESQPRECESSIONRTCDYNAMOSPHERETRANSPORT_HPP
