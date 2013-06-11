/** \file BoussinesqBetaSphGTransport.hpp
 *  \brief Implementation of the transport equation for the Boussinesq beta model with spherical gravity
 */

#ifndef BOUSSINESQBETASPHGTRANSPORT_HPP
#define BOUSSINESQBETASPHGTRANSPORT_HPP

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
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/IBoussinesqBetaSphGScalarEquation.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * \brief Implementation of the transport equation for the Boussinesq beta model with spherical gravity
    */
   class BoussinesqBetaSphGTransport: public IBoussinesqBetaSphGScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         BoussinesqBetaSphGTransport(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~BoussinesqBetaSphGTransport();

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

      private:
   };

}
}

#endif // BOUSSINESQBETASPHGTRANSPORT_HPP
