/** \file BoussinesqPerBetaCylGTransport.hpp
 *  \brief Implementation of the transport equation for the Boussinesq beta model with cylindrical gravity with periodic radius
 */

#ifndef BOUSSINESQPERBETACYLGTRANSPORT_HPP
#define BOUSSINESQPERBETACYLGTRANSPORT_HPP

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
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/IBoussinesqPerBetaCylGScalarEquation.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * \brief Implementation of the transport equation for the Boussinesq beta model with cylindrical gravity with periodic radius
    */
   class BoussinesqPerBetaCylGTransport: public IBoussinesqPerBetaCylGScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         BoussinesqPerBetaCylGTransport(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~BoussinesqPerBetaCylGTransport();

         /**
          * @brief Compute the nonlinear interaction term
          *
          * @param rNLComp Nonlinear term component
          */
         virtual void computeNonlinear(Datatypes::PhysicalScalarType& rNLComp) const;

      protected:
         /**
          * @brief Set variable requirements
          */
         virtual void setRequirements();

      private:
   };

}
}

#endif // BOUSSINESQPERBETACYLGTRANSPORT_HPP
