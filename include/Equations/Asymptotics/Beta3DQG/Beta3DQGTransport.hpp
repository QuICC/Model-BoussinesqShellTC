/** \file Beta3DQGTransport.hpp
 *  \brief Implementation of the transport equation for the 3DQG beta model
 */

#ifndef BETA3DQGTRANSPORT_HPP
#define BETA3DQGTRANSPORT_HPP

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
#include "Equations/Asymptotics/Beta3DQG/IBeta3DQGScalarEquation.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * \brief Implementation of the transport equation for the 3DQG beta model
    */
   class Beta3DQGTransport: public IBeta3DQGScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         Beta3DQGTransport(SharedIEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~Beta3DQGTransport();

         /**
          * @brief Compute the nonlinear interaction term
          *
          * @param rNLComp Nonlinear term component
          */
         virtual void computeNonlinear(Datatypes::PhysicalScalarType& rNLComp) const;

         /**
          * @brief Compute the linear term
          *
          * @param rRHS    RHS of timestepping equation
          */
         virtual void computeLinear(Datatypes::SpectralScalarType& rRHS) const;

      protected:
         /**
          * @brief Set variable requirements
          */
         virtual void setRequirements();

      private:
   };

}
}

#endif // BETA3DQGTRANSPORT_HPP
