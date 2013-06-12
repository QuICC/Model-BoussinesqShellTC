/** \file BoussinesqPerBetaCylGStreamfunction.hpp
 *  \brief Implementation of the streamfunction equation for the Boussinesq beta model with cylindrical gravity with periodic radius
 */

#ifndef BOUSSINESQPERBETACYLGSTREAMFUNCTION_HPP
#define BOUSSINESQPERBETACYLGSTREAMFUNCTION_HPP

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
    * \brief Implementation of the streamfunction equation for the Boussinesq beta model with cylindrical gravity with periodic radius
    */
   class BoussinesqPerBetaCylGStreamfunction: public IBoussinesqPerBetaCylGScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         BoussinesqPerBetaCylGStreamfunction(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~BoussinesqPerBetaCylGStreamfunction();

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

#endif // BOUSSINESQPERBETACYLGSTREAMFUNCTION_HPP
