/** \file RandomScalarState.hpp
 *  \brief Implementation of the diagnostic equation to generate a random scalar state
 */

#ifndef RANDOMSCALARSTATE_HPP
#define RANDOMSCALARSTATE_HPP

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
#include "Equations/IScalarDEquation.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * @brief Implementation of the diagnostic equation to generate a random scalar state
    */
   class RandomScalarState: public IScalarDEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          * @paarm name       Name of the field
          */
         RandomScalarState(SharedEquationParameters spEqParams, const PhysicalNames::Id name);

         /**
          * @brief Simple empty destructor
          */
         virtual ~RandomScalarState();

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
          * @brief Set coupling information
          */
         virtual void setCoupling();

      private:
   };

}
}

#endif // RANDOMSCALARSTATE_HPP
