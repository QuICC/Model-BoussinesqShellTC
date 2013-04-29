/** \file BoussinesqBetaSphGVertical.hpp
 *  \brief Implementation of the vertical velocity equation for the Boussinesq beta model with spherical gravity
 */

#ifndef BOUSSINESQBETASPHGVERTICAL_HPP
#define BOUSSINESQBETASPHGVERTICAL_HPP

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
#include "Equations/Asymptotics/Beta3DQG/IBoussinesqBetaSphGScalarEquation.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * \brief Implementation of the vertical velocity equation for the Boussinesq beta model with spherical gravity
    */
   class BoussinesqBetaSphGVertical: public IBoussinesqBetaSphGScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         BoussinesqBetaSphGVertical(SharedIEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~BoussinesqBetaSphGVertical();

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

#endif // BOUSSINESQBETASPHGVERTICAL_HPP
