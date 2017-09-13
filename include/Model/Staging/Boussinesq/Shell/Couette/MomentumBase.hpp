/**
 * @file MomentumBase.hpp
 * @brief Base implementation of the vector Navier-Stokes equation for the Boussinesq spherical Couette in a spherical shell
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SHELL_COUETTE_MOMENTUMBASE_HPP
#define QUICC_MODEL_BOUSSINESQ_SHELL_COUETTE_MOMENTUMBASE_HPP

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

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Shell {

namespace Couette {

   /**
    * @brief Implementation of the vector Navier-Stokes equation for the Boussinesq spherical Couette in a spherical shell
    */
   class MomentumBase: public IVectorEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         MomentumBase(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~MomentumBase();

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

         /**
          * @brief Set the nonlinear integration components
          */
         virtual void setNLComponents();

         /**
          * @brief Storage for the cos(theta) grid values (if required)
          */
         Array mCosTheta;

         /**
          * @brief Storage for the sin(theta) grid values (if required)
          */
         Array mSinTheta;

      private:
   };

}
}
}
}
}

#endif // QUICC_MODEL_BOUSSINESQ_SHELL_COUETTE_MOMENTUMBASE_HPP
