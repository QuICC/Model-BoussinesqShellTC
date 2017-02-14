/**
 * @file Induction.hpp
 * @brief Implementation of the vector induction equation for the Boussinesq thermal convection dynamo sphere
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef QUICC_EQUATIONS_BOUSSINESQ_SPHERE_PRECESSIONDYNAMO_INDUCTION_HPP
#define QUICC_EQUATIONS_BOUSSINESQ_SPHERE_PRECESSIONDYNAMO_INDUCTION_HPP

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

namespace Sphere {

namespace PrecessionDynamo {

   /**
    * @brief Implementation of the vector induction equation for the Boussinesq thermal convection dynamo in a sphere
    */
   class Induction: public IVectorEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         Induction(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~Induction();

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

      private:
   };

}
}
}
}
}

#endif // QUICC_EQUATIONS_BOUSSINESQ_SPHERE_PRECESSIONDYNAMO_INDUCTION_HPP
