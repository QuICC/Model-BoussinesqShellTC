/**
 * @file BoussinesqDynamoShellInduction.hpp
 * @brief Implementation of the vector induction equation for the Boussinesq thermal convection dynamo spherical shell
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BOUSSINESQDYNAMOSHELLINDUCTION_HPP
#define BOUSSINESQDYNAMOSHELLINDUCTION_HPP

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

   /**
    * @brief Implementation of the vector induction equation for the Boussinesq thermal convection dynamo in a spherical shell
    */
   class BoussinesqDynamoShellInduction: public IVectorEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         BoussinesqDynamoShellInduction(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~BoussinesqDynamoShellInduction();

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

#endif // BOUSSINESQDYNAMOSHELLINDUCTION_HPP
