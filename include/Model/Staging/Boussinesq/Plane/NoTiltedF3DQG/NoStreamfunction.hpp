/**
 * @file NoStreamfunction.hpp
 * @brief Implementation of the non orthogonal streamfunction computation for the Boussinesq tilted F-plane 3DQG model 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef QUICC_EQUATIONS_BOUSSINESQ_PLANE_NOTILTEDF3DQG_NOSTREAMFUNTION_HPP
#define QUICC_EQUATIONS_BOUSSINESQ_PLANE_NOTILTEDF3DQG_NOSTREAMFUNTION_HPP

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
    * @brief Implementation of the non orthogonal streamfunction computation for the Boussinesq tilted F-plane 3DQG model
    */
   class NoStreamfunction: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         NoStreamfunction(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~NoStreamfunction();

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

#endif // QUICC_EQUATIONS_BOUSSINESQ_PLANE_NOTILTEDF3DQG_NOSTREAMFUNTION_HPP
