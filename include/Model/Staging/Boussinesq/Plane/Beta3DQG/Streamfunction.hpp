/**
 * @file BoussinesqBeta3DQGPerStreamfunction.hpp
 * @brief Implementation of the streamfunction equation for the periodic Boussinesq Beta 3DQG model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BOUSSINESQBETA3DQGPERSTREAMFUNCTION_HPP
#define BOUSSINESQBETA3DQGPERSTREAMFUNCTION_HPP

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
    * @brief Implementation of the streamfunction equation for the periodic Boussinesq Beta 3DQG model
    */
   class BoussinesqBeta3DQGPerStreamfunction: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         BoussinesqBeta3DQGPerStreamfunction(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~BoussinesqBeta3DQGPerStreamfunction();

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

#endif // BOUSSINESQBETA3DQGPERSTREAMFUNCTION_HPP
