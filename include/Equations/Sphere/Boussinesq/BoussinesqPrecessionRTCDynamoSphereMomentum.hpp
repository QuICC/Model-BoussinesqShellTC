/**
 * @file BoussinesqPrecessionRTCDynamoSphereMomentum.hpp
 * @brief Implementation of the vector Navier-Stokes equation for the Boussinesq combined precession and thermal convection dynamo sphere
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BOUSSINESQPRECESSIONRTCDYNAMOSPHEREMOMENTUM_HPP
#define BOUSSINESQPRECESSIONRTCDYNAMOSPHEREMOMENTUM_HPP

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
    * @brief Implementation of the vector Navier-Stokes equation for the Boussinesq combined precession and thermal convection dynamo in a sphere
    */
   class BoussinesqPrecessionRTCDynamoSphereMomentum: public IVectorEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         BoussinesqPrecessionRTCDynamoSphereMomentum(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~BoussinesqPrecessionRTCDynamoSphereMomentum();

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
         /**
          * @brief Storage for the r grid values
          */
         Array mR;

         /**
          * @brief Storage for the theta grid values
          */
         Array mTheta;

         /**
          * @brief Storage for the phi grid values
          */
         Array mPhi;
   };

}
}

#endif // BOUSSINESQPRECESSIONRTCDYNAMOSPHEREMOMENTUM_HPP
