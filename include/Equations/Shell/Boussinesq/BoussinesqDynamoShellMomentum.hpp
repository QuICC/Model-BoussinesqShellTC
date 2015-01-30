/**
 * @file BoussinesqDynamoShellMomentum.hpp
 * @brief Implementation of the vector Navier-Stokes equation for the Boussinesq thermal convection dynamo spherical shell
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BOUSSINESQDYNAMOSHELLMOMENTUM_HPP
#define BOUSSINESQDYNAMOSHELLMOMENTUM_HPP

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

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * @brief Implementation of the vector Navier-Stokes equation for the Boussinesq thermal convection dynamo in a spherical shell
    */
   class BoussinesqDynamoShellMomentum: public IVectorEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         BoussinesqDynamoShellMomentum(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~BoussinesqDynamoShellMomentum();

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

#endif // BOUSSINESQDYNAMOSHELLMOMENTUM_HPP
