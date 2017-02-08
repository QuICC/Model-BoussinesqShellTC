/**
 * @file BoussinesqPrecessionSphereMomentum.hpp
 * @brief Implementation of the vector Navier-Stokes equation for the Boussinesq precession in a sphere
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BOUSSINESQPRECESSIONSPHEREMOMENTUM_HPP
#define BOUSSINESQPRECESSIONSPHEREMOMENTUM_HPP

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
    * @brief Implementation of the vector Navier-Stokes equation for the Boussinesq precession in a sphere
    */
   class BoussinesqPrecessionSphereMomentum: public IVectorEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         BoussinesqPrecessionSphereMomentum(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~BoussinesqPrecessionSphereMomentum();

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

#endif // BOUSSINESQPRECESSIONSPHEREMOMENTUM_HPP
