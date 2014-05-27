/**
 * @file BoussinesqShellVelocity.hpp
 * @brief Implementation of the Navier-Stokes equation for the Boussinesq spherical shell model 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BOUSSINESQSHELLVELOCITY_HPP
#define BOUSSINESQSHELLVELOCITY_HPP

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
    * @brief Implementation of the Navier-Stokes equation for the Boussinesq spherical shell model 
    */
   class BoussinesqShellVelocity: public IVectorEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param pyName     Python script name
          * @param spEqParams  Shared equation parameters
          */
         BoussinesqShellVelocity(const std::string& pyName, SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~BoussinesqShellVelocity();

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

#endif // BOUSSINESQSHELLVELOCITY_HPP
