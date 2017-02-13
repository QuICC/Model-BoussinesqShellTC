/**
 * @file BoussinesqTiltedFPlane3DQGVelocityZ.hpp
 * @brief Implementation of the upright vertical velocity equation for the Boussinesq tilted F-plane 3DQG model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BOUSSINESQTILTEDFPLANE3DQGVELOCITYZ_HPP
#define BOUSSINESQTILTEDFPLANE3DQGVELOCITYZ_HPP

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
    * @brief Implementation of the upright vertical velocity equation for the Boussinesq tilted F-plane 3DQG model
    */
   class BoussinesqTiltedFPlane3DQGVelocityZ: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         BoussinesqTiltedFPlane3DQGVelocityZ(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~BoussinesqTiltedFPlane3DQGVelocityZ();
         
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

#endif // BOUSSINESQTILTEDFPLANE3DQGVELOCITYZ_HPP
