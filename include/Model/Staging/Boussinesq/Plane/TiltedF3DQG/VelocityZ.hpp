/**
 * @file VelocityZ.hpp
 * @brief Implementation of the upright vertical velocity equation for the Boussinesq tilted F-plane 3DQG model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef QUICC_EQUATIONS_BOUSSINESQ_PLANE_TILTEDF3DQG_VELOCITYZ_HPP
#define QUICC_EQUATIONS_BOUSSINESQ_PLANE_TILTEDF3DQG_VELOCITYZ_HPP

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

namespace Boussinesq {

namespace Plane {

namespace TiltedF3DQG {

   /**
    * @brief Implementation of the upright vertical velocity equation for the Boussinesq tilted F-plane 3DQG model
    */
   class VelocityZ: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         VelocityZ(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~VelocityZ();
         
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
}
}
}

#endif // QUICC_EQUATIONS_BOUSSINESQ_PLANE_TILTEDF3DQG_VELOCITYZ_HPP
