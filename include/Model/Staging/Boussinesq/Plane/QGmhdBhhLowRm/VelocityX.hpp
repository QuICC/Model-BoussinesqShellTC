/**
 * @file VelocityX.hpp
 * @brief Implementation of the vertical vorticity computation for the Boussinesq F-plane 3DQG model 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef QUICC_EQUATIONS_BOUSSINESQ_PLANE_QGMHDBHHLOWRM_VELOCITYX_HPP
#define QUICC_EQUATIONS_BOUSSINESQ_PLANE_QGMHDBHHLOWRM_VELOCITYX_HPP

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

namespace QGmhdBhhLowRm {

   /**
    * @brief Implementation of the vertical vorticity computation for the Boussinesq F-plane 3DQG model
    */
   class VelocityX: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         VelocityX(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~VelocityX();
         
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

#endif // QUICC_EQUATIONS_BOUSSINESQ_PLANE_QGMHDBHHLOWRM_VELOCITYX_HPP
