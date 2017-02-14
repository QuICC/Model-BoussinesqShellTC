/**
 * @file NoVorticityZ.hpp
 * @brief Implementation of the non orthogonal vertical vorticity computation for the Boussinesq tilted F-plane 3DQG model 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef QUICC_EQUATIONS_BOUSSINESQ_PLANE_NOTILTEDF3DQG_NOVORTICITYZ_HPP
#define QUICC_EQUATIONS_BOUSSINESQ_PLANE_NOTILTEDF3DQG_NOVORTICITYZ_HPP

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

namespace NoTiltedF3DQG {

   /**
    * @brief Implementation of the non orthogonal vertical vorticity computation for the Boussinesq tilted F-plane 3DQG model
    */
   class NoVorticityZ: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         NoVorticityZ(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~NoVorticityZ();
         
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

#endif // QUICC_EQUATIONS_BOUSSINESQ_PLANE_NOTILTEDF3DQG_NOVORTICITYZ_HPP
