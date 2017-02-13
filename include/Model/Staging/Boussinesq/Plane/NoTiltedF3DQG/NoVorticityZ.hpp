/**
 * @file BoussinesqNoTiltedFPlane3DQGNoVorticityZ.hpp
 * @brief Implementation of the non orthogonal vertical vorticity computation for the Boussinesq tilted F-plane 3DQG model 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BOUSSINESQNOTILTEDFPLANE3DQGNOVORTICITYZ_HPP
#define BOUSSINESQNOTILTEDFPLANE3DQGNOVORTICITYZ_HPP

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
    * @brief Implementation of the non orthogonal vertical vorticity computation for the Boussinesq tilted F-plane 3DQG model
    */
   class BoussinesqNoTiltedFPlane3DQGNoVorticityZ: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         BoussinesqNoTiltedFPlane3DQGNoVorticityZ(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~BoussinesqNoTiltedFPlane3DQGNoVorticityZ();
         
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

#endif // BOUSSINESQNOTILTEDFPLANE3DQGNOVORTICITYZ_HPP
