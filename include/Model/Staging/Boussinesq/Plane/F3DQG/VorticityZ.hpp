/**
 * @file VorticityZ.hpp
 * @brief Implementation of the vertical vorticity computation for the Boussinesq F-plane 3DQG model 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef QUICC_EQUATIONS_BOUSSINESQ_PLANE_F3DQG_VORTICITYZ_HPP
#define QUICC_EQUATIONS_BOUSSINESQ_PLANE_F3DQG_VORTICITYZ_HPP

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
    * @brief Implementation of the vertical vorticity computation for the Boussinesq F-plane 3DQG model
    */
   class VorticityZ: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         VorticityZ(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~VorticityZ();
         
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

#endif // QUICC_EQUATIONS_BOUSSINESQ_PLANE_F3DQG_VORTICITYZ_HPP
