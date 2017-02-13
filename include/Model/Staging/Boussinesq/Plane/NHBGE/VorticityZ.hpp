/**
 * @file BoussinesqFPlaneNHBGEVorticityZ.hpp
 * @brief Implementation of the vertical vorticity computation for the Boussinesq F-plane NHBGE model 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BOUSSINESQFPLANENHBGEVORTICITYZ_HPP
#define BOUSSINESQFPLANENHBGEVORTICITYZ_HPP

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
    * @brief Implementation of the vertical vorticity computation for the Boussinesq F-plane NHBGE model
    */
   class BoussinesqFPlaneNHBGEVorticityZ: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         BoussinesqFPlaneNHBGEVorticityZ(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~BoussinesqFPlaneNHBGEVorticityZ();
         
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

#endif // BOUSSINESQFPLANENHBGEVORTICITYZ_HPP
