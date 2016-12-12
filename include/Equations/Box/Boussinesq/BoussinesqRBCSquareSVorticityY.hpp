/**
 * @file BoussinesqRBCSquareSVorticityY.hpp
 * @brief Implementation of the vorticity equation for Rayleigh-Benard convection in a square cavity (2D) (streamfunction formulation)
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BOUSSINESQRBCSQUARESVORTICITYY_HPP
#define BOUSSINESQRBCSQUARESVORTICITYY_HPP

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

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * @brief Implementation of the vorticity equation for Rayleigh-Benard convection in a square cavity (2D) (streamfunction formulation)
    */
   class BoussinesqRBCSquareSVorticityY: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams    Shared equation parameters
          */
         BoussinesqRBCSquareSVorticityY(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~BoussinesqRBCSquareSVorticityY();

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

#endif // BOUSSINESQRBCSQUARESVORTICITYY_HPP
