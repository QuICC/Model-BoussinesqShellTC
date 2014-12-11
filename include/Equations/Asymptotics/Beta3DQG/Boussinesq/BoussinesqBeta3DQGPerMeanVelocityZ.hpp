/**
 * @file BoussinesqBeta3DQGPerMeanVelocityZ.hpp
 * @brief Implementation of the mean velocity component for the periodic Boussinesq Beta 3DQG model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BOUSSINESQBETA3DQGPERMEANVELOCITYZ_HPP
#define BOUSSINESQBETA3DQGPERMEANVELOCITYZ_HPP

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
    * @brief Implementation of the kinetic energy for the periodic Boussinesq Beta 3DQG model
    */
   class BoussinesqBeta3DQGPerMeanVelocityZ: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams    Shared equation parameters
          */
         BoussinesqBeta3DQGPerMeanVelocityZ(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~BoussinesqBeta3DQGPerMeanVelocityZ();

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

#endif // BOUSSINESQBETA3DQGPERMEANVELOCITYZ_HPP
