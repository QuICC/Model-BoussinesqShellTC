/**
 * @file BoussinesqBeta3DQGVorticity.hpp
 * @brief Implementation of the voriticity computation for the Boussinesq Beta 3DQG model 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BOUSSINESQBETA3DQGVORTICITY_HPP
#define BOUSSINESQBETA3DQGVORTICITY_HPP

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
    * @brief Implementation of the voriticity computation for the Boussinesq Beta 3DQG model
    */
   class BoussinesqBeta3DQGVorticity: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param pyName     Python script name
          * @param spEqParams Shared equation parameters
          */
         BoussinesqBeta3DQGVorticity(const std::string& pyName, SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~BoussinesqBeta3DQGVorticity();
         
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

#endif // BOUSSINESQBETA3DQGVORTICITY_HPP
