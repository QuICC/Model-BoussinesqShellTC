/**
 * @file BoussinesqTiltedFPlane3DQGStreamfunction.hpp
 * @brief Implementation of the streamfunction equation for the Boussinesq tilted F-plane 3DQG model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BOUSSINESQTILTEDFPLANE3DQGSTREAMFUNCTION_HPP
#define BOUSSINESQTILTEDFPLANE3DQGSTREAMFUNCTION_HPP

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
    * @brief Implementation of the streamfunction equation for the Boussinesq tilted F-plane 3DQG model
    */
   class BoussinesqTiltedFPlane3DQGStreamfunction: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         BoussinesqTiltedFPlane3DQGStreamfunction(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~BoussinesqTiltedFPlane3DQGStreamfunction();
         
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

#endif // BOUSSINESQTILTEDFPLANE3DQGSTREAMFUNCTION_HPP
