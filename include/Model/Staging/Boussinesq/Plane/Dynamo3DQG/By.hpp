/**
 * @file BoussinesqDynamo3DQGBy.hpp
 * @brief Implementation of the vertical vorticity computation for the Boussinesq F-plane 3DQG model 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BOUSSINESQDYNAMO3DQGBY_HPP
#define BOUSSINESQDYNAMO3DQGBY_HPP

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
   class BoussinesqDynamo3DQGBy: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         BoussinesqDynamo3DQGBy(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~BoussinesqDynamo3DQGBy();
         
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

#endif // BOUSSINESQDYNAMO3DQGBY_HPP
