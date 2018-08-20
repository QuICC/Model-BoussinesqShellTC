/**
 * @file fjz.hpp
 * @brief Implementation of the vertical current computation for the Boussinesq F-plane 3DQG model 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef QUICC_EQUATIONS_BOUSSINESQ_PLANE_QGMHDBHH_FJZ_HPP
#define QUICC_EQUATIONS_BOUSSINESQ_PLANE_QGMHDBHH_FJZ_HPP

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

namespace QGmhdBhh {

   /**
    * @brief Implementation of the vertical current computation for the Boussinesq F-plane 3DQG model
    */
   class fjz: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         fjz(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~fjz();
         
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

#endif // QUICC_EQUATIONS_BOUSSINESQ_PLANE_QGMHDBHH_FJZ_HPP
