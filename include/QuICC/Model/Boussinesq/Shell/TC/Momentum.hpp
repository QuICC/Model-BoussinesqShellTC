/**
 * @file Momentum.hpp
 * @brief Implementation of the vector Navier-Stokes equation for the Boussinesq thermal convection spherical shell
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SHELL_TC_MOMENTUM_HPP
#define QUICC_MODEL_BOUSSINESQ_SHELL_TC_MOMENTUM_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Equations/IVectorEquation.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Shell {

namespace TC {

   /**
    * @brief Implementation of the vector Navier-Stokes equation for the Boussinesq thermal convection in a spherical shell
    */
   class Momentum: public IVectorEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         Momentum(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme);

         /**
          * @brief Simple empty destructor
          */
         virtual ~Momentum();

         /**
          * @brief Initialize nonlinear interaction kernel
          */
         virtual void initNLKernel(const bool force = false);

      protected:
         /**
          * @brief Set variable requirements
          */
         virtual void setRequirements();

         /**
          * @brief Set the equation coupling information
          */
         virtual void setCoupling();

         /**
          * @brief Set the nonlinear integration components
          */
         virtual void setNLComponents();

      private:
   };

}
}
}
}
}

#endif // QUICC_MODEL_BOUSSINESQ_SHELL_TC_MOMENTUM_HPP
