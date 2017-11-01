/**
 * @file Momentum.hpp
 * @brief Implementation of the vector Navier-Stokes equation for the Boussinesq spherical Couette in a spherical shell
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SHELL_COUETTE_MOMENTUM_HPP
#define QUICC_MODEL_BOUSSINESQ_SHELL_COUETTE_MOMENTUM_HPP

/// Define small macros allowing to convert to string
#define MAKE_STR_X( _P ) # _P
#define MAKE_STR( _P ) MAKE_STR_X( _P )

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
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Shell/Couette/MomentumBase.hpp )

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Shell {

namespace Couette {

   /**
    * @brief Implementation of the vector Navier-Stokes equation for the Boussinesq spherical Couette in a spherical shell
    */
   class Momentum: public MomentumBase
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         Momentum(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~Momentum();

         /**
          * @brief Set inhomogeneous boundary condition
          */
         Datatypes::SpectralScalarType::PointType boundaryValue(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const;

      protected:

      private:
   };

}
}
}
}
}

#endif // QUICC_MODEL_BOUSSINESQ_SHELL_COUETTE_MOMENTUM_HPP
