/** 
 * @file Momentum.cpp
 * @brief Source of the implementation of the vector Navier-Stokes equation in the Boussinesq spherical Couette in a spherical shell model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

/// Define small macros allowing to convert to string
#define MAKE_STR_X( _P ) # _P
#define MAKE_STR( _P ) MAKE_STR_X( _P )

// Configuration includes
//
#include "TypeSelectors/TransformSelector.hpp"

// System includes
//

// External includes
//

// Class include
//
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Shell/Free/Momentum.hpp )

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "PhysicalOperators/Cross.hpp"
#include "PhysicalOperators/SphericalCoriolis.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Shell {

namespace Free {

   Momentum::Momentum(SharedEquationParameters spEqParams)
      : MomentumBase(spEqParams)
   {
   }

   Momentum::~Momentum()
   {
   }

   Datatypes::SpectralScalarType::PointType Momentum::boundaryValue(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const
   {  
      if(compId == FieldComponents::Spectral::TOR)
      {

         return Datatypes::SpectralScalarType::PointType(0.0);
      } else
      {
         throw std::logic_error("Poloidal component should not set inhomogeneous boundary condition");

         return Datatypes::SpectralScalarType::PointType(0.0);
      }
   }

}
}
}
}
}
