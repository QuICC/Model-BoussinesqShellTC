/** 
 * @file CartesianExactStateIds.cpp
 * @brief Source of the implementation of the equation to generate an exact scalar solution in cartesian geometries
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Generator/States/CartesianExactStateIds.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   const MHDFloat CartesianExactStateIds::PCOS = 99999;

   const MHDFloat CartesianExactStateIds::PSIN = -99999;

   MHDFloat CartesianExactStateIds::cos(const MHDFloat amplitude, const MHDFloat mode, const MHDFloat theta)
   {
      return amplitude*std::cos(mode*theta);
   }

   MHDFloat CartesianExactStateIds::sin(const MHDFloat amplitude, const MHDFloat mode, const MHDFloat theta)
   {
      return amplitude*std::sin(mode*theta);
   }

   MHDFloat CartesianExactStateIds::poly(const MHDFloat amplitude, const MHDFloat mode, const MHDFloat x)
   {
      MHDFloat val;

      if(mode == CartesianExactStateIds::PCOS)
      {
         val = CartesianExactStateIds::cos(amplitude,mode,Math::PI*(x-1)/2.0);
      } else if(mode == CartesianExactStateIds::PSIN)
      {
         val = CartesianExactStateIds::sin(amplitude,mode,Math::PI*(x-1)/2.0);
      } else
      {
         val = amplitude*std::pow(x,mode);
      }

      return val;
   }

}
}
