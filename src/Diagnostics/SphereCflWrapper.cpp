/** 
 * @file SphereCflWrapper.cpp
 * @brief Source of the CFL constraint wrapper in a full sphere
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Debug includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Diagnostics/SphereCflWrapper.hpp"

// Project includes
//

namespace QuICC {

namespace Diagnostics {

   SphereCflWrapper::SphereCflWrapper(const SharedIVelocityWrapper spVelocity)
      : ISphericalCflWrapper(spVelocity)
   {
   }

   SphereCflWrapper::~SphereCflWrapper()
   {
   }

   MHDFloat SphereCflWrapper::effectiveMaxL(const MHDFloat r) const
   {
      MHDFloat l = 0;
      for(; l < static_cast<MHDFloat>(this->mspVelocity->spRes()->sim()->dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL)); l++)
      {
         if(r < this->jacobiRoot(l))
         {
            break;
         }
      }

      return l;
   }

   MHDFloat SphereCflWrapper::jacobiRoot(const MHDFloat l) const
   {
      MHDFloat rN = this->mspVelocity->spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL) - 1;
      if(l == 0)
      {
         return std::sin(((-0.5) - 1.8557571*std::pow(0.5,1./3.))/(2.0*rN));
      } else
      {
         return std::sin(((l - 0.5) + 1.8557571*std::pow(l-0.5,1./3.))/(2.0*(rN + l/2.0)));
      }
   }

}
}
