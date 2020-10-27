/** 
 * @file ShellCflWrapper.cpp
 * @brief Source of the CFL constraint wrapper in a spherical shell
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
#include "Diagnostics/ShellCflWrapper.hpp"

// Project includes
//

namespace QuICC {

namespace Diagnostics {

   ShellCflWrapper::ShellCflWrapper(const SharedIVectorWrapper spVelocity)
      : ISphericalCflWrapper(spVelocity)
   {
   }

   ShellCflWrapper::~ShellCflWrapper()
   {
   }

   MHDFloat ShellCflWrapper::effectiveMaxL(const MHDFloat r) const
   {
      return static_cast<MHDFloat>(this->mspVelocity->spRes()->sim()->dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL) - 1);
   }

}
}
