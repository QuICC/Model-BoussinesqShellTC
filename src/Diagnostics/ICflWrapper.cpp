/** 
 * @file ICflWrapper.cpp
 * @brief Source of the interface for the CFL constraint
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
#include "Diagnostics/ICflWrapper.hpp"

// Project includes
//

namespace QuICC {

namespace Diagnostics {

   ICflWrapper::ICflWrapper(const SharedIVectorWrapper spVelocity)
      : mspVelocity(spVelocity)
   {
   }

   ICflWrapper::ICflWrapper(const SharedIVectorWrapper spVelocity, const SharedIVectorWrapper spMagnetic)
      : mspVelocity(spVelocity), mspMagnetic(spMagnetic)
   {
   }

   ICflWrapper::~ICflWrapper()
   {
   }

}
}
