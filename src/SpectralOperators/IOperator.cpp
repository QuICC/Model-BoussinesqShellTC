/** 
 * @file IOperator.cpp
 * @brief Source of the interface for spectral operators with quasi inverses
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "SpectralOperators/IOperator.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Spectral {

   IOperator::IOperator(const int basisN)
      : UnitOperator(basisN)
   {
   }

   IOperator::~IOperator()
   {
   }

}
}
