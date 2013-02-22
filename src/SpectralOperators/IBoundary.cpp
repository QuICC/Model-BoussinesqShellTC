/** \file IBoundary.cpp
 *  \brief Source of the interface for spectral boundary operators
 */

// System includes
//
#include <assert.h>

// External includes
//

// Class include
//
#include "SpectralOperators/IBoundary.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Spectral {

   IBoundary::IBoundary(const int basisN)
      : mBasisN(basisN)
   {
   }

   IBoundary::~IBoundary()
   {
   }

   void IBoundary::reset(const int basisN)
   {
      // Set dimension
      this->mBasisN = basisN;
   }

   int IBoundary::basisN() const
   {
      return this->mBasisN;
   }

}
}
