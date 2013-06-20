/** \file PolySetup.cpp
 *  \brief Source of polynomial transform setup class
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "PolynomialTransforms/PolySetup.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   PolySetup::PolySetup(const int size, const int howmany, const int specSize, const std::vector<ArrayI>& fast, const ArrayI& slow)
      : TransformSetup(size, howmany, specSize), mFast(fast), mSlow(slow)
   {
   }

   PolySetup::~PolySetup()
   {
   }

   void PolySetup::init()
   {
   }

}
}
