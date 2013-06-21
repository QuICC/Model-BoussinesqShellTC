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

   PolySetup::PolySetup(const int size, const int howmany, const int specSize, const std::vector<ArrayI>& fast, const ArrayI& slow, const ArrayI& mult)
      : TransformSetup(size, howmany, specSize), mFast(fast), mSlow(slow), mMult(mult)
   {
   }

   PolySetup::~PolySetup()
   {
   }

   const std::vector<ArrayI>& PolySetup::fast() const
   {
      return this->mFast;
   }

   const ArrayI& PolySetup::slow() const
   {
      return this->mSlow;
   }

   const ArrayI& PolySetup::mult() const
   {
      return this->mMult;
   }

}
}
