/** 
 * @file PolySetup.cpp
 * @brief Source of polynomial transform setup class
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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

namespace QuICC {

namespace Transform {

   PolySetup::PolySetup(const int size, const int howmany, const int specSize, const std::vector<ArrayI>& fast, const ArrayI& slow, const ArrayI& mult, const int padSize)
      : TransformSetup(size, howmany, specSize), mFast(fast), mSlow(slow), mMult(mult), mPadSize(padSize)
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

   int PolySetup::padSize() const
   {
      return this->mPadSize;
   }

}
}
