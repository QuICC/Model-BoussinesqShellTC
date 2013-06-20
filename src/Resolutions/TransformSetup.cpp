/** \file TransformSetup.cpp
 *  \brief Source of basic transform setup class
 */

// System includes
//

// External includes
//

// Class include
//
#include "Resolutions/TransformSetup.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

   TransformSetup::TransformSetup(const int size, const int howmany, const int specSize)
      : mFwdSize(size), mHowmany(howmany), mSpecSize(specSize)
   {
   }

   TransformSetup::~TransformSetup()
   {
   }

   int TransformSetup::fwdSize() const
   {
      return this->mFwdSize;
   }

   int TransformSetup::howmany() const
   {
      return this->mHowmany;
   }

   int TransformSetup::specSize() const
   {
      return this->mSpecSize;
   }

}
}
