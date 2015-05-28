/** 
 * @file TransformSetup.cpp
 * @brief Source of basic transform setup class
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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
      : mFwdSize(size), mHowmany(1), mSpecSize(specSize)
   {
      this->mHowmany(0) = howmany;
   }

   TransformSetup::TransformSetup(const int size, const ArrayI& howmany, const int specSize)
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
      return this->mHowmany.sum();
   }

   int TransformSetup::howmany(const int parity) const
   {
      return this->mHowmany(parity);
   }

   int TransformSetup::specSize() const
   {
      return this->mSpecSize;
   }

   void TransformSetup::setBoxScale(const MHDFloat boxScale)
   {
   }

}
}
