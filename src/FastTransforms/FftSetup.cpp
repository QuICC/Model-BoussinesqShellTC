/** \file FftSetup.cpp
 *  \brief Source of FFT setup class
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "FastTransforms/FftSetup.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

   FftSetup::FftSetup(const int size, const int howmany, const int specSize, const bool isMixed)
      : mFwdSize(size), mBwdSize(0), mHowmany(howmany), mSpecSize(specSize), mIsMixed(isMixed), mScale(-1)
   {
      // Set the backward size
      if(this->mIsMixed)
      {
         this->mBwdSize = this->mFwdSize/2 + 1;
      } else
      {
         this->mBwdSize = this->mFwdSize;
      }
   }

   FftSetup::~FftSetup()
   {
   }

   void FftSetup::setScale(const MHDFloat scale)
   {
      this->mScale = scale;
   }

   bool FftSetup::isMixed() const
   {
      return this->mIsMixed;
   }

   int FftSetup::fwdSize() const
   {
      return this->mFwdSize;
   }

   int FftSetup::bwdSize() const
   {
      return this->mBwdSize;
   }

   int FftSetup::howmany() const
   {
      return this->mHowmany;
   }

   int FftSetup::specSize() const
   {
      return this->mSpecSize;
   }

   int FftSetup::padSize() const
   {
      return this->mBwdSize - this->mSpecSize;
   }

   MHDFloat FftSetup::scale() const
   {
      // Assert to make sure scale is initialised
      assert(this->mScale > 0.0);

      return this->mScale;
   }

}
}
