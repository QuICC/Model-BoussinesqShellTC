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
#include "Exceptions/Exception.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   FftSetup::FftSetup(const int size, const int howmany, const int specSize, const FftSetup::Type type)
      : mFwdSize(size), mBwdSize(0), mHowmany(howmany), mSpecSize(specSize), mType(type), mScale(-1)
   {
      // Set the backward size
      if(this->mType == FftSetup::MIXED)
      {
         this->mBwdSize = this->mFwdSize/2 + 1;
      } else if(this->mType == FftSetup::EQUAL)
      {
         this->mBwdSize = this->mFwdSize;
      } else if(this->mType == FftSetup::COMPONENT)
      {
         this->mBwdSize = this->mFwdSize;
      } else
      {
         throw Exception("Unknown FFT setup type requested");
      }
   }

   FftSetup::~FftSetup()
   {
   }

   void FftSetup::setScale(const MHDFloat scale)
   {
      this->mScale = scale;
   }

   FftSetup::Type FftSetup::type() const
   {
      return this->mType;
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
