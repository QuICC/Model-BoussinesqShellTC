/** 
 * @file FftSetup.cpp
 * @brief Source of FFT setup class
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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
      : TransformSetup(size, howmany, specSize), mBwdSize(0), mType(type), mScale(-1)
   {
      // Set the backward size
      if(this->mType == FftSetup::MIXED)
      {
         this->mBwdSize = this->mFwdSize/2 + 1;
      } else if(this->mType == FftSetup::REAL)
      {
         this->mBwdSize = this->mFwdSize;
      } else if(this->mType == FftSetup::COMPLEX)
      {
         this->mBwdSize = this->mFwdSize;
      } else if(this->mType == FftSetup::COMPONENT)
      {
         this->mBwdSize = this->mFwdSize;
      } else
      {
         throw Exception("Unknown FFT setup type requested");
      }

      // Safety assert
      assert(this->mBwdSize > this->mSpecSize);
      assert(this->mFwdSize >= this->mBwdSize);
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

   int FftSetup::bwdSize() const
   {
      return this->mBwdSize;
   }

   int FftSetup::padSize() const
   {
      if(this->mType == FftSetup::COMPLEX)
      {
         return this->mBwdSize - 2*this->mSpecSize;
      } else
      {
         return this->mBwdSize - this->mSpecSize;
      }
   }

   MHDFloat FftSetup::scale() const
   {
      // Assert to make sure scale is initialised
      assert(this->mScale > 0.0);

      return this->mScale;
   }

}
}
