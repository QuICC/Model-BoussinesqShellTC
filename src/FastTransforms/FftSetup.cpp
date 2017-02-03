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

namespace QuICC {

namespace Transform {

   FftSetup::FftSetup(const int size, const int howmany, const int specSize, const FftSetup::Type type)
      : TransformSetup(size, howmany, specSize), mBwdSize(0), mType(type), mScale(-1), mBoxScale(-1), mEvenBlocks(0,0), mOddBlocks(0,0)
   {
      this->setBackwardSize();

      // Safety assert
      assert(this->mBwdSize > this->mSpecSize);
      assert(this->mFwdSize >= this->mBwdSize);
   }

   FftSetup::FftSetup(const int size, const ArrayI& howmany, const MatrixI& evenBlocks, const MatrixI& oddBlocks, const int specSize, const FftSetup::Type type)
      : TransformSetup(size, howmany, specSize), mBwdSize(0), mType(type), mScale(-1), mBoxScale(-1), mEvenBlocks(evenBlocks), mOddBlocks(oddBlocks)
   {
      this->setBackwardSize();

      // Safety assert
      assert(this->mBwdSize > this->mSpecSize);
      assert(this->mFwdSize >= this->mBwdSize);
   }

   FftSetup::FftSetup(const int size, const int howmany, const MatrixI& idBlocks, const int specSize, const FftSetup::Type type)
      : TransformSetup(size, howmany, specSize), mBwdSize(0), mType(type), mScale(-1), mBoxScale(-1), mIdBlocks(idBlocks)
   {
      this->setBackwardSize();

      // Safety assert
      assert(this->mBwdSize > this->mSpecSize);
      assert(this->mFwdSize >= this->mBwdSize);
   }

   FftSetup::~FftSetup()
   {
   }

   void FftSetup::setBackwardSize()
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
   }

   void FftSetup::setScale(const MHDFloat scale)
   {
      this->mScale = scale;
   }

   void FftSetup::setBoxScale(const MHDFloat boxScale)
   {
      this->mBoxScale = boxScale;
   }

   FftSetup::Type FftSetup::type() const
   {
      return this->mType;
   }

   int FftSetup::bwdSize() const
   {
      return this->mBwdSize;
   }

   const MatrixI& FftSetup::evenBlocks() const
   {
      return this->mEvenBlocks;
   }

   const MatrixI& FftSetup::oddBlocks() const
   {
      return this->mOddBlocks;
   }

   const MatrixI& FftSetup::idBlocks() const
   {
      return this->mIdBlocks;
   }

   int FftSetup::padSize() const
   {
      if(this->mType == FftSetup::COMPLEX)
      {
         return this->mBwdSize - this->mSpecSize;
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

   MHDFloat FftSetup::boxScale() const
   {
      // Assert to make sure scale is initialised
      assert(this->mBoxScale > 0.0);

      return this->mBoxScale;
   }

}
}
