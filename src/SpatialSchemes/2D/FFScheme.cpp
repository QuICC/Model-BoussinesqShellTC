/** \file FFScheme.cpp
 *  \brief Source of the Fourier + Fourier scheme implementation
 */

// System includes
//

// External includes
//

// Class include
//
#include "Base/SpatialSchemes/2D/FFScheme.hpp"

// Project includes
//

namespace GeoMHDiSCC {

   const int FFScheme::DIMENSIONS = 2;

   SharedFFTSetup FFScheme::spSetup1D(SharedResolution spRes)
   {
      // Get size of FFT transform
      int size = spRes->cpu()->dim(0)->dimFwd();

      // Get spectral size of the FFT
      int specSize = spRes->sim()->dim(Dimensions::Space::SPECTRAL,0);

      // Get number of transforms
      int howmany = 0;
      for(int i = 0; i < spRes->cpu()->dim(0)->dim3D(); i++)
      {
         howmany += spRes->cpu()->dim(0)->dim2D(i);
      }

      return SharedFFTSetup(new FFTSetup(size, howmany, specSize));
   }

   SharedFFTSetup FFScheme::spSetup2D(SharedResolution spRes)
   {
      // Get size of FFT transform
      int size = spRes->cpu()->dim(1)->dimFwd();

      // Get spectral size of the FFT
      int specSize = spRes->sim()->dim(Dimensions::Space::SPECTRAL,1);

      // Get number of transforms
      int howmany = 0;
      for(int i = 0; i < spRes->cpu()->dim(1)->dim3D(); i++)
      {
         howmany += spRes->cpu()->dim(1)->dim2D(i);
      }

      return SharedFFTSetup(new FFTSetup(size, howmany, specSize));
   }

   FFScheme::FFScheme(const ArrayI& dim)
      : Regular2DScheme(dim)
   {
   }

   void FFScheme::setDimensions()
   {
      //
      // Initialise first transform
      //

      // Initialise forward dimension of first transform
      this->mDimensions.at(0)(0) = this->nX();

      // Initialise backward dimension of first transform
      this->mDimensions.at(0)(1) = this->nI();

      // Initialise second dimension of first transform
      this->mDimensions.at(0)(2) = this->nI();

      //
      // Initialise second transform
      //

      // Initialise forward dimension of second transform
      this->mDimensions.at(1)(0) = this->nY();

      // Initialise backward dimension of first transform
      this->mDimensions.at(1)(1) = this->nI();

      // Initialise second dimension of first transform
      this->mDimensions.at(1)(2) = this->nX();
   }

   void FFScheme::setCosts(const int shift)
   {
      // Set first costs
      this->mCosts(0) = 1.0;

      // Set second costs
      this->mCosts(1) = 1.0;
   }

   void FFScheme::setScalings(const int shift)
   {
      // Set first scaling
      this->mScalings(0) = 1.0;

      // Set second scaling
      this->mScalings(1) = 1.0;
   }

   void FFScheme::setMemory(const int shift)
   {
      // Set first memory footprint
      this->mMemory(0) = 1.0;

      // Set second memory footprint
      this->mMemory(1) = 1.0;
   }

   bool FFScheme::applicable() const
   {
      return true;
   }

   Array FFScheme::loadWeights()
   {
      Array weights(this->DIMENSIONS);

      weights = this->mCosts.array() * this->mScalings.array();

      return weights;
   }

   double FFScheme::memoryScore(SharedResolution spRes)
   {
      #ifdef GEOMHDISCC_MEMORYUSAGE_LIMITED
         return this->mMemory.prod();
      #else
         return 1.0;
      #endif //GEOMHDISCC_MEMORYUSAGE_LIMITED
   }

}
