/** \file TFTScheme.cpp
 *  \brief Source of the Chebyshev(FFT) + Fourier + Chebyshev(FFT) scheme implementation
 */

// System includes
//
#include <set>

// External includes
//

// Class include
//
#include "SpatialSchemes/3D/TFTScheme.hpp"

// Project includes
//

namespace GeoMHDiSCC {

   const int TFTScheme::DIMENSIONS = 3;

   Transform::SharedFftSetup TFTScheme::spSetup1D(SharedResolution spRes)
   {
      // Get size of FFT transform
      int size = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DATF1D>();

      // Get spectral size of the FFT
      int specSize = spRes->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

      // Get number of transforms
      int howmany = 0;
      for(int i = 0; i < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); i++)
      {
         howmany += spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(i);
      }

      return Transform::SharedFftSetup(new Transform::FftSetup(size, howmany, specSize, false));
   }

   Transform::SharedFftSetup TFTScheme::spSetup2D(SharedResolution spRes)
   {
      // Get size of FFT transform
      int size = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DATF1D>();

      // Get spectral size of the FFT
      int specSize = spRes->sim()->dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL);

      // Get number of transforms
      int howmany = 0;
      for(int i = 0; i < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); i++)
      {
         howmany += spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(i);
      }

      return Transform::SharedFftSetup(new Transform::FftSetup(size, howmany, specSize, true));
   }

   Transform::SharedFftSetup TFTScheme::spSetup3D(SharedResolution spRes)
   {
      // Get size of FFT transform
      int size = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DATF1D>();

      // Get spectral size of the FFT
      int specSize = spRes->sim()->dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL);

      // Get number of transforms
      int howmany = 0;
      for(int i = 0; i < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); i++)
      {
         howmany += spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(i);
      }

      return Transform::SharedFftSetup(new Transform::FftSetup(size, howmany, specSize, false));
   }

   TFTScheme::TFTScheme(const ArrayI& dim)
      : Regular3DScheme(dim)
   {
   }

   TFTScheme::~TFTScheme()
   {
   }

   void TFTScheme::setDimensions()
   {
      //
      // Initialise first transform
      //

      // Initialise forward dimension of first transform
      this->mDimensions.at(0)(0) = this->nX();

      // Initialise backward dimension of first transform
      this->mDimensions.at(0)(1) = this->nI();

      // Initialise second dimension of first transform
      this->mDimensions.at(0)(2) = this->nK();

      // Initialise third dimension of first transform
      this->mDimensions.at(0)(3) = this->nJ();

      //
      // Initialise second transform
      //

      // Initialise forward dimension of second transform
      this->mDimensions.at(1)(0) = this->nY();

      // Initialise backward dimension of second transform
      this->mDimensions.at(1)(1) = this->nY()/2 + 1;

      // Initialise second dimension of second transform
      this->mDimensions.at(1)(2) = this->nX();

      // Initialise third dimension of second transform
      this->mDimensions.at(1)(3) = this->nK();

      //
      // Initialise third transform
      //

      // Initialise forward dimension of third transform
      this->mDimensions.at(2)(0) = this->nZ();

      // Initialise backward dimension of third transform
      this->mDimensions.at(2)(1) = this->nZ();

      // Initialise second dimension of third transform
      this->mDimensions.at(2)(2) = this->nY();

      // Initialise third dimension of third transform
      this->mDimensions.at(2)(3) = this->nX();
   }

   void TFTScheme::setCosts(const int shift)
   {
      // Set first costs
      this->mCosts(0) = 1.0;

      // Set second costs
      this->mCosts(1) = 1.0;

      // Set third costs
      this->mCosts(2) = 1.0;
   }

   void TFTScheme::setScalings(const int shift)
   {
      // Set first scaling
      this->mScalings(0) = 1.0;

      // Set second scaling
      this->mScalings(1) = 1.0;

      // Set third scaling
      this->mScalings(2) = 1.0;
   }

   void TFTScheme::setMemory(const int shift)
   {
      // Set first memory footprint
      this->mMemory(0) = 1.0;

      // Set second memory footprint
      this->mMemory(1) = 1.0;

      // Set third memory footprint
      this->mMemory(2) = 1.0;
   }

   bool TFTScheme::applicable() const
   {
      return true;
   }

   Array TFTScheme::loadWeights()
   {
      Array weights(this->DIMENSIONS);

      weights = this->mCosts.array() * this->mScalings.array();

      return weights;
   }

   double TFTScheme::memoryScore(SharedResolution spRes)
   {
      #ifdef GEOMHDISCC_MEMORYUSAGE_LIMITED
         return this->mMemory.prod();
      #else
         return 1.0;
      #endif //GEOMHDISCC_MEMORYUSAGE_LIMITED
   }
}
