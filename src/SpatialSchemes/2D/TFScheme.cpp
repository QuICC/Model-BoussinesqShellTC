/** 
 * @file TFScheme.cpp
 * @brief Source of the Chebyshev(FFT) + Fourier scheme implementation
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// External includes
//

// Class include
//
#include "SpatialSchemes/2D/TFScheme.hpp"

// Project includes
//
#include "TypeSelectors/FftSelector.hpp"

namespace GeoMHDiSCC {

namespace Schemes {
   
   std::string TFScheme::type()
   {
      return "TF";
   }

   void TFScheme::addTransformSetups(SharedResolution spRes) const
   {
      // Add setup for first transform
      Transform::SharedFftSetup  spS1D = this->spSetup1D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA1D, spS1D);

      // Add setup for second transform
      Transform::SharedFftSetup  spS2D = this->spSetup2D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA2D, spS2D);
   }

   Transform::SharedFftSetup TFScheme::spSetup1D(SharedResolution spRes) const
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

      return Transform::SharedFftSetup(new Transform::FftSetup(size, howmany, specSize, Transform::FftSetup::COMPONENT));
   }

   Transform::SharedFftSetup TFScheme::spSetup2D(SharedResolution spRes) const
   {
      // Get size of FFT transform
      int size = spRes->cpu()->dim(Dimensions::Transform::TRA2D)->dim<Dimensions::Data::DATF1D>();

      // Get spectral size of the FFT
      int specSize = spRes->sim()->dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL);

      // Get number of transforms
      int howmany = 0;
      for(int i = 0; i < spRes->cpu()->dim(Dimensions::Transform::TRA2D)->dim<Dimensions::Data::DAT3D>(); i++)
      {
         howmany += spRes->cpu()->dim(Dimensions::Transform::TRA2D)->dim<Dimensions::Data::DAT2D>(i);
      }

      return Transform::SharedFftSetup(new Transform::FftSetup(size, howmany, specSize, Transform::FftSetup::MIXED));
   }

   TFScheme::TFScheme(const ArrayI& dim)
      : IRegular2DScheme(dim)
   {
   }

   TFScheme::~TFScheme()
   {
   }

   void TFScheme::setDimensions()
   {
      //
      // Set transform space sizes
      //
      ArrayI traSize(2);
      traSize(0) = this->mI + 1;
      traSize(1) = this->mJ + 1;
      this->setTransformSpace(traSize);

      //
      // Compute sizes
      //

      // Get standard dealiased FFT size
      int nX = Transform::FftToolsType::dealiasCosFft(this->mI+1);
      // Check for optimised FFT sizes
      nX = Transform::FftToolsType::optimizeFft(nX);

      // Get standard dealiased FFT size
      int nY = Transform::FftToolsType::dealiasMixedFft(this->mJ+1);
      // Check for optimised FFT sizes
      nY = Transform::FftToolsType::optimizeFft(nY);

      //
      // Initialise first transform
      //

      // Initialise forward dimension of first transform
      this->setDimension(nX, Dimensions::Transform::TRA1D, Dimensions::Data::DATF1D);

      // Initialise backward dimension of first transform
      this->setDimension(nX, Dimensions::Transform::TRA1D, Dimensions::Data::DATB1D);

      // Initialise second dimension of first transform
      this->setDimension(traSize(1), Dimensions::Transform::TRA1D, Dimensions::Data::DAT2D);

      //
      // Initialise second transform
      //

      // Initialise forward dimension of second transform
      this->setDimension(nY, Dimensions::Transform::TRA2D, Dimensions::Data::DATF1D);

      // Initialise backward dimension of second transform
      this->setDimension(nY/2 + 1, Dimensions::Transform::TRA2D, Dimensions::Data::DATB1D);

      // Initialise second dimension of second transform
      this->setDimension(nX, Dimensions::Transform::TRA2D, Dimensions::Data::DAT2D);
   }

   void TFScheme::setCosts()
   {
      // Set first transform cost
      this->setCost(1.0, Dimensions::Transform::TRA1D);

      // Set second transform cost
      this->setCost(1.0, Dimensions::Transform::TRA2D);
   }

   void TFScheme::setScalings()
   {
      // Set first transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA1D);

      // Set second transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA2D);
   }

   void TFScheme::setMemoryScore()
   {
      // Set first transform memory footprint
      this->setMemory(1.0, Dimensions::Transform::TRA1D);

      // Set second transform memory footprint
      this->setMemory(1.0, Dimensions::Transform::TRA2D);
   }

   bool TFScheme::applicable() const
   {
      return true;
   }

}
}
