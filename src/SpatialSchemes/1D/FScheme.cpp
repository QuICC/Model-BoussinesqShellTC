/** 
 * @file FScheme.cpp
 * @brief Source of the Fourier scheme implementation
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// External includes
//

// Class include
//
#include "SpatialSchemes/1D/FScheme.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"
#include "FastTransforms/FftwTools.hpp"

namespace GeoMHDiSCC {

namespace Schemes {
   
   std::string FScheme::type()
   {
      return "F";
   }

   void FScheme::addTransformSetups(SharedResolution spRes) const
   {
      // Add setup for first transform
      Transform::SharedFftSetup  spS1D = this->spSetup1D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA1D, spS1D);
   }

   Transform::SharedFftSetup FScheme::spSetup1D(SharedResolution spRes) const
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

      return Transform::SharedFftSetup(new Transform::FftSetup(size, howmany, specSize, Transform::FftSetup::MIXED));
   }

   FScheme::FScheme(const ArrayI& dim)
      : IRegular1DScheme(dim)
   {
   }

   FScheme::~FScheme()
   {
   }

   void FScheme::setDimensions()
   {
      //
      // Compute sizes
      //

      // Get standard dealiased FFT size
      int nX = Transform::FftwTools::dealiasMixedFft(this->mI+1);
      // Check for optimised FFT sizes
      nX = Transform::FftwTools::optimizeFft(nX);

      //
      // Initialise first transform
      //

      // Initialise forward dimension of first transform
      this->setDimension(nX, Dimensions::Transform::TRA1D, Dimensions::Data::DATF1D);

      // Initialise backward dimension of first transform
      this->setDimension(nX/2 + 1, Dimensions::Transform::TRA1D, Dimensions::Data::DATB1D);
   }

   void FScheme::setCosts()
   {
      // Set first transform cost
      this->setCost(1.0, Dimensions::Transform::TRA1D);
   }

   void FScheme::setScalings()
   {
      // Set first transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA1D);
   }

   void FScheme::setMemoryScore()
   {
      // Set first transform memory footprint
      this->setMemory(1.0, Dimensions::Transform::TRA1D);
   }

   bool FScheme::applicable() const
   {
      return true;
   }
}
}
