/** \file FFFScheme.cpp
 *  \brief Source of the Fourier + Fourier + Fourier scheme implementation
 */

// System includes
//
#include <set>

// External includes
//

// Class include
//
#include "SpatialSchemes/3D/FFFScheme.hpp"

// Project includes
//

namespace GeoMHDiSCC {
   
   std::string FFFScheme::type()
   {
      return "FFF";
   }

   void FFFScheme::addTransformSetups(SharedResolution spRes) const
   {
      // Add setup for first transform
      Transform::SharedFftSetup  spS1D = this->spSetup1D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA1D, spS1D);

      // Add setup for second transform
      Transform::SharedFftSetup  spS2D = this->spSetup2D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA2D, spS2D);

      // Add setup for third transform
      Transform::SharedFftSetup  spS3D = this->spSetup3D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA3D, spS3D);
   }

   Transform::SharedFftSetup FFFScheme::spSetup1D(SharedResolution spRes) const
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

      return Transform::SharedFftSetup(new Transform::FftSetup(size, howmany, specSize, Transform::FftSetup::EQUAL));
   }

   Transform::SharedFftSetup FFFScheme::spSetup2D(SharedResolution spRes) const
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

      return Transform::SharedFftSetup(new Transform::FftSetup(size, howmany, specSize, Transform::FftSetup::EQUAL));
   }

   Transform::SharedFftSetup FFFScheme::spSetup3D(SharedResolution spRes) const
   {
      // Get size of FFT transform
      int size = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DATF1D>();

      // Get spectral size of the FFT
      int specSize = spRes->sim()->dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL);

      // Get number of transforms
      int howmany = 0;
      for(int i = 0; i < spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>(); i++)
      {
         howmany += spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(i);
      }

      return Transform::SharedFftSetup(new Transform::FftSetup(size, howmany, specSize, Transform::FftSetup::MIXED));
   }

   FFFScheme::FFFScheme(const ArrayI& dim)
      : IRegular3DScheme(dim)
   {
   }

   FFFScheme::~FFFScheme()
   {
   }

   void FFFScheme::setDimensions()
   {
      //
      // Compute sizes
      //

      // Get standard dealiased FFT size
      int nX = Transform::FftwTools::dealiasFft(this->mI+1);
      // Check for optimised FFT sizes
      nX = Transform::FftwTools::optimizeFft(nX);

      // Get standard dealiased FFT size
      int nY = Transform::FftwTools::dealiasFft(this->mJ+1);
      // Check for optimised FFT sizes
      nY = Transform::FftwTools::optimizeFft(nY);

      // Get standard dealiased FFT size
      int nZ = Transform::FftwTools::dealiasMixedFft(this->mK+1);
      // Check for optimised FFT sizes
      nZ = Transform::FftwTools::optimizeFft(nZ);

      //
      // Initialise first transform
      //

      // Initialise forward dimension of first transform
      this->setDimension(nX, Dimensions::Transform::TRA1D, Dimensions::Data::DATF1D);

      // Initialise backward dimension of first transform
      this->setDimension(nX, Dimensions::Transform::TRA1D, Dimensions::Data::DATB1D);

      // Initialise second dimension of first transform
      this->setDimension(this->mK + 1, Dimensions::Transform::TRA1D, Dimensions::Data::DAT2D);

      // Initialise third dimension of first transform
      this->setDimension(this->mJ + 1, Dimensions::Transform::TRA1D, Dimensions::Data::DAT3D);

      //
      // Initialise second transform
      //

      // Initialise forward dimension of second transform
      this->setDimension(nY, Dimensions::Transform::TRA2D, Dimensions::Data::DATF1D);

      // Initialise backward dimension of second transform
      this->setDimension(nY, Dimensions::Transform::TRA2D, Dimensions::Data::DATB1D);

      // Initialise second dimension of second transform
      this->setDimension(nX, Dimensions::Transform::TRA2D, Dimensions::Data::DAT2D);

      // Initialise third dimension of second transform
      this->setDimension(this->mK + 1, Dimensions::Transform::TRA2D, Dimensions::Data::DAT3D);

      //
      // Initialise third transform
      //

      // Initialise forward dimension of third transform
      this->setDimension(nZ, Dimensions::Transform::TRA3D, Dimensions::Data::DATF1D);

      // Initialise backward dimension of third transform
      this->setDimension(nZ/2 + 1, Dimensions::Transform::TRA3D, Dimensions::Data::DATB1D);

      // Initialise second dimension of third transform
      this->setDimension(nY, Dimensions::Transform::TRA3D, Dimensions::Data::DAT2D);

      // Initialise third dimension of third transform
      this->setDimension(nX, Dimensions::Transform::TRA3D, Dimensions::Data::DAT3D);
   }

   void FFFScheme::setCosts()
   {
      // Set first transform cost
      this->setCost(1.0, Dimensions::Transform::TRA1D);

      // Set second transform cost
      this->setCost(1.0, Dimensions::Transform::TRA2D);

      // Set third transform cost
      this->setCost(1.0, Dimensions::Transform::TRA3D);
   }

   void FFFScheme::setScalings()
   {
      // Set first transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA1D);

      // Set second transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA2D);

      // Set third transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA3D);
   }

   void FFFScheme::setMemoryScore()
   {
      // Set first transform memory footprint
      this->setMemory(1.0, Dimensions::Transform::TRA1D);

      // Set second transform memory footprint
      this->setMemory(1.0, Dimensions::Transform::TRA2D);

      // Set third transform memory footprint
      this->setMemory(1.0, Dimensions::Transform::TRA3D);
   }

   bool FFFScheme::applicable() const
   {
      return true;
   }
}
