/** \file SLFScheme.cpp
 *  \brief Source of the spherical Chebyshev(FFT) + Spherical Harmonics (Associated Legendre(poly) + Fourrier) scheme implementation
 */

// System includes
//
#include <set>

// External includes
//

// Class include
//
#include "SpatialSchemes/3D/SLFScheme.hpp"

// Project includes
//
#include "PolynomialTransforms/PolynomialTools.hpp"

namespace GeoMHDiSCC {

namespace Schemes {
   
   std::string SLFScheme::type()
   {
      return "SLF";
   }

   void SLFScheme::addTransformSetups(SharedResolution spRes) const
   {
      // Add setup for first transform
      Transform::SharedFftSetup  spS1D = this->spSetup1D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA1D, spS1D);

      // Add setup for second transform
      Transform::SharedPolySetup  spS2D = this->spSetup2D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA2D, spS2D);

      // Add setup for third transform
      Transform::SharedFftSetup  spS3D = this->spSetup3D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA3D, spS3D);
   }

   Transform::SharedFftSetup SLFScheme::spSetup1D(SharedResolution spRes) const
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

   Transform::SharedPolySetup SLFScheme::spSetup2D(SharedResolution spRes) const
   {
      // Get physical size of polynomial transform
      int size = spRes->cpu()->dim(Dimensions::Transform::TRA2D)->dim<Dimensions::Data::DATF1D>();

      // Get spectral size of the polynomial transform
      int specSize = spRes->sim()->dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL);

      // Storage for the list of indexes
      std::vector<ArrayI>  fast;
      fast.reserve(spRes->cpu()->dim(Dimensions::Transform::TRA2D)->dim<Dimensions::Data::DAT3D>());
      ArrayI  slow(spRes->cpu()->dim(Dimensions::Transform::TRA2D)->dim<Dimensions::Data::DAT3D>());

      // Multiplier from second dimension 
      ArrayI mult(spRes->cpu()->dim(Dimensions::Transform::TRA2D)->dim<Dimensions::Data::DAT3D>());

      // Get number of transforms and list of indexes
      int howmany = 0;
      for(int i = 0; i < spRes->cpu()->dim(Dimensions::Transform::TRA2D)->dim<Dimensions::Data::DAT3D>(); i++)
      {
         howmany += spRes->cpu()->dim(Dimensions::Transform::TRA2D)->dim<Dimensions::Data::DAT2D>(i);

         slow(i) = spRes->cpu()->dim(Dimensions::Transform::TRA2D)->idx<Dimensions::Data::DAT3D>(i);

         fast.push_back(ArrayI(spRes->cpu()->dim(Dimensions::Transform::TRA2D)->dim<Dimensions::Data::DATB1D>(i)));
         for(int j = 0; j < fast.at(i).size(); j++)
         {
            fast.at(i)(j) = spRes->cpu()->dim(Dimensions::Transform::TRA2D)->idx<Dimensions::Data::DATB1D>(j,i);
         }

         mult(i) = spRes->cpu()->dim(Dimensions::Transform::TRA2D)->dim<Dimensions::Data::DAT2D>(i);
      }

      return Transform::SharedPolySetup(new Transform::PolySetup(size, howmany, specSize, fast, slow, mult));
   }

   Transform::SharedFftSetup SLFScheme::spSetup3D(SharedResolution spRes) const
   {
      // Get size of FFT transform
      int size = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DATF1D>();

      // Get spectral size of the FFT
      int specSize = spRes->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      // Get number of transforms
      int howmany = 0;
      for(int i = 0; i < spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>(); i++)
      {
         howmany += spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(i);
      }

      return Transform::SharedFftSetup(new Transform::FftSetup(size, howmany, specSize, Transform::FftSetup::MIXED));
   }

   SLFScheme::SLFScheme(const ArrayI& dim)
      : IRegularSHScheme(dim)
   {
   }

   SLFScheme::~SLFScheme()
   {
   }

   void SLFScheme::setDimensions()
   {
      //
      // Compute sizes
      //

      // Get standard dealiased FFT size
      int nR = Transform::FftwTools::dealiasFft(this->mI+1);
      // Check for optimised FFT sizes
      nR = Transform::FftwTools::optimizeFft(nR);

      // Get dealiased associated Legendre transform size
      int nTh = Transform::PolynomialTools::dealias(this->mL+1);

      // Get standard dealiased FFT size
      int nPh = Transform::FftwTools::dealiasMixedFft(this->mM+1);
      // Check for optimised FFT sizes
      nPh = Transform::FftwTools::optimizeFft(nPh);

      //
      // Initialise first transform
      //

      // Initialise forward dimension of first transform
      this->setDimension(nR, Dimensions::Transform::TRA1D, Dimensions::Data::DATF1D);

      // Initialise backward dimension of first transform
      this->setDimension(nR, Dimensions::Transform::TRA1D, Dimensions::Data::DATB1D);

      // Initialise second dimension of first transform
      this->setDimension(this->mL + 1, Dimensions::Transform::TRA1D, Dimensions::Data::DAT2D);

      // Initialise third dimension of first transform
      this->setDimension(this->mM + 1, Dimensions::Transform::TRA1D, Dimensions::Data::DAT3D);

      //
      // Initialise second transform
      //

      // Initialise forward dimension of second transform
      this->setDimension(nTh, Dimensions::Transform::TRA2D, Dimensions::Data::DATF1D);

      // Initialise backward dimension of second transform
      this->setDimension(this->mL + 1, Dimensions::Transform::TRA2D, Dimensions::Data::DATB1D);

      // Initialise second dimension of second transform
      this->setDimension(nR, Dimensions::Transform::TRA2D, Dimensions::Data::DAT2D);

      // Initialise third dimension of second transform
      this->setDimension(this->mM + 1, Dimensions::Transform::TRA2D, Dimensions::Data::DAT3D);

      //
      // Initialise third transform
      //

      // Initialise forward dimension of third transform
      this->setDimension(nPh, Dimensions::Transform::TRA3D, Dimensions::Data::DATF1D);

      // Initialise backward dimension of third transform
      this->setDimension(nPh/2 + 1, Dimensions::Transform::TRA3D, Dimensions::Data::DATB1D);

      // Initialise second dimension of third transform
      this->setDimension(nTh, Dimensions::Transform::TRA3D, Dimensions::Data::DAT2D);

      // Initialise third dimension of third transform
      this->setDimension(nR, Dimensions::Transform::TRA3D, Dimensions::Data::DAT3D);
   }

   void SLFScheme::setCosts()
   {
      // Set first transform cost
      this->setCost(1.0, Dimensions::Transform::TRA1D);

      // Set second transform cost
      this->setCost(1.0, Dimensions::Transform::TRA2D);

      // Set third transform cost
      this->setCost(1.0, Dimensions::Transform::TRA3D);
   }

   void SLFScheme::setScalings()
   {
      // Set first transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA1D);

      // Set second transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA2D);

      // Set third transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA3D);
   }

   void SLFScheme::setMemoryScore()
   {
      // Set first transform memory footprint
      this->setMemory(1.0, Dimensions::Transform::TRA1D);

      // Set second transform memory footprint
      this->setMemory(1.0, Dimensions::Transform::TRA2D);

      // Set third transform memory footprint
      this->setMemory(1.0, Dimensions::Transform::TRA3D);
   }

   bool SLFScheme::applicable() const
   {
      return true;
   }
}
}