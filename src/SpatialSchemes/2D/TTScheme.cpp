/** 
 * @file TTScheme.cpp
 * @brief Source of the Chebyshev(FFT) + Chebyshev(FFT) scheme implementation
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include "Framework/FrameworkMacro.h"

// External includes
//
#include <set>
#include <vector>

// Class include
//
#include "SpatialSchemes/2D/TTScheme.hpp"

// Project includes
//
#include "Resolutions/Tools/RegularIndexCounter.hpp"

namespace GeoMHDiSCC {

namespace Schemes {
   
   std::string TTScheme::type()
   {
      return "TT";
   }

   void TTScheme::tuneResolution(SharedResolution spRes, const Parallel::SplittingDescription& descr)
   {
      TTScheme::tuneMpiResolution(descr);
      
      // Create spectral space sub communicators
      #if defined QUICC_MPI && defined QUICC_MPISPSOLVE
         // Initialise the ranks with local rank
         std::set<int>  ranks;
         for(int cpu = 0; cpu < spRes->nCpu(); ++cpu)
         {
            ranks.insert(cpu);
         }

         FrameworkMacro::initSubComm(FrameworkMacro::SPECTRAL, 1);

         FrameworkMacro::setSubComm(FrameworkMacro::SPECTRAL, 0, ranks);
      #endif //defined QUICC_MPI && defined QUICC_MPISPSOLVE
   }

   void TTScheme::addTransformSetups(SharedResolution spRes) const
   {
      // Add setup for first transform
      Transform::SharedFftSetup  spS1D = this->spSetup1D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA1D, spS1D);

      // Add setup for second transform
      Transform::SharedFftSetup  spS2D = this->spSetup2D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA2D, spS2D);
   }

   Transform::SharedFftSetup TTScheme::spSetup1D(SharedResolution spRes) const
   {
      // Get size of FFT transform
      int size = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DATF1D>();

      // Get spectral size of the FFT
      int specSize = spRes->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

      // Get number of transforms
      int howmany = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>();

      return Transform::SharedFftSetup(new Transform::FftSetup(size, howmany, specSize, Transform::FftSetup::REAL));
   }

   Transform::SharedFftSetup TTScheme::spSetup2D(SharedResolution spRes) const
   {
      // Get size of FFT transform
      int size = spRes->cpu()->dim(Dimensions::Transform::TRA2D)->dim<Dimensions::Data::DATF1D>();

      // Get spectral size of the FFT
      int specSize = spRes->sim()->dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL);

      // Get number of transforms
      int howmany = spRes->cpu()->dim(Dimensions::Transform::TRA2D)->dim<Dimensions::Data::DAT2D>();

      return Transform::SharedFftSetup(new Transform::FftSetup(size, howmany, specSize, Transform::FftSetup::REAL));
   }

   TTScheme::TTScheme(const ArrayI& dim)
      : IRegular2DScheme(dim)
   {
   }

   TTScheme::~TTScheme()
   {
   }

   void TTScheme::setDimensions()
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
      int nX = Transform::Fft::ToolsSelector::dealiasCosFft(this->mI+1);
      // Check for optimised FFT sizes
      nX = Transform::Fft::ToolsSelector::optimizeFft(nX);

      // Get standard dealiased FFT size
      int nY = Transform::Fft::ToolsSelector::dealiasCosFft(this->mJ+1);
      // Check for optimised FFT sizes
      nY = Transform::Fft::ToolsSelector::optimizeFft(nY);

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
      this->setDimension(nY, Dimensions::Transform::TRA2D, Dimensions::Data::DATB1D);

      // Initialise second dimension of second transform
      this->setDimension(nX, Dimensions::Transform::TRA2D, Dimensions::Data::DAT2D);
   }

   void TTScheme::setCosts()
   {
      // Set first transform cost
      this->setCost(1.0, Dimensions::Transform::TRA1D);

      // Set second transform cost
      this->setCost(1.0, Dimensions::Transform::TRA2D);
   }

   void TTScheme::setScalings()
   {
      // Set first transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA1D);

      // Set second transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA2D);
   }

   void TTScheme::setMemoryScore()
   {
      // Set first transform memory footprint
      this->setMemory(1.0, Dimensions::Transform::TRA1D);

      // Set second transform memory footprint
      this->setMemory(1.0, Dimensions::Transform::TRA2D);
   }

   bool TTScheme::applicable() const
   {
      return true;
   }

}
}
