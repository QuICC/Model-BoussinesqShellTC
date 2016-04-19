/** 
 * @file SLFmScheme.cpp
 * @brief Source of the spherical Chebyshev(FFT) + Spherical Harmonics (Associated Legendre(poly) + Fourrier) scheme implementation with m spectral ordering
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//
#include "Framework/FrameworkMacro.h"

// System includes
//
#include <set>

// External includes
//

// Class include
//
#include "SpatialSchemes/3D/SLFmScheme.hpp"

// Project includes
//
#include "PolynomialTransforms/PolynomialTools.hpp"

namespace GeoMHDiSCC {

namespace Schemes {
   
   std::string SLFmScheme::type()
   {
      return "SLFm";
   }

   void SLFmScheme::tuneResolution(SharedResolution spRes, const Parallel::SplittingDescription& descr)
   {
      SLFmScheme::tuneMpiResolution(descr);
      
      // Create spectral space sub communicators
      #if defined GEOMHDISCC_MPI && defined GEOMHDISCC_MPISPSOLVE

         // MPI error code
         int ierr;

         // Extract modes available on local rank
         std::vector<std::set<int> >  ranks;
         ArrayI modes(spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>());
         std::map<int, int>  mapModes;
         int k_ = 0;
         for(int k = 0; k < spRes->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL); ++k)
         {
            ranks.push_back(std::set<int>());
            if(k_ < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>() && k == spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k_))
            {
               ranks.back().insert(FrameworkMacro::id());
               modes(k_) = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k_);
               mapModes.insert(std::make_pair(spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k_),k));
               k_++;
            }
         }

         // Share modes information with CPUs in first transform group
         int traId = 0;
         ArrayI tmp;
         for(int commCpu = 0; commCpu < FrameworkMacro::transformCpus(0).size(); ++commCpu)
         {
            ArrayI iData(2);
            iData.setConstant(-1);
            if(commCpu == FrameworkMacro::transformId(traId))
            {
               // Send the size and global rank
               iData(0) = modes.size();
               iData(1) = FrameworkMacro::id();
               FrameworkMacro::syncTransform(traId);
               ierr = MPI_Bcast(iData.data(), iData.size(), MPI_INT, commCpu, FrameworkMacro::transformComm(traId));
               FrameworkMacro::check(ierr, 815);

               // Send modes
               FrameworkMacro::syncTransform(traId);
               ierr = MPI_Bcast(modes.data(), modes.size(), MPI_INT, commCpu, FrameworkMacro::transformComm(traId));
               FrameworkMacro::check(ierr, 817);
            } else
            {
               // Get size and global rank
               FrameworkMacro::syncTransform(traId);
               ierr = MPI_Bcast(iData.data(), iData.size(), MPI_INT, commCpu, FrameworkMacro::transformComm(traId));
               FrameworkMacro::check(ierr, 818);

               // Receive modes
               tmp.resize(iData(0));
               FrameworkMacro::syncTransform(traId);
               ierr = MPI_Bcast(tmp.data(), tmp.size(), MPI_INT, commCpu, FrameworkMacro::transformComm(traId));
               FrameworkMacro::check(ierr, 820);

               // Expand modes to rank map
               std::map<int,int>::iterator mapIt;
               for(int i = 0; i < iData(0); i++)
               {
                  mapIt = mapModes.find(tmp(i));
                  if(mapIt != mapModes.end())
                  {
                     ranks.at(mapIt->second).insert(iData(1));
                  }
               }
            }
         }

         // Synchronize
         FrameworkMacro::synchronize();

         // Initialize storage for SPECTRAL sub communicators
         FrameworkMacro::initSubComm(FrameworkMacro::SPECTRAL, spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>());

         std::set<int>  subRanks;
         int i_ = 0;
         for(size_t i = 0; i < ranks.size(); i++)
         {
            subRanks.clear();
            for(int cpu = 0; cpu < spRes->nCpu(); ++cpu)
            {
               int size;
               if(cpu == FrameworkMacro::id())
               {
                  size = ranks.at(i).size();
                  FrameworkMacro::synchronize();
                  ierr = MPI_Bcast(&size, 1, MPI_INT, cpu, MPI_COMM_WORLD);
                  FrameworkMacro::check(ierr, 821);

                  if(size > 0)
                  {
                     tmp.resize(size);
                     int j = 0;
                     for(std::set<int>::iterator sIt = ranks.at(i).begin(); sIt != ranks.at(i).end(); ++sIt)
                     {
                        tmp(j) = *sIt;
                        ++j;
                        subRanks.insert(*sIt);
                     }
                     FrameworkMacro::synchronize();
                     ierr = MPI_Bcast(tmp.data(), size, MPI_INT, cpu, MPI_COMM_WORLD);
                     FrameworkMacro::check(ierr, 822);
                  }
               } else
               {
                  // Get size
                  FrameworkMacro::synchronize();
                  ierr = MPI_Bcast(&size, 1, MPI_INT, cpu, MPI_COMM_WORLD);
                  FrameworkMacro::check(ierr, 823);

                  // Receive ranks
                  if(size > 0)
                  {
                     tmp.resize(size);
                     FrameworkMacro::synchronize();
                     ierr = MPI_Bcast(tmp.data(), tmp.size(), MPI_INT, cpu, MPI_COMM_WORLD);
                     FrameworkMacro::check(ierr, 824);

                     for(int j = 0; j < size; ++j)
                     {
                        subRanks.insert(tmp(j));
                     }
                  }
               }
            }

            FrameworkMacro::setSubComm(FrameworkMacro::SPECTRAL, i_, subRanks);

            if(i_ < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>() && i == static_cast<size_t>(spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(i_)))
            {
               i_++;
            }
         }

         // Synchronize
         FrameworkMacro::synchronize();
      #endif //defined GEOMHDISCC_MPI && defined GEOMHDISCC_MPISPSOLVE
   }

   void SLFmScheme::addTransformSetups(SharedResolution spRes) const
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

   Transform::SharedFftSetup SLFmScheme::spSetup1D(SharedResolution spRes) const
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

   Transform::SharedPolySetup SLFmScheme::spSetup2D(SharedResolution spRes) const
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

   Transform::SharedFftSetup SLFmScheme::spSetup3D(SharedResolution spRes) const
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

   SLFmScheme::SLFmScheme(const ArrayI& dim)
      : IRegularSHmScheme(dim)
   {
   }

   SLFmScheme::~SLFmScheme()
   {
   }

   void SLFmScheme::setDimensions()
   {
      //
      // Set transform space sizes
      //
      ArrayI traSize(3);
      traSize(0) = this->mI + 1;
      traSize(1) = this->mL + 1;
      traSize(2) = this->mM + 1;
      this->setTransformSpace(traSize);

      //
      // Compute sizes
      //

      // Get standard dealiased FFT size
      int nR = Transform::Fft::ToolsSelector::dealiasCosFft(this->mI+1+8);
      // Check for optimised FFT sizes
      nR = Transform::Fft::ToolsSelector::optimizeFft(nR);

      // Get dealiased associated Legendre transform size
      int nTh = Transform::PolynomialTools::dealias(this->mL+1);

      // Get standard dealiased FFT size
      int nPh = Transform::Fft::ToolsSelector::dealiasMixedFft(this->mM+1);
      // Check for optimised FFT sizes
      nPh = Transform::Fft::ToolsSelector::optimizeFft(nPh);

      //
      // Initialise first transform
      //

      // Initialise forward dimension of first transform
      this->setDimension(nR, Dimensions::Transform::TRA1D, Dimensions::Data::DATF1D);

      // Initialise backward dimension of first transform
      this->setDimension(nR, Dimensions::Transform::TRA1D, Dimensions::Data::DATB1D);

      // Initialise second dimension of first transform
      this->setDimension(traSize(1), Dimensions::Transform::TRA1D, Dimensions::Data::DAT2D);

      // Initialise third dimension of first transform
      this->setDimension(traSize(2), Dimensions::Transform::TRA1D, Dimensions::Data::DAT3D);

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
      this->setDimension(traSize(2), Dimensions::Transform::TRA2D, Dimensions::Data::DAT3D);

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

   void SLFmScheme::setCosts()
   {
      // Set first transform cost
      this->setCost(1.0, Dimensions::Transform::TRA1D);

      // Set second transform cost
      this->setCost(1.0, Dimensions::Transform::TRA2D);

      // Set third transform cost
      this->setCost(1.0, Dimensions::Transform::TRA3D);
   }

   void SLFmScheme::setScalings()
   {
      // Set first transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA1D);

      // Set second transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA2D);

      // Set third transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA3D);
   }

   void SLFmScheme::setMemoryScore()
   {
      // Set first transform memory footprint
      this->setMemory(1.0, Dimensions::Transform::TRA1D);

      // Set second transform memory footprint
      this->setMemory(1.0, Dimensions::Transform::TRA2D);

      // Set third transform memory footprint
      this->setMemory(1.0, Dimensions::Transform::TRA3D);
   }

   bool SLFmScheme::applicable() const
   {
      return true;
   }
}
}
