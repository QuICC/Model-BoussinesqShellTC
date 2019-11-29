/** 
 * @file WFTScheme.cpp
 * @brief Source of the cylindrical Worland(poly) + Fourier + Chebyshev(FFT) scheme implementation
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include <set>

// External includes
//

// Class include
//
#include "SpatialSchemes/3D/WFTScheme.hpp"

// Project includes
//
#include "PolynomialTransforms/PolynomialTools.hpp"

namespace QuICC {

namespace Schemes {
   
   std::string WFTScheme::type()
   {
      return "WFT";
   }

   void WFTScheme::tuneResolution(SharedResolution spRes, const Parallel::SplittingDescription& descr)
   {
      WFTScheme::tuneMpiResolution(descr);
      
      // Create spectral space sub communicators
      #if defined QUICC_MPI && defined QUICC_MPISPSOLVE
         // MPI error code
         int ierr;

         // Get world group
         MPI_Group world;
         MPI_Group group;
         ierr = MPI_Comm_group(MPI_COMM_WORLD, &world);
         FrameworkMacro::check(ierr, 811);

         // Create minimial MPI group
         ierr = MPI_Group_incl(world, FrameworkMacro::transformCpus(0).size(), FrameworkMacro::transformCpus(0).data(), &group);
         FrameworkMacro::check(ierr, 812);

         // Create minimial MPI communicator
         MPI_Comm comm;
         ierr = MPI_Comm_create(MPI_COMM_WORLD, group, &comm);
         FrameworkMacro::check(ierr, 813);

         // Initialise the ranks with local rank
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

         // Loop over all cpus
         int commId;
         int globalCpu = FrameworkMacro::id();
         ierr = MPI_Comm_rank(comm, &commId); 
         FrameworkMacro::check(ierr, 814);
         ArrayI tmp;
         for(int commCpu = 0; commCpu < FrameworkMacro::transformCpus(0).size(); ++commCpu)
         {
            int size;
            if(commCpu == commId)
            {
               // Send the size
               size = modes.size();
               ierr = MPI_Bcast(&size, 1, MPI_INT, commCpu, comm);
               FrameworkMacro::check(ierr, 815);
               MPI_Barrier(comm);

               // Send global CPU rank 
               globalCpu = FrameworkMacro::id();
               ierr = MPI_Bcast(&globalCpu, 1, MPI_INT, commCpu, comm);
               FrameworkMacro::check(ierr, 816);
               MPI_Barrier(comm);

               // Send modes
               ierr = MPI_Bcast(modes.data(), modes.size(), MPI_INT, commCpu, comm);
               FrameworkMacro::check(ierr, 817);
               MPI_Barrier(comm);
            } else
            {
               // Get size
               ierr = MPI_Bcast(&size, 1, MPI_INT, commCpu, comm);
               FrameworkMacro::check(ierr, 818);
               MPI_Barrier(comm);

               // Get global CPU rank 
               ierr = MPI_Bcast(&globalCpu, 1, MPI_INT, commCpu, comm);
               FrameworkMacro::check(ierr, 819);
               MPI_Barrier(comm);

               // Receive modes
               tmp.resize(size);
               ierr = MPI_Bcast(tmp.data(), tmp.size(), MPI_INT, commCpu, comm);
               FrameworkMacro::check(ierr, 820);
               MPI_Barrier(comm);

               std::map<int,int>::iterator mapIt;
               for(int i = 0; i < size; i++)
               {
                  mapIt = mapModes.find(tmp(i));
                  if(mapIt != mapModes.end())
                  {
                     ranks.at(mapIt->second).insert(globalCpu);
                  }
               }
            }

            // Synchronize
            FrameworkMacro::synchronize();
         }

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
                  ierr = MPI_Bcast(&size, 1, MPI_INT, cpu, MPI_COMM_WORLD);
                  FrameworkMacro::check(ierr, 821);
                  FrameworkMacro::synchronize();

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
                     ierr = MPI_Bcast(tmp.data(), size, MPI_INT, cpu, MPI_COMM_WORLD);
                     FrameworkMacro::check(ierr, 822);
                     FrameworkMacro::synchronize();
                  }
               } else
               {
                  // Get size
                  ierr = MPI_Bcast(&size, 1, MPI_INT, cpu, MPI_COMM_WORLD);
                  FrameworkMacro::check(ierr, 823);
                  FrameworkMacro::synchronize();

                  // Receive ranks
                  if(size > 0)
                  {
                     tmp.resize(size);
                     ierr = MPI_Bcast(tmp.data(), tmp.size(), MPI_INT, cpu, MPI_COMM_WORLD);
                     FrameworkMacro::check(ierr, 824);
                     FrameworkMacro::synchronize();

                     for(int j = 0; j < size; ++j)
                     {
                        subRanks.insert(tmp(j));
                     }
                  }
               }

               // Synchronize
               FrameworkMacro::synchronize();
            }

            FrameworkMacro::setSubComm(FrameworkMacro::SPECTRAL, i_, subRanks);

            if(i_ < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>() && i == static_cast<size_t>(spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(i_)))
            {
               i_++;
            }
         }
         
         // Free communicator
         ierr = MPI_Comm_free(&comm);
         FrameworkMacro::check(ierr, 825);
      #endif //defined QUICC_MPI && defined QUICC_MPISPSOLVE
   }

   void WFTScheme::addTransformSetups(SharedResolution spRes) const
   {
      // Add setup for first transform
      Transform::SharedPolySetup  spS1D = this->spSetup1D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA1D, spS1D);

      // Add setup for second transform
      Transform::SharedFftSetup  spS2D = this->spSetup2D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA2D, spS2D);

      // Add setup for third transform
      Transform::SharedFftSetup  spS3D = this->spSetup3D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA3D, spS3D);
   }

   Transform::SharedPolySetup WFTScheme::spSetup1D(SharedResolution spRes) const
   {
      // Get physical size of polynomial transform
      int size = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DATF1D>();

      // Get spectral size of the polynomial transform
      int specSize = spRes->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

      // Storage for the list of indexes
      std::vector<ArrayI>  fast;
      fast.reserve(spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>());
      ArrayI  slow(spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>());

      // Multiplier from second dimension 
      ArrayI mult(spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>());

      // Get number of transforms and list of indexes
      int howmany = 0;
      for(int i = 0; i < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); i++)
      {
         howmany += spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(i);

         slow(i) = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(i);

         fast.push_back(ArrayI(spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DATB1D>(i)));
         for(int j = 0; j < fast.at(i).size(); j++)
         {
            fast.at(i)(j) = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DATB1D>(j,i);
         }

         mult(i) = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(i);
      }

      int padSize = this->dim(Dimensions::Transform::TRA1D, Dimensions::Data::DATB1D) - specSize;

      return Transform::SharedPolySetup(new Transform::PolySetup(size, howmany, specSize, fast, slow, mult, padSize));
   }

   Transform::SharedFftSetup WFTScheme::spSetup2D(SharedResolution spRes) const
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

   Transform::SharedFftSetup WFTScheme::spSetup3D(SharedResolution spRes) const
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

      return Transform::SharedFftSetup(new Transform::FftSetup(size, howmany, specSize, Transform::FftSetup::REAL));
   }

   WFTScheme::WFTScheme(const ArrayI& dim)
      : IRegular3DScheme(dim)
   {
   }

   WFTScheme::~WFTScheme()
   {
   }

   void WFTScheme::setDimensions()
   {
      //
      // Set transform space sizes
      //
      ArrayI traSize(3);
      traSize(0) = this->mI + 1;
      traSize(1) = this->mJ + 1;
      traSize(2) = this->mK + 1;
      this->setTransformSpace(traSize);

      //
      // Compute sizes
      //

      // Get dealiased Worland size
      int nR = Transform::PolynomialTools::dealias(this->mI+this->mJ/2+8);

      // Get mixed dealiased FFT size
      int nTh = Transform::Fft::ToolsSelector::dealiasMixedFft(this->mJ+1);
      // Check for optimised FFT sizes
      nTh = Transform::Fft::ToolsSelector::optimizeFft(nTh);

      // Get standard dealiased FFT size
      int nZ = Transform::Fft::ToolsSelector::dealiasCosFft(this->mK+1);
      // Check for optimised FFT sizes
      nZ = Transform::Fft::ToolsSelector::optimizeFft(nZ);

      //
      // Initialise first transform
      //

      // Initialise forward dimension of first transform
      this->setDimension(nR, Dimensions::Transform::TRA1D, Dimensions::Data::DATF1D);

      // Initialise backward dimension of first transform
      this->setDimension(this->mI+1, Dimensions::Transform::TRA1D, Dimensions::Data::DATB1D);

      // Initialise second dimension of first transform
      this->setDimension(traSize(2), Dimensions::Transform::TRA1D, Dimensions::Data::DAT2D);

      // Initialise third dimension of first transform
      this->setDimension(traSize(1), Dimensions::Transform::TRA1D, Dimensions::Data::DAT3D);

      //
      // Initialise second transform
      //

      // Initialise forward dimension of second transform
      this->setDimension(nTh, Dimensions::Transform::TRA2D, Dimensions::Data::DATF1D);

      // Initialise backward dimension of second transform
      this->setDimension(nTh/2 + 1, Dimensions::Transform::TRA2D, Dimensions::Data::DATB1D);

      // Initialise second dimension of second transform
      this->setDimension(nR, Dimensions::Transform::TRA2D, Dimensions::Data::DAT2D);

      // Initialise third dimension of second transform
      this->setDimension(traSize(2), Dimensions::Transform::TRA2D, Dimensions::Data::DAT3D);

      //
      // Initialise third transform
      //

      // Initialise forward dimension of third transform
      this->setDimension(nZ, Dimensions::Transform::TRA3D, Dimensions::Data::DATF1D);

      // Initialise backward dimension of third transform
      this->setDimension(nZ, Dimensions::Transform::TRA3D, Dimensions::Data::DATB1D);

      // Initialise second dimension of third transform
      this->setDimension(nTh, Dimensions::Transform::TRA3D, Dimensions::Data::DAT2D);

      // Initialise third dimension of third transform
      this->setDimension(nR, Dimensions::Transform::TRA3D, Dimensions::Data::DAT3D);
   }

   void WFTScheme::setCosts()
   {
      // Set first transform cost
      this->setCost(1.0, Dimensions::Transform::TRA1D);

      // Set second transform cost
      this->setCost(1.0, Dimensions::Transform::TRA2D);

      // Set third transform cost
      this->setCost(1.0, Dimensions::Transform::TRA3D);
   }

   void WFTScheme::setScalings()
   {
      // Set first transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA1D);

      // Set second transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA2D);

      // Set third transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA3D);
   }

   void WFTScheme::setMemoryScore()
   {
      // Set first transform memory footprint
      this->setMemory(1.0, Dimensions::Transform::TRA1D);

      // Set second transform memory footprint
      this->setMemory(1.0, Dimensions::Transform::TRA2D);

      // Set third transform memory footprint
      this->setMemory(1.0, Dimensions::Transform::TRA3D);
   }

   bool WFTScheme::applicable() const
   {
      return true;
   }
}
}
