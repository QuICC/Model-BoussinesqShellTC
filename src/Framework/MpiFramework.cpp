/** 
 * @file MpiFramework.cpp
 * @brief Source of the implementation of an MPI framework
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//
#include "Profiler/ProfilerMacro.h"

// System includes
//
#include <mpi.h>

// External includes
//

// Class include
//
#include "Framework/MpiFramework.hpp"

// Project includes
//
#include "Base/Precision.hpp"
#include "IoHdf5/Hdf5File.hpp"

namespace GeoMHDiSCC {

   void MpiFramework::init()
   {
      // Initialise the profiler if needed
      ProfilerMacro_init();

      // Initialise the precision framework
      Precision::init();

      // Initialise MPI system
      MPI_Init(0, 0);

      // Set MPI rank of local CPU
      MPI_Comm_rank(MPI_COMM_WORLD, &MpiFramework::mCpuId);
   }

   void MpiFramework::setup(const int nCpu)
   {
      // Set the number of CPUs
      MpiFramework::mNCpu = nCpu;

      // Get MPI size
      int size;
      MPI_Comm_size(MPI_COMM_WORLD, &size);

      // Check that the framework setup is right
      MpiFramework::checkFramework(size);
   }

   void MpiFramework::synchronize()
   {
      // Create MPI barrier to force synchronisation
      MPI_Barrier(MPI_COMM_WORLD);
   }

   void MpiFramework::syncSubComm(const MpiFramework::SubCommId id, const int idx)
   {
      assert(MpiFramework::mSubComm.find(id) != MpiFramework::mSubComm.end());

      // Create MPI barrier to force synchronisation
      MPI_Barrier(MpiFramework::mSubComm.find(id)->second.at(idx));
   }

   MPI_Group MpiFramework::getSubGroup(const MpiFramework::SubCommId id, const int idx)
   {
      assert(MpiFramework::mSubGroup.find(id) != MpiFramework::mSubGroup.end());

      return MpiFramework::mSubGroup.find(id)->second.at(idx);
   }

   MPI_Comm MpiFramework::getSubComm(const MpiFramework::SubCommId id, const int idx)
   {
      assert(MpiFramework::mSubComm.find(id) != MpiFramework::mSubComm.end());

      return MpiFramework::mSubComm.find(id)->second.at(idx);
   }

   void MpiFramework::finalize()
   {
      // Make sure all finished and are synchronised
      MPI_Barrier(MPI_COMM_WORLD);

      // Finalize HDF5 MPI data
      IoHdf5::finalizeHdf5();

      // Finalise MPI system
      MPI_Finalize();
   }

   void MpiFramework::setSubComm(const MpiFramework::SubCommId id, const std::vector<std::vector<int> >& ranks)
   {
      assert(ranks.size() > 0);

      MPI_Group   world;
      MPI_Comm_group(MPI_COMM_WORLD, &world);

      MPI_Group group;
      MPI_Comm comm;
      mSubGroup.insert(std::make_pair(SPECTRAL, std::vector<MPI_Group>()));
      mSubComm.insert(std::make_pair(SPECTRAL, std::vector<MPI_Comm>()));

      ArrayI curRanks;
      for(size_t i = 0; i < ranks.size(); ++i)
      {
         assert(ranks.at(i).size() > 0);
         // Create array of ranks
         curRanks.resize(ranks.at(i).size());
         for(int j = 0; j < ranks.at(i).size(); j++)
         {
            curRanks(j) = ranks.at(i).at(j);
         }

         // Create spectral group
         MPI_Group_incl(world, curRanks.size(), curRanks.data(), &group);
         mSubGroup.find(SPECTRAL)->second.push_back(group);

         // Create spectral communicator
         MPI_Comm_create(MPI_COMM_WORLD, group, &comm);
         mSubComm.find(SPECTRAL)->second.push_back(comm);
      }
   }

   std::map<MpiFramework::SubCommId,std::vector<MPI_Group> > MpiFramework::mSubGroup = std::map<MpiFramework::SubCommId,std::vector<MPI_Group> >();

   std::map<MpiFramework::SubCommId,std::vector<MPI_Comm> > MpiFramework::mSubComm = std::map<MpiFramework::SubCommId,std::vector<MPI_Comm> >();

}
