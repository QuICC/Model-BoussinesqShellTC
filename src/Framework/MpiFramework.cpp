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

      // Initialise local sub communicator
      MPI_Group   world;
      MPI_Comm_group(MPI_COMM_WORLD, &world);

      // Create local group
      MPI_Group group;
      MPI_Group_incl(world, 1, &MpiFramework::mCpuId, &group);
      mSubGroup.insert(std::make_pair(LOCAL, group));

      // Create local communicator
      MPI_Comm comm;
      MPI_Comm_create(MPI_COMM_WORLD, MpiFramework::getSubGroup(LOCAL), &comm);
      mSubComm.insert(std::make_pair(LOCAL, comm));
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

   void MpiFramework::syncSubComm(const MpiFramework::SubCommId id)
   {
      assert(MpiFramework::mSubComm.find(id) != MpiFramework::mSubComm.end());

      // Create MPI barrier to force synchronisation
      MPI_Barrier(MpiFramework::mSubComm.find(id)->second);
   }

   MPI_Group MpiFramework::getSubGroup(const MpiFramework::SubCommId id)
   {
      assert(MpiFramework::mSubGroup.find(id) != MpiFramework::mSubGroup.end());

      return MpiFramework::mSubGroup.find(id)->second;
   }

   MPI_Comm MpiFramework::getSubComm(const MpiFramework::SubCommId id)
   {
      assert(MpiFramework::mSubComm.find(id) != MpiFramework::mSubComm.end());

      return MpiFramework::mSubComm.find(id)->second;
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

   void MpiFramework::setSpectralComm(const std::vector<int>& ranks)
   {
      if(ranks.size() > 0)
      {
         MPI_Group   world;
         MPI_Comm_group(MPI_COMM_WORLD, &world);

         // Create spectral group
         MPI_Group group;
         MPI_Group_incl(world, ranks.size(), &ranks.front(), &group);
         mSubGroup.insert(std::make_pair(SPECTRAL, group));

         // Create spectral communicator
         MPI_Comm comm;
         MPI_Comm_create(MPI_COMM_WORLD, MpiFramework::getSubGroup(SPECTRAL), &comm);
         mSubComm.insert(std::make_pair(SPECTRAL, comm));
      
      // Spectral MPI group and communicator are WORLD
      } else
      {
         MPI_Group   world;
         MPI_Comm_group(MPI_COMM_WORLD, &world);

         mSubGroup.insert(std::make_pair(SPECTRAL, world));

         mSubComm.insert(std::make_pair(SPECTRAL, MPI_COMM_WORLD));
      }
   }

   std::map<MpiFramework::SubCommId,MPI_Group> MpiFramework::mSubGroup = std::map<MpiFramework::SubCommId,MPI_Group>();

   std::map<MpiFramework::SubCommId,MPI_Comm> MpiFramework::mSubComm = std::map<MpiFramework::SubCommId,MPI_Comm>();

}
