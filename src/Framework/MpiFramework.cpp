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

   void MpiFramework::addTransformComm(const ArrayI& ids)
   {
      MpiFramework::mTransformCpus.push_back(ids);

      MPI_Group world;
      MPI_Group group;
      MPI_Comm_group(MPI_COMM_WORLD, &world);
      #if defined GEOMHDISCC_MPIIMPL_MVAPICH || defined GEOMHDISCC_MPIIMPL_MPICH
         MPI_Group_incl(world, ids.size(), const_cast<int*>(ids.data()), &group);
      #else
         MPI_Group_incl(world, ids.size(), ids.data(), &group);
      #endif //defined GEOMHDISCC_MPIIMPL_MVAPICH || defined GEOMHDISCC_MPIIMPL_MPICH
      MPI_Comm comm;
      MPI_Comm_create(MPI_COMM_WORLD, group, &comm);

      MpiFramework::mTransformGroups.push_back(group);

      MpiFramework::mTransformComms.push_back(comm);

      int rank;
      MPI_Comm_rank(comm, &rank);
      MpiFramework::mTransformIds.push_back(rank);
   }

   void MpiFramework::syncTransform(const int traId)
   {
      MPI_Barrier(MpiFramework::mTransformComms.at(traId));
   }

   int MpiFramework::transformId(const int traId)
   {
      return MpiFramework::mTransformIds.at(traId);
   }

   const ArrayI& MpiFramework::transformCpus(const int traId)
   {
      return MpiFramework::mTransformCpus.at(traId);
   }

   MPI_Group MpiFramework::transformGroup(const int traId)
   {
      return MpiFramework::mTransformGroups.at(traId);
   }

   MPI_Comm MpiFramework::transformComm(const int traId)
   {
      return MpiFramework::mTransformComms.at(traId);
   }

   void MpiFramework::syncSubComm(const MpiFramework::SubCommId id, const int idx)
   {
      assert(MpiFramework::mSubComm.find(id) != MpiFramework::mSubComm.end());
      assert(MpiFramework::mSubComm.find(id)->second.size() > static_cast<size_t>(idx));

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

   void MpiFramework::initSubComm(const MpiFramework::SubCommId id, const int size)
   {
      // Add group and communicator
      mSubGroup.insert(std::make_pair(id, std::vector<MPI_Group>()));
      mSubComm.insert(std::make_pair(id, std::vector<MPI_Comm>()));

      for(int i = 0; i < size; i++)
      {
         mSubGroup.find(id)->second.push_back(MPI_GROUP_NULL);
         mSubComm.find(id)->second.push_back(MPI_COMM_NULL);
      }
   }

   void MpiFramework::setSubComm(const MpiFramework::SubCommId id, const int idx, const std::set<int>& ranks)
   {
      assert(ranks.size() > 0);
      assert(mSubGroup.find(id) != mSubGroup.end());
      assert(mSubComm.find(id) != mSubComm.end());

      MPI_Group   world;
      MPI_Comm_group(MPI_COMM_WORLD, &world);

      MPI_Group group;
      MPI_Comm comm;

      // Create array of ranks
      ArrayI curRanks(ranks.size());
      int j = 0;
      for(std::set<int>::iterator sIt = ranks.begin(); sIt != ranks.end(); ++sIt)
      {
         curRanks(j) = *sIt;
         j++;
      }

      // Create sub group
      MPI_Group_incl(world, curRanks.size(), curRanks.data(), &group);
      // Create sub communicator
      MPI_Comm_create(MPI_COMM_WORLD, group, &comm);

      if(comm != MPI_COMM_NULL)
      {
         assert(mSubGroup.find(id)->second.size() > static_cast<size_t>(idx));
         assert(mSubComm.find(id)->second.size() > static_cast<size_t>(idx));

         mSubGroup.find(SPECTRAL)->second.at(idx) = group;
         mSubComm.find(SPECTRAL)->second.at(idx) = comm;
      }

      MpiFramework::synchronize();
   }

   void MpiFramework::sleep(const MHDFloat seconds)
   {
      // Sleep MPI process for given amount of seconds
      MHDFloat start = MPI_Wtime();
      MHDFloat tmp = 0.0;
      while(MPI_Wtime() - start < seconds)
      {
         // Do Nothing ... just waiting
         tmp += 1.0;
      }
   }

   std::vector<int>  MpiFramework::mTransformIds = std::vector<int>();

   std::vector<ArrayI>  MpiFramework::mTransformCpus = std::vector<ArrayI>();

   std::vector<MPI_Group>  MpiFramework::mTransformGroups = std::vector<MPI_Group>();

   std::vector<MPI_Comm>  MpiFramework::mTransformComms = std::vector<MPI_Comm>();

   std::map<MpiFramework::SubCommId,std::vector<MPI_Group> > MpiFramework::mSubGroup = std::map<MpiFramework::SubCommId,std::vector<MPI_Group> >();

   std::map<MpiFramework::SubCommId,std::vector<MPI_Comm> > MpiFramework::mSubComm = std::map<MpiFramework::SubCommId,std::vector<MPI_Comm> >();

}
