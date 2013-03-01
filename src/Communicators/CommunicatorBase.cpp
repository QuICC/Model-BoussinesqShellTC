/** \file CommunicatorBase.cpp
 *  \brief Source of the implementation of the communicator base
 */

// Configuration includes
//
#include "StorageProfiler/StorageProfilerMacro.h"

// System includes
//

// External includes
//

// Class include
//
#include "Communicators/CommunicatorBase.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Parallel {

   CommunicatorBase::CommunicatorBase()
   {
   }

   CommunicatorBase::~CommunicatorBase()
   {
      // Cleanup the communication buffers
      this->cleanupBuffers();
   }

   void CommunicatorBase::createBuffers(const int nBuffers)
   {
      for(int i = 0; i < nBuffers; i++)
      {
         // Create buffer
         this->mBuffers.push_back(std::vector<char *>());
      }
   }

   void CommunicatorBase::allocateBuffers(const int idx, const std::vector<int>& sizes, const int maxPacks)
   {
      // Create CPU group buffers
      for(unsigned int id = 0; id < sizes.size(); ++id)
      {
         this->mBuffers.at(idx).push_back(new char[maxPacks*sizes.at(id)]);
      }

      #ifdef GEOMHDISCC_STORAGEPROFILE
         MHDFloat mem = 0.0;

         for(unsigned int id = 0; id < sizes.size(); ++id)
         {
            mem += maxPacks*sizes.at(id);
         }

         StorageProfilerMacro_update(StorageProfilerMacro::MPI, mem);

         #ifdef GEOMHDISCC_STORAGEPROFILER_DETAILED
            StorageProfilerMacro_update(StorageProfilerMacro::MPIBUFFERS, mem);
         #endif // GEOMHDISCC_STORAGEPROFILER_DETAILED
      #endif // GEOMHDISCC_STORAGEPROFILE
   }

   void CommunicatorBase::allocateBuffers(const int idx, const std::vector<int>& aSizes, const int maxAPacks, const std::vector<int>& bSizes, const int maxBPacks)
   {
      // Create CPU group buffers
      for(unsigned int id = 0; id < std::min(aSizes.size(), bSizes.size()); ++id)
      {
         this->mBuffers.at(idx).push_back(new char[std::max(maxAPacks*aSizes.at(id), maxBPacks*bSizes.at(id))]);
      }

      // Deal with different number of CPUs in groups
      if(aSizes.size() > bSizes.size())
      {
         for(unsigned int id = bSizes.size(); id < aSizes.size(); ++id)
         {
            this->mBuffers.at(idx).push_back(new char[maxAPacks*aSizes.at(id)]);
         }
      } else if(bSizes.size() > aSizes.size())
      {
         for(unsigned int id = aSizes.size(); id < bSizes.size(); ++id)
         {
            this->mBuffers.at(idx).push_back(new char[maxBPacks*bSizes.at(id)]);
         }
      }

      #ifdef GEOMHDISCC_STORAGEPROFILE
         MHDFloat mem = 0.0;

         for(unsigned int id = 0; id < std::min(aSizes.size(), bSizes.size()); ++id)
         {
            mem += std::max(maxAPacks*aSizes.at(id), maxBPacks*bSizes.at(id));
         }

         // Deal with different number of CPUs in groups
         if(aSizes.size() > bSizes.size())
         {
            for(unsigned int id = bSizes.size(); id < aSizes.size(); ++id)
            {
               mem += maxAPacks*aSizes.at(id);
            }
         } else if(bSizes.size() > aSizes.size())
         {
            for(unsigned int id = aSizes.size(); id < bSizes.size(); ++id)
            {
               mem += maxBPacks*bSizes.at(id);
            }
         }

         StorageProfilerMacro_update(StorageProfilerMacro::MPI, mem);

         #ifdef GEOMHDISCC_STORAGEPROFILER_DETAILED
            StorageProfilerMacro_update(StorageProfilerMacro::MPIBUFFERS, mem);
         #endif // GEOMHDISCC_STORAGEPROFILER_DETAILED
      #endif // GEOMHDISCC_STORAGEPROFILE
   }

   void CommunicatorBase::cleanupBuffers()
   {
      // loop over all buffers
      for(unsigned int n = 0; n < this->mBuffers.size(); ++n)
      {
         // Free the buffers memory
         for(unsigned int id = 0; id < this->mBuffers.at(n).size(); ++id)
         {
            delete[] this->mBuffers.at(n).at(id);
         }
      }
   }

}
}
