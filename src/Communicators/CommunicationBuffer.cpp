/** 
 * @file CommunicationBuffer.cpp
 * @brief Source of the implementation of the "raw" communication buffer
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//
#include "StorageProfiler/StorageProfilerMacro.h"

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "Communicators/CommunicationBuffer.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Parallel {

   CommunicationBuffer::CommunicationBuffer()
   {
   }

   CommunicationBuffer::~CommunicationBuffer()
   {
      // Cleanup the communication buffers
      this->cleanupBuffers();
   }

   void CommunicationBuffer::allocate(const std::vector<int>& sizes, const int maxPacks)
   {
      // Create CPU group buffers
      for(size_t id = 0; id < sizes.size(); ++id)
      {
         this->mData.push_back(new char[maxPacks*sizes.at(id)]);
      }

      #ifdef GEOMHDISCC_STORAGEPROFILE
         MHDFloat mem = 0.0;

         for(size_t id = 0; id < sizes.size(); ++id)
         {
            mem += maxPacks*sizes.at(id);
         }

         StorageProfilerMacro_update(StorageProfilerMacro::MPI, mem);

         #ifdef GEOMHDISCC_STORAGEPROFILER_DETAILED
            StorageProfilerMacro_update(StorageProfilerMacro::MPIBUFFERS, mem);
         #endif // GEOMHDISCC_STORAGEPROFILER_DETAILED
      #endif // GEOMHDISCC_STORAGEPROFILE
   }

   void CommunicationBuffer::allocateMax(const std::vector<int>& aSizes, const int maxAPacks, const std::vector<int>& bSizes, const int maxBPacks)
   {
      // Create CPU group buffers
      for(size_t id = 0; id < std::min(aSizes.size(), bSizes.size()); ++id)
      {
         this->mData.push_back(new char[std::max(maxAPacks*aSizes.at(id), maxBPacks*bSizes.at(id))]);
      }

      // Deal with different number of CPUs in groups
      if(aSizes.size() > bSizes.size())
      {
         for(size_t id = bSizes.size(); id < aSizes.size(); ++id)
         {
            this->mData.push_back(new char[maxAPacks*aSizes.at(id)]);
         }
      } else if(bSizes.size() > aSizes.size())
      {
         for(size_t id = aSizes.size(); id < bSizes.size(); ++id)
         {
            this->mData.push_back(new char[maxBPacks*bSizes.at(id)]);
         }
      }

      #ifdef GEOMHDISCC_STORAGEPROFILE
         MHDFloat mem = 0.0;

         for(size_t id = 0; id < std::min(aSizes.size(), bSizes.size()); ++id)
         {
            mem += std::max(maxAPacks*aSizes.at(id), maxBPacks*bSizes.at(id));
         }

         // Deal with different number of CPUs in groups
         if(aSizes.size() > bSizes.size())
         {
            for(size_t id = bSizes.size(); id < aSizes.size(); ++id)
            {
               mem += maxAPacks*aSizes.at(id);
            }
         } else if(bSizes.size() > aSizes.size())
         {
            for(size_t id = aSizes.size(); id < bSizes.size(); ++id)
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

   char* CommunicationBuffer::at(const int id)
   {
      // Safety assert
      assert(this->mData.size() > static_cast<size_t>(id));

      return this->mData.at(id);
   }

   void CommunicationBuffer::cleanupBuffers()
   {
      // Free the buffers memory
      for(size_t id = 0; id < this->mData.size(); ++id)
      {
         delete[] this->mData.at(id);
      }
   }

}
}
