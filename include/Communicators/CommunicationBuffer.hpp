/**
 * @file CommunicationBuffer.hpp
 * @brief Implementation of a "raw" communication buffer 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef COMMUNICATIONBUFFER_HPP
#define COMMUNICATIONBUFFER_HPP

// Configuration includes
//
#include "StorageProfiler/StorageProfilerMacro.h"
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//
#include <cassert>
#include <vector>

// External includes
//

// Project includes
//

namespace GeoMHDiSCC {

namespace Parallel {

   /**
    * 
    * @brief Implementation of a "raw" communication buffer
    */ 
   template <typename TData> class CommunicationBuffer
   {
      public:
         /**
         * @brief Constructor
         */
         CommunicationBuffer();

         /**
         * @brief Destructor
         */
         ~CommunicationBuffer();

         /**
          * @brief Allocate buffers
          *
          * @param sizes      Sizes per node
          * @param maxPacks   Maximum number of packs
          */
         void allocate(const std::vector<int>& sizes, const int maxPacks);

         /**
          * @brief Allocate buffers to fit both requested sizes
          *
          * @param aSizes      Sizes per node
          * @param maxAPacks   Maximum number of packs
          * @param bSizes      Sizes per node
          * @param maxBPacks   Maximum number of packs
          */
         void allocateMax(const std::vector<int>& aSizes, const int maxAPacks, const std::vector<int>& bSizes, const int maxBPacks);

         /**
          * @brief Get pointer to raw buffer storage
          */
         TData* at(const int id);
         
      protected:

      private:
         /**
          * @brief Cleanup the communication buffers
          */
         void cleanupBuffers();

         /**
          * @brief MPI communication buffers
          */
         std::vector<TData *> mData;
   };

   template <typename TData> CommunicationBuffer<TData>::CommunicationBuffer()
   {
   }

   template <typename TData> CommunicationBuffer<TData>::~CommunicationBuffer()
   {
      // Cleanup the communication buffers
      this->cleanupBuffers();
   }

   template <typename TData> void CommunicationBuffer<TData>::allocate(const std::vector<int>& sizes, const int maxPacks)
   {
      // Create CPU group buffers
      for(size_t id = 0; id < sizes.size(); ++id)
      {
         this->mData.push_back(new TData[maxPacks*sizes.at(id)]);
      }

      #ifdef GEOMHDISCC_STORAGEPROFILE
         MHDFloat mem = 0.0;

         for(size_t id = 0; id < sizes.size(); ++id)
         {
            mem += Debug::MemorySize<TData>::BYTES*maxPacks*sizes.at(id);
         }

         StorageProfilerMacro_update(StorageProfilerMacro::MPI, mem);

         #ifdef GEOMHDISCC_STORAGEPROFILER_DETAILED
            StorageProfilerMacro_update(StorageProfilerMacro::MPIBUFFERS, mem);
         #endif // GEOMHDISCC_STORAGEPROFILER_DETAILED
      #endif // GEOMHDISCC_STORAGEPROFILE
   }

   template <typename TData> void CommunicationBuffer<TData>::allocateMax(const std::vector<int>& aSizes, const int maxAPacks, const std::vector<int>& bSizes, const int maxBPacks)
   {
      // Create CPU group buffers
      for(size_t id = 0; id < std::min(aSizes.size(), bSizes.size()); ++id)
      {
         this->mData.push_back(new TData[std::max(maxAPacks*aSizes.at(id), maxBPacks*bSizes.at(id))]);
      }

      // Deal with different number of CPUs in groups
      if(aSizes.size() > bSizes.size())
      {
         for(size_t id = bSizes.size(); id < aSizes.size(); ++id)
         {
            this->mData.push_back(new TData[maxAPacks*aSizes.at(id)]);
         }
      } else if(bSizes.size() > aSizes.size())
      {
         for(size_t id = aSizes.size(); id < bSizes.size(); ++id)
         {
            this->mData.push_back(new TData[maxBPacks*bSizes.at(id)]);
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

   template <typename TData> TData* CommunicationBuffer<TData>::at(const int id)
   {
      // Safety assert
      assert(this->mData.size() > static_cast<size_t>(id));

      return this->mData.at(id);
   }

   template <typename TData> void CommunicationBuffer<TData>::cleanupBuffers()
   {
      // Free the buffers memory
      for(size_t id = 0; id < this->mData.size(); ++id)
      {
         delete[] this->mData.at(id);
      }
   }

}
}

#endif // COMMUNICATIONBUFFER_HPP
