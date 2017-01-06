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
         TData* data();

         /**
          * @brief Get pointer to raw sub buffer start
          */
         TData* at(const int id);

         /**
          * @brief Get zero positions for sub buffer
          */
         int* zero();

         /**
          * @brief Get current position in sub buffer
          */
         int& pos(const int id);

         /**
          * @brief Reset positions in sub buffers
          */
         void resetPositions();

         /**
          * @brief Get total available storage
          */
         int total() const;
         
      protected:

      private:
         /**
          * @brief Cleanup the communication buffers
          */
         void cleanupBuffers();

         /**
          * @brief Total memory in buffer
          */
         int mTotal;

         /**
          * @brief MPI communication buffers
          */
         TData * mData;

         /**
          * @brief Start position for sub buffers
          */
         std::vector<int>  mZero;

         /**
          * @brief Current position in sub buffers
          */
         std::vector<int>  mPos;
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
      this->mTotal = 0;
      for(std::vector<int>::const_iterator it = sizes.begin(); it != sizes.end(); ++it)
      {
         // Create zero position and initialize position
         this->mZero.push_back(this->mTotal);
         this->mPos.push_back(0);

         this->mTotal += (*it)*maxPacks;
      }

      // Allocate large buffer
      this->mData = new TData[this->mTotal];

      #ifdef QUICC_STORAGEPROFILE
         MHDFloat mem = 0.0;

         mem += Debug::MemorySize<TData>::BYTES*this->mTotal;

         StorageProfilerMacro_update(StorageProfilerMacro::MPI, mem);

         #ifdef QUICC_STORAGEPROFILER_DETAILED
            StorageProfilerMacro_update(StorageProfilerMacro::MPIBUFFERS, mem);
         #endif // QUICC_STORAGEPROFILER_DETAILED
      #endif // QUICC_STORAGEPROFILE
   }

   template <typename TData> void CommunicationBuffer<TData>::allocateMax(const std::vector<int>& aSizes, const int maxAPacks, const std::vector<int>& bSizes, const int maxBPacks)
   {
      // Create CPU group buffers
      this->mTotal = 0;
      for(size_t id = 0; id < std::min(aSizes.size(), bSizes.size()); ++id)
      {
         // Create zero position and initialize position
         this->mZero.push_back(this->mTotal);
         this->mPos.push_back(0);

         this->mTotal += std::max(maxAPacks*aSizes.at(id), maxBPacks*bSizes.at(id));
      }

      // Deal with different number of CPUs in groups
      if(aSizes.size() > bSizes.size())
      {
         for(size_t id = bSizes.size(); id < aSizes.size(); ++id)
         {
            // Create zero position and initialize position
            this->mZero.push_back(this->mTotal);
            this->mPos.push_back(0);

            this->mTotal += maxAPacks*aSizes.at(id);
         }
      } else if(bSizes.size() > aSizes.size())
      {
         for(size_t id = aSizes.size(); id < bSizes.size(); ++id)
         {
            // Create zero position and initialize position
            this->mZero.push_back(this->mTotal);
            this->mPos.push_back(0);

            this->mTotal += maxBPacks*bSizes.at(id);
         }
      }

      this->mData = new TData[this->mTotal];

      #ifdef QUICC_STORAGEPROFILE
         MHDFloat mem = 0.0;

         mem += Debug::MemorySize<TData>::BYTES*this->mTotal;

         StorageProfilerMacro_update(StorageProfilerMacro::MPI, mem);

         #ifdef QUICC_STORAGEPROFILER_DETAILED
            StorageProfilerMacro_update(StorageProfilerMacro::MPIBUFFERS, mem);
         #endif // QUICC_STORAGEPROFILER_DETAILED
      #endif // QUICC_STORAGEPROFILE
   }

   template <typename TData> inline int CommunicationBuffer<TData>::total() const
   {
      return this->mTotal;
   }

   template <typename TData> inline TData* CommunicationBuffer<TData>::data()
   {
      return this->mData;
   }

   template <typename TData> inline int* CommunicationBuffer<TData>::zero()
   {
      return &this->mZero[0];
   }

   template <typename TData> inline TData* CommunicationBuffer<TData>::at(const int id)
   {
      return (this->mData + this->mZero.at(id));
   }

   template <typename TData> void CommunicationBuffer<TData>::resetPositions()
   {
      for(std::vector<int>::iterator it = this->mPos.begin(); it != this->mPos.end(); ++it)
      {
         *it = 0;
      }
   }

   template <typename TData> inline int& CommunicationBuffer<TData>::pos(const int id)
   {
      // Safety assert
      assert(this->mPos.size() > static_cast<size_t>(id));

      return this->mPos.at(id);
   }

   template <typename TData> void CommunicationBuffer<TData>::cleanupBuffers()
   {
      // Free the buffers memory
      delete[] this->mData;
   }

}
}

#endif // COMMUNICATIONBUFFER_HPP
