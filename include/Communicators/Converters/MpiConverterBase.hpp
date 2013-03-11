/** \file MpiConverterBase.hpp
 *  \brief Templated implementation of the base of a MPI data converter.
 */

#ifndef MPICONVERTERBASE_HPP
#define MPICONVERTERBASE_HPP

// Debug includes
//
#include "StorageProfiler/StorageProfilerMacro.h"
#include "Profiler/ProfilerMacro.h"

// Configuration includes
//
#include "Framework/FrameworkMacro.h"

// System includes
//
#include <cassert>
#include <set>
#include <map>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MpiTypes.hpp"
#include "Resolutions/Resolution.hpp"
#include "Communicators/Converters/IConverter.hpp"
#include "Communicators/CommunicationBuffer.hpp"

namespace GeoMHDiSCC {

namespace Parallel {

   /**
    * \brief Templated implementation of the base of a MPI data converter.
    */
   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> class MpiConverterBase: public IConverter<TFwdA,TBwdA,TFwdB,TBwdB>
   {
      public:
         /**
          * @brief Constructor
          */
         MpiConverterBase();

         /**
          * @brief Destructor
          */
         virtual ~MpiConverterBase();

         /**
          * @brief Set the communication buffers
          *
          * @brief spFwd Forward communication buffers
          * @brief spBwd Backward communication buffers
          */
         void setBuffers(SharedCommunicationBuffer spFwd, SharedCommunicationBuffer spBwd);

         /**
          * @brief Get forward buffer sizes
          */
         const std::vector<int> & fwdSizes() const;

         /**
          * @brief Get backward buffer sizes
          */
         const std::vector<int> & bwdSizes() const;
         
      protected:
         /**
          * @brief Reset Receive positions
          */
         void resetRecvPositions();

         /**
          * @brief Reset Send positions
          */
         void resetSendPositions();

         /**
          * @brief Get a pointer to the receive backward requests
          *
          * @param size Pack size of the requested request
          */
         MPI_Request * pRecvBRequests(const int size);

         /**
          * @brief Get a pointer to the receive forward requests
          *
          * @param size Pack size of the requested request
          */
         MPI_Request * pRecvFRequests(const int size);

         /**
          * @brief Get a pointer to the send backward requests
          *
          * @param size Pack size of the requested request
          */
         MPI_Request * pSendBRequests(const int size);

         /**
          * @brief Get a pointer to the send forward requests
          *
          * @param size Pack size of the requested request
          */
         MPI_Request * pSendFRequests(const int size);

         /**
          * @brief Get ring recv source for id
          *
          * @param id   ID of the CPU
          * @param ref  ID of the reference CPU
          * @param size Size of the CPU group
          */
         int  recvSrc(const int id, const int ref, const int size) const;

         /**
          * @brief Get ring send destination for id
          *
          * @param id   ID of the CPU
          * @param ref  ID of the reference CPU
          * @param size Size of the CPU group
          */
         int  sendDest(const int id, const int ref, const int size) const;

         /**
          * @brief Initialise the positions
          */
         void initPositions();

         /**
          * @brief Setup the MPI communication requests requests
          */
         void setupRequests();

         /**
          * @brief Cleanup the MPI communication requests
          */
         void cleanupRequests();

         /**
          * @brief The number of packs in the "previous/active" forward send
          */
         int   mActiveFSendPacks;

         /**
          * @brief The number of packs in the "previous/active" backward send
          */
         int   mActiveBSendPacks;

         /**
          * @brief Storage for the receive position pointers
          */
         std::vector<int>  mRecvPositions;

         /**
          * @brief Storage for the send position pointers
          */
         std::vector<int>  mSendPositions;

         /**
          * @brief Storage for the non blocking communication requests: Recv F
          */
         std::map<int, std::vector<MPI_Request> >  mRecvFRequests;

         /**
          * @brief Storage for the non blocking communication requests: Recv B
          */
         std::map<int, std::vector<MPI_Request> >  mRecvBRequests;

         /**
          * @brief Storage for the non blocking communication requests: Send F
          */
         std::map<int, std::vector<MPI_Request> >  mSendFRequests;

         /**
          * @brief Storage for the non blocking communication requests: Send B
          */
         std::map<int, std::vector<MPI_Request> >  mSendBRequests;



         /**
          * @brief Size of the forward packet
          *
          * @param id ID of the node
          */
         int sizeFPacket(const int id) const;

         /**
          * @brief Size of the backward packet
          *
          * @param id ID of the node
          */
         int sizeBPacket(const int id) const;

         /**
          * @brief Get size of the forward CPU group
          */
         int nFCpu() const;

         /**
          * @brief Get size of the backward CPU group
          */
         int nBCpu() const;

         /**
          * @brief Get global MPI rank of CPU from forward CPU group
          *
          * @param id CPU group id
          */
         int fCpu(const int id) const;

         /**
          * @brief Get global MPI rank of CPU from backward CPU group
          *
          * @param id CPU group id
          */
         int bCpu(const int id) const;

         /**
          * @brief Sending communication status
          */
         bool  mIsSending;

         /**
          * @brief Receiving communication status
          */
         bool  mIsReceiving;

         /**
          * @brief List of CPU ranks involved in the forward conversion
          */
         std::vector<int>  mFCpuGroup;

         /**
          * @brief List of CPU ranks involved in the backward conversion
          */
         std::vector<int>  mBCpuGroup;

         /**
          * @brief Forward communication
          */
         SharedCommunicationBuffer mspFBuffers;

         /**
          * @brief Backward communication buffer pointer
          */
         SharedCommunicationBuffer mspBBuffers;

         /**
          * @brief List of the forward buffer sizes
          */
         std::vector<int>  mFSizes;

         /**
          * @brief List of the backward buffer sizes
          */
         std::vector<int>  mBSizes;

         /**
          * @brief Communication packs counter
          */
         int mPacks;

         /**
          * @brief Possible forward transform packs
          */
         ArrayI   mForwardPacks;

         /**
          * @brief Possible backward transform packs
          */
         ArrayI   mBackwardPacks;

      private:
   };

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> inline int MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::nFCpu() const
   {
      return this->mFCpuGroup.size();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> inline int MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::nBCpu() const
   {
      return this->mBCpuGroup.size();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> inline int MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::fCpu(const int id) const
   {
      return this->mFCpuGroup.at(id);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> inline int MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::bCpu(const int id) const
   {
      return this->mBCpuGroup.at(id);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> inline int MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::sizeFPacket(const int id) const
   {
      return this->mPacks*this->mFSizes.at(id);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> inline int MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::sizeBPacket(const int id) const
   {
      return this->mPacks*this->mBSizes.at(id);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> inline MPI_Request * MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::pRecvBRequests(const int size)
   {
      return &(this->mRecvBRequests[size].front());
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> inline MPI_Request * MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::pRecvFRequests(const int size)
   {
      return &(this->mRecvFRequests[size].front());
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> inline MPI_Request * MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::pSendBRequests(const int size)
   {
      return &(this->mSendBRequests[size].front());
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> inline MPI_Request * MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::pSendFRequests(const int size)
   {
      return &(this->mSendFRequests[size].front());
   }
      
   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::MpiConverterBase()
      : IConverter<TFwdA, TBwdA, TFwdB, TBwdB>(), mIsSending(false), mIsReceiving(false), mPacks(0), mActiveFSendPacks(0), mActiveBSendPacks(0)
   {
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::~MpiConverterBase()
   {
      // Cleanup the requests memory
      this->cleanupRequests();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> void MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::initPositions()
   {
      // Get maximum position size
      int maxSize = std::max(this->nFCpu(), this->nBCpu());

      // Initialise the position values
      for(int i = 0; i < maxSize; ++i)
      {
         this->mRecvPositions.push_back(0);

         this->mSendPositions.push_back(0);
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> void MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::resetRecvPositions()
   {
      // Create position iterator
      std::vector<int>::iterator it;

      // Reset all positions to zero
      for(it = this->mRecvPositions.begin(); it != this->mRecvPositions.end(); ++it)
      {
         (*it) = 0;
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> void MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::resetSendPositions()
   {
      // Create position iterator
      std::vector<int>::iterator it;

      // Reset all positions to zero
      for(it = this->mSendPositions.begin(); it != this->mSendPositions.end(); ++it)
      {
         (*it) = 0;
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> int MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::sendDest(const int id, const int ref, const int size) const
   {
      // Create send ring
      return ((id + 1 + ref) % size);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> int MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::recvSrc(const int id, const int ref, const int size) const
   {
      // Create recv ring
      return ((size - 1 - id + ref) % size);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> void MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::setupRequests()
   {
      // Storage for global location flags
      int dest;
      int src;
      int tag;
      // Storage for CPU group location flags
      int grpMe;
      int grpDest;
      int grpSrc;

      // Storage for the number of packs
      int packs;

      // Initialise forward transform requests
      for(int k = 0; k < this->mForwardPacks.size(); ++k)
      {
         // Get the pack size
         packs = this->mForwardPacks(k);

         // Initialise receive forward with empty requests
         this->mRecvFRequests.insert(std::make_pair<int, std::vector<MPI_Request> >(packs, std::vector<MPI_Request>()));
         for(int id = 0; id < this->nFCpu(); ++id)
         {
            this->mRecvFRequests[packs].push_back(MPI_REQUEST_NULL);
         }

         // Initialise send backward with empty requests
         this->mSendBRequests.insert(std::make_pair<int, std::vector<MPI_Request> >(packs, std::vector<MPI_Request>()));
         for(int id = 0; id < this->nBCpu(); ++id)
         {
            this->mSendBRequests[packs].push_back(MPI_REQUEST_NULL);
         }

         // Create receive forward requests
         for(int id = 0; id < this->nFCpu(); ++id)
         {
            // Get CPU group index of local node
            grpMe = (*std::find(this->mFCpuGroup.begin(), this->mFCpuGroup.end(), FrameworkMacro::id()));

            // Get source index in CPU group
            grpSrc = this->recvSrc(id, grpMe, this->nFCpu());

            // Get global MPI source rank
            src = this->fCpu(grpSrc);

            // Set MPI tag
            tag = src;

            //Safety asserts
            assert(grpSrc < this->mpFBuffers->size());
            assert(grpSrc < this->mFSizes.size());
            assert(grpSrc < this->mRecvFRequests[packs].size());

            // initialise the Recv request
            MPI_Recv_init(this->mpFBuffers->at(grpSrc), packs*this->mFSizes.at(grpSrc), MPI_PACKED, src, tag, MPI_COMM_WORLD, &(this->mRecvFRequests[packs].at(grpSrc)));
         }

         // Create send backward requests
         for(int id = 0; id < this->nBCpu(); ++id)
         {
            // Set MPI tag
            tag = FrameworkMacro::id();
            // Get CPU group index of local node
            grpMe = (*std::find(this->mBCpuGroup.begin(), this->mBCpuGroup.end(), tag));
            // Get destination index in CPU group
            grpDest = this->sendDest(id, grpMe, this->nBCpu());
            // Get global MPI destination rank
            dest = this->bCpu(grpDest);

            //Safety asserts
            assert(grpDest < this->mpBBuffers->size());
            assert(grpDest < this->mBSizes.size());
            assert(grpDest < this->mSendBRequests[packs].size());

            // initialise the Send request
            MPI_Send_init(this->mpBBuffers->at(grpDest), packs*this->mBSizes.at(grpDest), MPI_PACKED, dest, tag, MPI_COMM_WORLD, &(this->mSendBRequests[packs].at(grpDest)));
         }
      }

      // Initialise backward transform requests
      for(int k = 0; k < this->mBackwardPacks.size(); ++k)
      {
         // Get the pack size
         packs = this->mBackwardPacks(k);

         // Initialise receive backward with empty requests
         this->mRecvBRequests.insert(std::make_pair<int, std::vector<MPI_Request> >(packs, std::vector<MPI_Request>()));
         for(int id = 0; id < this->nBCpu(); ++id)
         {
            this->mRecvBRequests[packs].push_back(MPI_REQUEST_NULL);
         }

         // Initialise send forward with empty requests
         this->mSendFRequests.insert(std::make_pair<int, std::vector<MPI_Request> >(packs, std::vector<MPI_Request>()));
         for(int id = 0; id < this->nFCpu(); ++id)
         {
            this->mSendFRequests[packs].push_back(MPI_REQUEST_NULL);
         }

         // Create receive backward requests
         for(int id = 0; id < this->nBCpu(); ++id)
         {
            // Get CPU group index of local node
            grpMe = (*std::find(this->mBCpuGroup.begin(), this->mBCpuGroup.end(), FrameworkMacro::id()));
            // Get source index in CPU group
            grpSrc = this->recvSrc(id, grpMe, this->nBCpu());
            // Get global MPI source rank
            src = this->bCpu(grpSrc);
            // Set MPI tag
            tag = src;
            // initialise the Recv request
            MPI_Recv_init(this->mpBBuffers->at(grpSrc), packs*this->mBSizes.at(grpSrc), MPI_PACKED, src, tag, MPI_COMM_WORLD, &(this->mRecvBRequests[packs].at(grpSrc)));
         }

         // Create send forward requests
         for(int id = 0; id < this->nFCpu(); ++id)
         {
            // Set MPI tag
            tag = FrameworkMacro::id();
            // Get CPU group index of local node
            grpMe = (*std::find(this->mFCpuGroup.begin(), this->mFCpuGroup.end(), tag));
            // Get destination index in CPU group
            grpDest = this->sendDest(id, grpMe, this->nFCpu());
            // Get global MPI destination rank
            dest = this->fCpu(grpDest);
            // initialise the Send request
            MPI_Send_init(this->mpFBuffers->at(grpDest), packs*this->mFSizes.at(grpDest), MPI_PACKED, dest, tag, MPI_COMM_WORLD, &(this->mSendFRequests[packs].at(grpDest)));
         }
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> void MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::cleanupRequests()
   {
      // Create iterator
      std::map<int, std::vector<MPI_Request> >::iterator it;

      // Free requests from Recv B
      for(it = this->mRecvBRequests.begin(); it != this->mRecvBRequests.end(); it++)
      {
         for(unsigned int i = 0; i < (*it).second.size(); ++i)
         {
            MPI_Request_free(&((*it).second.at(i)));
         }
      }

      // Free requests from Recv F
      for(it = this->mRecvFRequests.begin(); it != this->mRecvFRequests.end(); it++)
      {
         for(unsigned int i = 0; i < (*it).second.size(); ++i)
         {
            MPI_Request_free(&((*it).second.at(i)));
         }
      }

      // Free requests from Send B
      for(it = this->mSendBRequests.begin(); it != this->mSendBRequests.end(); it++)
      {
         for(unsigned int i = 0; i < (*it).second.size(); ++i)
         {
            MPI_Request_free(&((*it).second.at(i)));
         }
      }

      // Free requests from Send F
      for(it = this->mSendFRequests.begin(); it != this->mSendFRequests.end(); it++)
      {
         for(unsigned int i = 0; i < (*it).second.size(); ++i)
         {
            MPI_Request_free(&((*it).second.at(i)));
         }
      }
   }

}
}

#endif // MPICONVERTERBASE_HPP
