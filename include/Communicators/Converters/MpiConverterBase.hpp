/**
 * @file MpiConverterBase.hpp
 * @brief Templated implementation of the base of a MPI data converter 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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
#include "Enums/TransformDirection.hpp"
#include "Resolutions/Resolution.hpp"
#include "Communicators/Converters/IConverter.hpp"
#include "Communicators/CommunicationBuffer.hpp"

namespace GeoMHDiSCC {

namespace Parallel {

   /**
    * @brief Templated implementation of the base of a MPI data converter.
    */
   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> class MpiConverterBase: public IConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>
   {
      public:
         #if defined GEOMHDISCC_MPIPACK_MANUAL
            typedef SharedPtrMacro<CommunicationBuffer<typename TFwdA::PointType> > SharedFwdBufferType;
            typedef SharedPtrMacro<CommunicationBuffer<typename TBwdB::PointType> > SharedBwdBufferType;
         #else
            typedef SharedPtrMacro<CommunicationBuffer<char> > SharedFwdBufferType;
            typedef SharedPtrMacro<CommunicationBuffer<char> > SharedBwdBufferType;
         #endif //defined GEOMHDISCC_MPIPACK_MANUAL

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
         void setBuffers(SharedFwdBufferType spFwd, SharedBwdBufferType spBwd);

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
          * @brief Reset Fwe buffer positions
          */
         void resetFwdPositions();

         /**
          * @brief Reset Bwd buffer positions
          */
         void resetBwdPositions();

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
          * @brief Setup the MPI communication requests requests
          */
         void setupRequests();

         /**
          * @brief Cleanup the MPI communication requests
          */
         void cleanupRequests();

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
          * @brief Get MPI rank of CPU from forward CPU group
          *
          * @param id CPU group id
          */
         int fCpu(const int id) const;

         /**
          * @brief Get MPI rank of CPU from backward CPU group
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
          * @brief Communication packs counter
          */
         int mPacks;

         /**
          * @brief Direction of operations
          */
         TransformDirection::Id   mDirection;

         /**
          * @brief Transform ID
          */
         Dimensions::Transform::Id mTraId;

         /**
          * @brief Direction of active operation
          */
         TransformDirection::Id   mActiveDirection;

         /**
          * @brief The number of packs in the "previous/active" send operations
          */
         int   mActiveSend;

         /**
          * @brief The number of packs in the "previous/active" receive operations
          */
         int   mActiveReceive;

         /**
          * @brief Storage for the receive position pointers
          */
         //std::vector<int>  mRecvPositions;

         /**
          * @brief Storage for the send position pointers
          */
         //std::vector<int>  mSendPositions;

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
         SharedFwdBufferType mspFBuffers;

         /**
          * @brief Backward communication buffer pointer
          */
         SharedBwdBufferType mspBBuffers;

         /**
          * @brief List of the forward buffer sizes
          */
         std::vector<int>  mFSizes;

         /**
          * @brief List of the backward buffer sizes
          */
         std::vector<int>  mBSizes;

         /**
          * @brief Possible forward transform packs
          */
         ArrayI   mForwardPacks;

         /**
          * @brief Possible backward transform packs
          */
         ArrayI   mBackwardPacks;

      private:
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
   };

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> inline const std::vector<int>& MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::fwdSizes() const
   {
      return this->mFSizes;
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> inline const std::vector<int>& MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::bwdSizes() const
   {
      return this->mBSizes;
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::setBuffers(typename MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::SharedFwdBufferType spFwd, typename MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::SharedBwdBufferType spBwd)
   {
      // Set the forward buffers
      this->mspFBuffers = spFwd;

      // Set the backward buffers
      this->mspBBuffers = spBwd;
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> inline int MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::nFCpu() const
   {
      return this->mFCpuGroup.size();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> inline int MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::nBCpu() const
   {
      return this->mBCpuGroup.size();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> inline int MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::fCpu(const int id) const
   {
      return this->mFCpuGroup.at(id);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> inline int MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::bCpu(const int id) const
   {
      return this->mBCpuGroup.at(id);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> inline int MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::sizeFPacket(const int id) const
   {
      return this->mPacks*this->mFSizes.at(id);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> inline int MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::sizeBPacket(const int id) const
   {
      return this->mPacks*this->mBSizes.at(id);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> inline MPI_Request * MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::pRecvBRequests(const int size)
   {
      return &(this->mRecvBRequests.at(size).front());
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> inline MPI_Request * MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::pRecvFRequests(const int size)
   {
      return &(this->mRecvFRequests.at(size).front());
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> inline MPI_Request * MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::pSendBRequests(const int size)
   {
      return &(this->mSendBRequests.at(size).front());
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> inline MPI_Request * MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::pSendFRequests(const int size)
   {
      return &(this->mSendFRequests.at(size).front());
   }
      
   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::MpiConverterBase()
      : IConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>(), mIsSending(false), mIsReceiving(false), mPacks(0), mActiveSend(0), mActiveReceive(0)
   {
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::~MpiConverterBase()
   {
      // Cleanup the requests memory
      this->cleanupRequests();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::resetFwdPositions()
   {
      this->mspFBuffers->resetPositions();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::resetBwdPositions()
   {
      this->mspBBuffers->resetPositions();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> int MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::sendDest(const int id, const int ref, const int size) const
   {
      // Create send ring
      return ((id + 1 + ref) % size);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> int MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::recvSrc(const int id, const int ref, const int size) const
   {
      // Create recv ring
      return ((size - 1 - id + ref) % size);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::setupRequests()
   {
      // Storage for global location flags
      int dest;
      int src;
      int tag;
      // Storage for CPU group location flags
      int grpMe;
      int grpDest;
      int grpSrc;

      // Shift tag to produce unique tag in 2D distribution
      // (Probably not required anymore but doesn't harm)
      int tagShift = 0;
      if(this->mTraId == Dimensions::Transform::TRA1D)
      {
         tagShift = 0;
      } else if(this->mTraId == Dimensions::Transform::TRA2D)
      {
         tagShift = FrameworkMacro::nCpu();
      } else
      {
         FrameworkMacro::abort(991);
      }

      // MPI error code
      int ierr;

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
            this->mRecvFRequests.at(packs).push_back(MPI_REQUEST_NULL);
         }

         // Initialise send backward with empty requests
         this->mSendBRequests.insert(std::make_pair<int, std::vector<MPI_Request> >(packs, std::vector<MPI_Request>()));
         for(int id = 0; id < this->nBCpu(); ++id)
         {
            this->mSendBRequests.at(packs).push_back(MPI_REQUEST_NULL);
         }

         // Create receive forward requests
         for(int id = 0; id < this->nFCpu(); ++id)
         {
            // Get CPU group index of local node
            grpMe = (*std::find(this->mFCpuGroup.begin(), this->mFCpuGroup.end(), FrameworkMacro::transformId(this->mTraId)));

            // Get source index in CPU group
            grpSrc = this->recvSrc(id, grpMe, this->nFCpu());

            // Get source MPI rank in group
            src = this->fCpu(grpSrc);

            // Set shifted MPI tag to make it unique
            tag = src + tagShift;

            //Safety asserts
            assert(static_cast<size_t>(grpSrc) < this->mFSizes.size());
            assert(static_cast<size_t>(grpSrc) < this->mRecvFRequests.at(packs).size());

            // initialise the Recv request
            #if defined GEOMHDISCC_MPIPACK_MANUAL
               ierr = MPI_Recv_init(this->mspFBuffers->at(grpSrc), packs*this->mFSizes.at(grpSrc), MpiTypes::template type<typename TFwdA::PointType>(), src, tag, FrameworkMacro::transformComm(this->mTraId), &(this->mRecvFRequests.at(packs).at(grpSrc)));
               FrameworkMacro::check(ierr, 981);
            #else
               ierr = MPI_Recv_init(this->mspFBuffers->at(grpSrc), packs*this->mFSizes.at(grpSrc), MPI_PACKED, src, tag, FrameworkMacro::transformComm(this->mTraId), &(this->mRecvFRequests.at(packs).at(grpSrc)));
               FrameworkMacro::check(ierr, 981);
            #endif //defined GEOMHDISCC_MPIPACK_MANUAL
         }

         // Create send backward requests
         for(int id = 0; id < this->nBCpu(); ++id)
         {
            // Get CPU group index of local node
            grpMe = (*std::find(this->mBCpuGroup.begin(), this->mBCpuGroup.end(), FrameworkMacro::transformId(this->mTraId)));

            // Set shifted MPI tag to make it unique
            tag = FrameworkMacro::transformId(this->mTraId) + tagShift;

            // Get destination index in CPU group
            grpDest = this->sendDest(id, grpMe, this->nBCpu());

            // Get destination MPI rank in group
            dest = this->bCpu(grpDest);

            //Safety asserts
            assert(static_cast<size_t>(grpDest) < this->mBSizes.size());
            assert(static_cast<size_t>(grpDest) < this->mSendBRequests.at(packs).size());

            // initialise the Send request
            #if defined GEOMHDISCC_MPIPACK_MANUAL
               ierr = MPI_Send_init(this->mspBBuffers->at(grpDest), packs*this->mBSizes.at(grpDest), MpiTypes::template type<typename TBwdB::PointType>(), dest, tag, FrameworkMacro::transformComm(this->mTraId), &(this->mSendBRequests.at(packs).at(grpDest)));
               FrameworkMacro::check(ierr, 982);
            #else
               ierr = MPI_Send_init(this->mspBBuffers->at(grpDest), packs*this->mBSizes.at(grpDest), MPI_PACKED, dest, tag, FrameworkMacro::transformComm(this->mTraId), &(this->mSendBRequests.at(packs).at(grpDest)));
               FrameworkMacro::check(ierr, 982);
            #endif //defined GEOMHDISCC_MPIPACK_MANUAL
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
            this->mRecvBRequests.at(packs).push_back(MPI_REQUEST_NULL);
         }

         // Initialise send forward with empty requests
         this->mSendFRequests.insert(std::make_pair<int, std::vector<MPI_Request> >(packs, std::vector<MPI_Request>()));
         for(int id = 0; id < this->nFCpu(); ++id)
         {
            this->mSendFRequests.at(packs).push_back(MPI_REQUEST_NULL);
         }

         // Create receive backward requests
         for(int id = 0; id < this->nBCpu(); ++id)
         {
            // Get CPU group index of local node
            grpMe = (*std::find(this->mBCpuGroup.begin(), this->mBCpuGroup.end(), FrameworkMacro::transformId(this->mTraId)));

            // Get source index in CPU group
            grpSrc = this->recvSrc(id, grpMe, this->nBCpu());

            // Get source MPI rank in group
            src = this->bCpu(grpSrc);

            // Set shifted MPI tag to make it unique
            tag = src + tagShift;

            // initialise the Recv request
            #if defined GEOMHDISCC_MPIPACK_MANUAL
               ierr = MPI_Recv_init(this->mspBBuffers->at(grpSrc), packs*this->mBSizes.at(grpSrc), MpiTypes::template type<typename TBwdB::PointType>(), src, tag, FrameworkMacro::transformComm(this->mTraId), &(this->mRecvBRequests.at(packs).at(grpSrc)));
            #else
               ierr = MPI_Recv_init(this->mspBBuffers->at(grpSrc), packs*this->mBSizes.at(grpSrc), MPI_PACKED, src, tag, FrameworkMacro::transformComm(this->mTraId), &(this->mRecvBRequests.at(packs).at(grpSrc)));
            #endif //defined GEOMHDISCC_MPIPACK_MANUAL
            FrameworkMacro::check(ierr, 983);
         }

         // Create send forward requests
         for(int id = 0; id < this->nFCpu(); ++id)
         {
            // Get CPU group index of local node
            grpMe = (*std::find(this->mFCpuGroup.begin(), this->mFCpuGroup.end(), FrameworkMacro::transformId(this->mTraId)));

            // Set shifted MPI tag to make it unique
            tag = FrameworkMacro::transformId(this->mTraId) + tagShift;

            // Get destination index in CPU group
            grpDest = this->sendDest(id, grpMe, this->nFCpu());

            // Get destination MPI rank in group
            dest = this->fCpu(grpDest);

            // initialise the Send request
            #if defined GEOMHDISCC_MPIPACK_MANUAL
               ierr = MPI_Send_init(this->mspFBuffers->at(grpDest), packs*this->mFSizes.at(grpDest), MpiTypes::template type<typename TFwdA::PointType>(), dest, tag, FrameworkMacro::transformComm(this->mTraId), &(this->mSendFRequests.at(packs).at(grpDest)));
            #else
               ierr = MPI_Send_init(this->mspFBuffers->at(grpDest), packs*this->mFSizes.at(grpDest), MPI_PACKED, dest, tag, FrameworkMacro::transformComm(this->mTraId), &(this->mSendFRequests.at(packs).at(grpDest)));
               FrameworkMacro::check(ierr, 984);
            #endif //defined GEOMHDISCC_MPIPACK_MANUAL
         }
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::cleanupRequests()
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
