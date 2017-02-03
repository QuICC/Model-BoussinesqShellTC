/**
 * @file MpiConverterSendRecv.hpp
 * @brief Implementation of the MPI data converter with send/recv communication
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef MPICONVERTERSENDRECV_HPP
#define MPICONVERTERSENDRECV_HPP

// Configuration includes
//
#include "Debug/DebuggerMacro.h"
#include "Framework/FrameworkMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Communicators/Converters/MpiConverterBase.hpp"
#include "Communicators/Converters/MpiConverterTools.hpp"
#include "Resolutions/Resolution.hpp"
#include "Timers/StageTimer.hpp"

namespace QuICC {

namespace Parallel {

   /**
    * @brief Implementation of the MPI data converter with send/recv communication.
    */
   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> class MpiConverterSendRecv: public MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>
   {
      public:
         /**
          * @brief Constructor
          */
         MpiConverterSendRecv();

         /**
          * @brief Destructor
          */
         virtual ~MpiConverterSendRecv();

         /**
          * @brief Finish the setup of the converter
          */
         virtual void setup();

         /**
          * @brief Setup upcoming communication
          *
          * @param packs Number of packets in communication packing
          */
         virtual void setupCommunication(const int packs, const TransformDirection::Id direction);

         /**
          * @brief Start persistent send for forward transform
          */
         virtual void initiateForwardSend();

         /**
          * @brief Post persistent receive for forward transform
          */
         virtual void prepareForwardReceive();

         /**
          * @brief Start persistent send for backward transform
          */
         virtual void initiateBackwardSend();

         /**
          * @brief Post persistent receive for backward transform
          */
         virtual void prepareBackwardReceive();

      #ifdef QUICC_STORAGEPROFILE
         /**
         * @brief Do storage profiling
         */
         virtual void profileStorage() const;
      #endif // QUICC_STORAGEPROFILE
         
      protected:
         /**
          * @brief Send forward data 
          *
          * @param data Data to send
          */
         void sendFwd(const TFwdA &data);

         /**
          * @brief Send backward data
          *
          * @param data Data to send
          */
         void sendBwd(const TBwdB &data);

         /**
          * @brief Receive forward data
          *
          * @param rData Storage for received data
          */
         void receiveFwd(TFwdA &rData);

         /**
          * @brief Receive backward data
          *
          * @param rData Storage for received data
          */
         void receiveBwd(TBwdB &rData);

      private:
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
          * @brief Make sure forward buffer is available
          */
         void syncFwdBuffer();

         /**
          * @brief Make sure backward buffer is available
          */
         void syncBwdBuffer();

         /**
          * @brief Make sure forward buffer used by send is available
          */
         void syncFwdSendBuffer();

         /**
          * @brief Make sure backward buffer used by send is available
          */
         void syncBwdSendBuffer();

         /**
          * @brief Make sure forward buffer used by receive is available
          */
         void syncFwdRecvBuffer();

         /**
          * @brief Make sure backward buffer used by receive is available
          */
         void syncBwdRecvBuffer();

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
          * @brief Sending communication status
          */
         bool  mIsSending;

         /**
          * @brief Receiving communication status
          */
         bool  mIsReceiving;

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

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> inline MPI_Request * MpiConverterSendRecv<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::pRecvBRequests(const int size)
   {
      return &(this->mRecvBRequests.at(size).front());
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> inline MPI_Request * MpiConverterSendRecv<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::pRecvFRequests(const int size)
   {
      return &(this->mRecvFRequests.at(size).front());
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> inline MPI_Request * MpiConverterSendRecv<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::pSendBRequests(const int size)
   {
      return &(this->mSendBRequests.at(size).front());
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> inline MPI_Request * MpiConverterSendRecv<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::pSendFRequests(const int size)
   {
      return &(this->mSendFRequests.at(size).front());
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> MpiConverterSendRecv<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::MpiConverterSendRecv()
      : mIsSending(false), mIsReceiving(false), mActiveSend(0), mActiveReceive(0)
   {
      // Check that all dimensions match
      Debug::StaticAssert< (TFwdA::FieldDimension == TBwdA::FieldDimension) >();
      Debug::StaticAssert< (TBwdA::FieldDimension == TFwdB::FieldDimension) >();
      Debug::StaticAssert< (TFwdB::FieldDimension == TBwdB::FieldDimension) >();

      // Check that the data type is the same
      Debug::StaticTypeAssert<typename TFwdA::PointType , typename TBwdB::PointType>();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> MpiConverterSendRecv<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::~MpiConverterSendRecv()
   {
      // Cleanup the requests memory
      this->cleanupRequests();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterSendRecv<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::setup()
   {
      // Initialize positions
      this->resetFwdPositions();
      this->resetBwdPositions();

      // Initialise the active packs
      this->mActiveSend = 0;
      this->mActiveReceive = 0;
      this->mActiveDirection = TransformDirection::FORWARD;

      // setup the communication requests
      this->setupRequests();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterSendRecv<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::setupCommunication(const int packs, const TransformDirection::Id direction)
   {
      // Store the number of packs in active transfer
      this->mActiveSend = this->mPacks;
      this->mActiveReceive = this->mPacks;
      this->mActiveDirection = this->mDirection;

      // Store the number of packs in the next communication in direction
      this->mPacks = packs;
      this->mDirection = direction;
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterSendRecv<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::prepareForwardReceive()
   {
      // Don't do anything if the number of packs is zero
      if(this->mPacks > 0)
      {
         // Make sure the buffer is free
         this->syncFwdBuffer();

         // Make sure calls are posted at the right moment
         int flag;
         int ierr = MPI_Testall(this->nFCpu(), this->pRecvFRequests(this->mPacks), &flag, MPI_STATUSES_IGNORE);
         FrameworkMacro::check(ierr, 712);

         // Can only happen if something went really wrong
         if(!flag)
         {
            FrameworkMacro::abort(982+10*static_cast<int>(this->mTraId));
         }

         // Prepost the receive calls
         ierr = MPI_Startall(this->nFCpu(), this->pRecvFRequests(this->mPacks));
         FrameworkMacro::check(ierr, 722);
         this->resetBwdPositions();
         this->mIsReceiving = true;
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterSendRecv<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::initiateBackwardSend()
   {
      // Don't do anything if the number of packs is zero
      if(this->mPacks > 0)
      {
      	// Synchronize 
     	   FrameworkMacro::syncTransform(this->mTraId);

         // Make sure calls are posted at the right moment
         int flag;
         int ierr = MPI_Testall(this->nBCpu(), this->pSendBRequests(this->mPacks), &flag, MPI_STATUSES_IGNORE);
         FrameworkMacro::check(ierr, 713);

         // Can only happen if something went really wrong
         if(!flag)
         {
            FrameworkMacro::abort(983+10*static_cast<int>(this->mTraId));
         }

         // Post non blocking send calls 
         ierr = MPI_Startall(this->nBCpu(), this->pSendBRequests(this->mPacks));
         FrameworkMacro::check(ierr, 723);
         this->resetFwdPositions();
         this->mIsSending = true;
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterSendRecv<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::prepareBackwardReceive()
   {
      // Don't do anything if the number of packs is zero
      if(this->mPacks > 0)
      {
         // Make sure the buffer is free
         this->syncBwdBuffer();

         // Make sure calls are posted at the right moment
         int flag;
         int ierr = MPI_Testall(this->nBCpu(), this->pRecvBRequests(this->mPacks), &flag, MPI_STATUSES_IGNORE);
         FrameworkMacro::check(ierr, 714);

         // Can only happen if something went really wrong
         if(!flag)
         {
            FrameworkMacro::abort(984+10*static_cast<int>(this->mTraId));
         }

         // Prepost the receive calls
         ierr = MPI_Startall(this->nBCpu(), this->pRecvBRequests(this->mPacks));
         FrameworkMacro::check(ierr, 724);
         this->resetFwdPositions();
         this->mIsReceiving = true;
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterSendRecv<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::initiateForwardSend()
   {
      // Don't do anything if the number of packs is zero
      if(this->mPacks > 0)
      {
      	// Synchronize 
     	   FrameworkMacro::syncTransform(this->mTraId);

         // Make sure calls are posted at the right moment
         int flag;
         int ierr = MPI_Testall(this->nFCpu(), this->pSendFRequests(this->mPacks), &flag, MPI_STATUSES_IGNORE);
         FrameworkMacro::check(ierr, 715);

         // Can only happen if something went really wrong
         if(!flag)
         {
            FrameworkMacro::abort(985+10*static_cast<int>(this->mTraId));
         }

         // Post non blocking send calls 
         ierr = MPI_Startall(this->nFCpu(), this->pSendFRequests(this->mPacks));
         FrameworkMacro::check(ierr, 725);
         this->resetBwdPositions();
         this->mIsSending = true;
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterSendRecv<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::syncFwdBuffer()
   {
      if(this->mActiveDirection == TransformDirection::FORWARD && this->mDirection == TransformDirection::BACKWARD)
      {
         this->syncFwdRecvBuffer();

      } else if(this->mActiveDirection == TransformDirection::BACKWARD && this->mDirection == TransformDirection::FORWARD)
      {
         this->syncFwdSendBuffer();

      } else if(this->mActiveDirection == TransformDirection::FORWARD && this->mDirection == this->mActiveDirection)
      {
         this->syncFwdRecvBuffer();

      } else if(this->mActiveDirection == TransformDirection::BACKWARD && this->mDirection == this->mActiveDirection)
      {
         this->syncFwdSendBuffer();
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterSendRecv<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::syncBwdBuffer()
   {
      if(this->mActiveDirection == TransformDirection::FORWARD && this->mDirection == TransformDirection::BACKWARD)
      {
         this->syncBwdSendBuffer();

      } else if(this->mActiveDirection == TransformDirection::BACKWARD && this->mDirection == TransformDirection::FORWARD)
      {
         this->syncBwdRecvBuffer();

      } else if(this->mActiveDirection == TransformDirection::FORWARD && this->mDirection == this->mActiveDirection)
      {
         this->syncBwdSendBuffer();

      } else if(this->mActiveDirection == TransformDirection::BACKWARD && this->mDirection == this->mActiveDirection)
      {
         this->syncBwdRecvBuffer();
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterSendRecv<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::syncFwdRecvBuffer()
   {
      // Make sure previous communication has finished
      if(this->mActiveReceive > 0)
      {
         int flag;
         // (depending on MPI implementation the double test (test+wait) is required for the expected result)
         int ierr = MPI_Testall(this->nFCpu(), this->pRecvFRequests(this->mActiveReceive), &flag, MPI_STATUSES_IGNORE);
         FrameworkMacro::check(ierr, 731);

         // If not all are ready yet wait for completion
         if(! flag)
         {
            ierr = MPI_Waitall(this->nFCpu(), this->pRecvFRequests(this->mActiveReceive), MPI_STATUSES_IGNORE);
            FrameworkMacro::check(ierr, 741);
         }

         // Clear active packs
         this->mActiveReceive = 0;
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterSendRecv<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::syncBwdRecvBuffer()
   {
      // Make sure previous communication has finished
      if(this->mActiveReceive > 0)
      {
         int flag;
         // (depending on MPI implementation the double test (test+wait) is required for the expected result)
         int ierr = MPI_Testall(this->nBCpu(), this->pRecvBRequests(this->mActiveReceive), &flag, MPI_STATUSES_IGNORE);
        	FrameworkMacro::check(ierr, 732);

         // If not all are ready yet wait for completion
         if(! flag)
         {
            ierr = MPI_Waitall(this->nBCpu(), this->pRecvBRequests(this->mActiveReceive), MPI_STATUSES_IGNORE);
        	   FrameworkMacro::check(ierr, 742);
         }

         // Clear active packs
         this->mActiveReceive = 0;
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterSendRecv<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::syncFwdSendBuffer()
   {
      // Make sure previous communication has finished
      if(this->mActiveSend > 0)
      {
         int flag;
         // (depending on MPI implementation the double test (test+wait) is required for the expected result)
         int ierr = MPI_Testall(this->nBCpu(), this->pSendFRequests(this->mActiveSend), &flag, MPI_STATUSES_IGNORE);
        	FrameworkMacro::check(ierr, 733);

         // If not all are ready yet wait for completion
         if(! flag)
         {
            ierr = MPI_Waitall(this->nBCpu(), this->pSendFRequests(this->mActiveSend), MPI_STATUSES_IGNORE);
        	   FrameworkMacro::check(ierr, 743);
         }

         // Clear active packs
         this->mActiveSend = 0;
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterSendRecv<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::syncBwdSendBuffer()
   {
      // Make sure previous communication has finished
      if(this->mActiveSend > 0)
      {
         int flag;
         // (depending on MPI implementation the double test (test+wait) is required for the expected result)
         int ierr = MPI_Testall(this->nFCpu(), this->pSendBRequests(this->mActiveSend), &flag, MPI_STATUSES_IGNORE);
        	FrameworkMacro::check(ierr, 734);

         // If not all are ready yet wait for completion
         if(! flag)
         {
            ierr = MPI_Waitall(this->nFCpu(), this->pSendBRequests(this->mActiveSend), MPI_STATUSES_IGNORE);
        	   FrameworkMacro::check(ierr, 744);
         }

         // Clear active packs
         this->mActiveSend = 0;
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterSendRecv<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::sendFwd(const TFwdA& data)
   {
      // Should be in backward direction
      assert(this->mDirection == TransformDirection::BACKWARD);

      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::FWDSENDWAIT);

      // Make sure the buffer is free
      this->syncFwdBuffer();

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::FWDSENDWAIT);

      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::FWDSENDCONV);

      // Pack data into send buffer
      for(int id = 0; id < this->nFCpu(); ++id)
      {
         #if defined QUICC_MPIPACK_MANUAL
            MpiConverterTools<TFwdA::FieldDimension>::template pack<typename TFwdA::PointType>(this->mspFBuffers->at(id), this->mspFBuffers->pos(id), data, this->mFTypes.at(id));
         #else
            int ierr = MPI_Pack(const_cast<typename TFwdA::PointType *>(data.data().data()), 1, this->mFTypes.at(id), this->mspFBuffers->at(id), this->sizeFPacket(id), &(this->mspFBuffers->pos(id)), FrameworkMacro::transformComm(this->mTraId));
            FrameworkMacro::check(ierr, 761);
         #endif //defined QUICC_MPIPACK_MANUAL
      }

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::FWDSENDCONV);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterSendRecv<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::sendBwd(const TBwdB& data)
   {
      // Should be in backward direction
      assert(this->mDirection == TransformDirection::FORWARD);

      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWDSENDWAIT);

      // Make sure the buffer is free
      this->syncBwdBuffer();

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWDSENDWAIT);

      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWDSENDCONV);

      // Pack data into send buffer
      for(int id = 0; id < this->nBCpu(); ++id)
      {
         #if defined QUICC_MPIPACK_MANUAL
            MpiConverterTools<TBwdB::FieldDimension>::template pack<typename TBwdB::PointType>(this->mspBBuffers->at(id), this->mspBBuffers->pos(id), data, this->mBTypes.at(id));
         #else
            int ierr = MPI_Pack(const_cast<typename TBwdB::PointType *>(data.data().data()), 1, this->mBTypes.at(id), this->mspBBuffers->at(id), this->sizeBPacket(id), &(this->mspBBuffers->pos(id)), FrameworkMacro::transformComm(this->mTraId));
            FrameworkMacro::check(ierr, 762);
         #endif //defined QUICC_MPIPACK_MANUAL
      }

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWDSENDCONV);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterSendRecv<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::receiveFwd(TFwdA &rData)
   {
      // Should be in backward direction
      assert(this->mDirection == TransformDirection::FORWARD);

      // Communication interface is receiving the data 
      if(this->mIsReceiving)
      {
         // Number of receive calls in total
         int keepWaiting = this->nFCpu();
         int count = 0;
         ArrayI   idx = ArrayI::Zero(this->nFCpu());

         // Wait until everything has been received
         while(keepWaiting != 0)
         {
            // Start detailed profiler
            DetailedProfilerMacro_start(ProfilerMacro::FWDRECVWAIT);

            // Wait for some of the requests to finish
            #ifdef QUICC_DEBUG
               MPI_Status stats[this->nFCpu()];
               int ierr = MPI_Waitsome(this->nFCpu(), this->pRecvFRequests(this->mPacks), &count, idx.data(), stats);
               FrameworkMacro::check(ierr, 771);
               DebuggerMacro_msg("Received FWD packs", 5);
            #else
               int ierr = MPI_Waitsome(this->nFCpu(), this->pRecvFRequests(this->mPacks), &count, idx.data(), MPI_STATUSES_IGNORE);
               FrameworkMacro::check(ierr, 771);
            #endif //QUICC_DEBUG

            // Stop detailed profiler
            DetailedProfilerMacro_stop(ProfilerMacro::FWDRECVWAIT);

            // Start detailed profiler
            DetailedProfilerMacro_start(ProfilerMacro::FWDRECVCONV);

            // Unpack already received data from receive buffer
            for(int id = 0; id < count; ++id)
            {
               DebuggerMacro_showValue("Tag: ", 6, stats[id].MPI_TAG);
               DebuggerMacro_showValue("-> From: ", 6, stats[id].MPI_SOURCE);
               int pos = idx(id);

               #if defined QUICC_MPIPACK_MANUAL
                  MpiConverterTools<TFwdA::FieldDimension>::template unpack<typename TFwdA::PointType>(rData, this->mFTypes.at(pos), this->mspFBuffers->at(pos), this->mspFBuffers->pos(pos));
               #else
                  int ierr = MPI_Unpack(this->mspFBuffers->at(pos), this->sizeFPacket(pos), &(this->mspFBuffers->pos(pos)), rData.rData().data(), 1, this->mFTypes.at(pos), FrameworkMacro::transformComm(this->mTraId));
                  FrameworkMacro::check(ierr, 763);
               #endif //defined QUICC_MPIPACK_MANUAL
            }

            // Stop detailed profiler
            DetailedProfilerMacro_stop(ProfilerMacro::FWDRECVCONV);

            // Update the number of missing receives
            keepWaiting -= count;
         }

         // Reset communication status
         this->mIsReceiving = false;

      // Data is here and just need to be unpacked
      } else
      {
         // Start detailed profiler
         DetailedProfilerMacro_start(ProfilerMacro::FWDRECVCONV);

         // Unpack data from receive buffer
         for(int id = 0; id < this->nFCpu(); ++id)
         {
            DebuggerMacro_msg("Unpacking FWD packs", 5);

            #if defined QUICC_MPIPACK_MANUAL
               MpiConverterTools<TFwdA::FieldDimension>::template unpack<typename TFwdA::PointType>(rData, this->mFTypes.at(id), this->mspFBuffers->at(id), this->mspFBuffers->pos(id));
            #else
               int ierr = MPI_Unpack(this->mspFBuffers->at(id), this->sizeFPacket(id), &(this->mspFBuffers->pos(id)), rData.rData().data(), 1, this->mFTypes.at(id), FrameworkMacro::transformComm(this->mTraId));
               FrameworkMacro::check(ierr, 764);
            #endif //defined QUICC_MPIPACK_MANUAL

         }

         // Stop detailed profiler
         DetailedProfilerMacro_stop(ProfilerMacro::FWDRECVCONV);
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterSendRecv<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::receiveBwd(TBwdB &rData)
   {
      // Should be in backward direction
      assert(this->mDirection == TransformDirection::BACKWARD);

      // Communication interface is receiving the data 
      if(this->mIsReceiving)
      {
         // Number of receive calls in total
         int keepWaiting = this->nBCpu();
         int count = 0;
         ArrayI   idx = ArrayI::Zero(this->nBCpu());

         // Wait until everything has been received
         while(keepWaiting != 0)
         {
            // Start detailed profiler
            DetailedProfilerMacro_start(ProfilerMacro::BWDRECVWAIT);

            // Wait for some of the requests to finish
            #ifdef QUICC_DEBUG
               MPI_Status stats[this->nBCpu()];
               int ierr = MPI_Waitsome(this->nBCpu(), this->pRecvBRequests(this->mPacks), &count, idx.data(), stats);
               FrameworkMacro::check(ierr, 772);
               DebuggerMacro_msg("Received BWD packs", 5);
            #else 
               int ierr = MPI_Waitsome(this->nBCpu(), this->pRecvBRequests(this->mPacks), &count, idx.data(), MPI_STATUSES_IGNORE);
               FrameworkMacro::check(ierr, 772);
            #endif //QUICC_DEBUG

            // Stop detailed profiler
            DetailedProfilerMacro_stop(ProfilerMacro::BWDRECVWAIT);

            // Start detailed profiler
            DetailedProfilerMacro_start(ProfilerMacro::BWDRECVCONV);

            // Unpack already received data from receive buffer
            for(int id = 0; id < count; ++id)
            {
               DebuggerMacro_showValue("Tag: ", 6, stats[id].MPI_TAG);
               DebuggerMacro_showValue("-> From: ", 6, stats[id].MPI_SOURCE);
               int pos = idx(id);

               #if defined QUICC_MPIPACK_MANUAL
                  MpiConverterTools<TBwdB::FieldDimension>::template unpack<typename TBwdB::PointType>(rData, this->mBTypes.at(pos), this->mspBBuffers->at(pos), this->mspBBuffers->pos(pos));
               #else
                  int ierr = MPI_Unpack(this->mspBBuffers->at(pos), this->sizeBPacket(pos), &(this->mspBBuffers->pos(pos)), rData.rData().data(), 1, this->mBTypes.at(pos), FrameworkMacro::transformComm(this->mTraId));
                  FrameworkMacro::check(ierr, 765);
               #endif //defined QUICC_MPIPACK_MANUAL
            }

            // Stop detailed profiler
            DetailedProfilerMacro_stop(ProfilerMacro::BWDRECVCONV);

            // Update the number of missing receives
            keepWaiting -= count;
         }

         // Reset communication status
         this->mIsReceiving = false;

      // Data is here and just need to be unpacked
      } else
      {
         // Start detailed profiler
         DetailedProfilerMacro_start(ProfilerMacro::BWDRECVCONV);

         // Unpack data from receive buffer
         for(int id = 0; id < this->nBCpu(); ++id)
         {
            DebuggerMacro_msg("Unpacking BWD packs", 5);

            #if defined QUICC_MPIPACK_MANUAL
               MpiConverterTools<TBwdB::FieldDimension>::template unpack<typename TBwdB::PointType>(rData, this->mBTypes.at(id), this->mspBBuffers->at(id), this->mspBBuffers->pos(id));
            #else
               int ierr = MPI_Unpack(this->mspBBuffers->at(id), this->sizeBPacket(id), &(this->mspBBuffers->pos(id)), rData.rData().data(), 1, this->mBTypes.at(id), FrameworkMacro::transformComm(this->mTraId));
               FrameworkMacro::check(ierr, 766);
            #endif //defined QUICC_MPIPACK_MANUAL
         }

         // Stop detailed profiler
         DetailedProfilerMacro_stop(ProfilerMacro::BWDRECVCONV);
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> int MpiConverterSendRecv<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::sendDest(const int id, const int ref, const int size) const
   {
      // Create send ring
      return ((id + 1 + ref) % size);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> int MpiConverterSendRecv<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::recvSrc(const int id, const int ref, const int size) const
   {
      // Create recv ring
      return ((size - 1 - id + ref) % size);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterSendRecv<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::setupRequests()
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
            #if defined QUICC_MPIPACK_MANUAL
               ierr = MPI_Recv_init(this->mspFBuffers->at(grpSrc), packs*this->mFSizes.at(grpSrc), MpiTypes::template type<typename TFwdA::PointType>(), src, tag, FrameworkMacro::transformComm(this->mTraId), &(this->mRecvFRequests.at(packs).at(grpSrc)));
               FrameworkMacro::check(ierr, 981);
            #else
               ierr = MPI_Recv_init(this->mspFBuffers->at(grpSrc), packs*this->mFSizes.at(grpSrc), MPI_PACKED, src, tag, FrameworkMacro::transformComm(this->mTraId), &(this->mRecvFRequests.at(packs).at(grpSrc)));
               FrameworkMacro::check(ierr, 981);
            #endif //defined QUICC_MPIPACK_MANUAL
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
            #if defined QUICC_MPIPACK_MANUAL
               ierr = MPI_Send_init(this->mspBBuffers->at(grpDest), packs*this->mBSizes.at(grpDest), MpiTypes::template type<typename TBwdB::PointType>(), dest, tag, FrameworkMacro::transformComm(this->mTraId), &(this->mSendBRequests.at(packs).at(grpDest)));
               FrameworkMacro::check(ierr, 982);
            #else
               ierr = MPI_Send_init(this->mspBBuffers->at(grpDest), packs*this->mBSizes.at(grpDest), MPI_PACKED, dest, tag, FrameworkMacro::transformComm(this->mTraId), &(this->mSendBRequests.at(packs).at(grpDest)));
               FrameworkMacro::check(ierr, 982);
            #endif //defined QUICC_MPIPACK_MANUAL
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
            #if defined QUICC_MPIPACK_MANUAL
               ierr = MPI_Recv_init(this->mspBBuffers->at(grpSrc), packs*this->mBSizes.at(grpSrc), MpiTypes::template type<typename TBwdB::PointType>(), src, tag, FrameworkMacro::transformComm(this->mTraId), &(this->mRecvBRequests.at(packs).at(grpSrc)));
            #else
               ierr = MPI_Recv_init(this->mspBBuffers->at(grpSrc), packs*this->mBSizes.at(grpSrc), MPI_PACKED, src, tag, FrameworkMacro::transformComm(this->mTraId), &(this->mRecvBRequests.at(packs).at(grpSrc)));
            #endif //defined QUICC_MPIPACK_MANUAL
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
            #if defined QUICC_MPIPACK_MANUAL
               ierr = MPI_Send_init(this->mspFBuffers->at(grpDest), packs*this->mFSizes.at(grpDest), MpiTypes::template type<typename TFwdA::PointType>(), dest, tag, FrameworkMacro::transformComm(this->mTraId), &(this->mSendFRequests.at(packs).at(grpDest)));
            #else
               ierr = MPI_Send_init(this->mspFBuffers->at(grpDest), packs*this->mFSizes.at(grpDest), MPI_PACKED, dest, tag, FrameworkMacro::transformComm(this->mTraId), &(this->mSendFRequests.at(packs).at(grpDest)));
               FrameworkMacro::check(ierr, 984);
            #endif //defined QUICC_MPIPACK_MANUAL
         }
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterSendRecv<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::cleanupRequests()
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

#ifdef QUICC_STORAGEPROFILE
   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterSendRecv<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::profileStorage() const
   {
      MHDFloat memTypes = 0.0;
      MHDFloat memComm = 0.0;

      // General communication storage
      memComm += 4.0*(this->mFCommSizes.size());
      memComm += 4.0*(this->mBCommSizes.size());

      // Requests communication storage
      memComm += 0.0;

      // MPI datatypes storage
      memTypes += 0.0;

      StorageProfilerMacro_update(StorageProfilerMacro::MPI, memTypes + memComm);
      #ifdef QUICC_STORAGEPROFILER_DETAILED
         StorageProfilerMacro_update(StorageProfilerMacro::MPITYPES, memTypes);

         StorageProfilerMacro_update(StorageProfilerMacro::MPICOMM, memComm);
      #endif // QUICC_STORAGEPROFILER_DETAILED
   }
#endif // QUICC_STORAGEPROFILE

}
}

#endif // MPICONVERTERSENDRECV_HPP
