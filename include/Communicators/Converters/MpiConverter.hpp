/**
 * @file MpiConverter.hpp
 * @brief Implementation of the FDSH MPI data converter 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef MPICONVERTER_HPP
#define MPICONVERTER_HPP

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
#include "StorageProviders/StoragePairProviderMacro.h"
#include "Resolutions/Resolution.hpp"
#include "Timers/StageTimer.hpp"

//#define GEOMHDISCC_MPICOMM_SENDRECV
//#define GEOMHDISCC_MPICOMM_ALLTOALL

namespace GeoMHDiSCC {

namespace Parallel {

   /**
    * @brief Implementation of the FDSH MPI data converter.
    */
   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> class MpiConverter: public MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>
   {
      public:
         /**
          * @brief Constructor
          */
         MpiConverter();

         /**
          * @brief Destructor
          */
         virtual ~MpiConverter();

         /**
          * @brief Initialise the packs
          *
          * @param spRes      Shared Resolution
          * @param fwdDim     Dimension index for forward transform
          * @param fwdTmp     TFwdA temporary
          * @param bwdTmp     TBwdB temporary
          * @param fwdPacks   Array of possible pack sizes for forward transform
          * @param bwdPacks   Array of possible pack sizes for backward transform
          */
         void init(SharedResolution spRes, const Dimensions::Transform::Id fwdDim, TFwdA& fwdTmp, TBwdB& bwdTmp, const ArrayI& fwdPacks, const ArrayI& bwdPacks);

         /**
          * @brief Finish the setup of the converter
          */
         virtual void setup();

         /**
          * @brief Convert data from TFwdA to TBwdB
          *
          * @param in      Input forward data
          * @param storage Storage provider
          */
         virtual void convertFwd(const TFwdA &in, StoragePairProviderMacro<TFwdB, TBwdB> &storage);

         /**
          * @brief Convert data from TBwdB to TFwdA
          *
          * @param in      Input backward data
          * @param storage Storage provider
          */
         virtual void convertBwd(const TBwdB &in, StoragePairProviderMacro<TFwdA, TBwdA> &storage);

         /**
          * @brief Get the converted data from TBwdA to TFwdB conversion
          *
          * @param storage Storage provider
          */
         virtual TFwdA& getFwd(StoragePairProviderMacro<TFwdA, TBwdA> &storage);

         /**
          * @brief Get the converted data from TFwdB to TBwdA conversion
          *
          * @param storage Storage provider
          */
         virtual TBwdB& getBwd(StoragePairProviderMacro<TFwdB, TBwdB> &storage);

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

      #ifdef GEOMHDISCC_STORAGEPROFILE
         /**
         * @brief Do storage profiling
         */
         virtual void profileStorage() const;
      #endif // GEOMHDISCC_STORAGEPROFILE
         
      protected:

      private:
         /**
          * @brief Initialise the datatypes
          *
          * @param spRes   Shared Resolution
          * @param fTmp    TFwdA temporary
          * @param bTmp    TBwdB temporary
          */
         void initTypes(SharedResolution spRes, TFwdA& fTmp, TBwdB& bTmp);

         /**
          * @brief Initialise the sizes and CPU lists
          */
         void initLists();

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
          * @brief Cleanup the data types
          */
         void cleanupTypes();

         #if defined GEOMHDISCC_MPIPACK_MANUAL
            /**
             * @brief Storage for the forward datatypes
             */
            std::vector<std::vector<typename MpiConverterTools<TFwdA::FieldDimension>::Coordinate> >  mFTypes;

            /**
             * @brief Storage for the backward datatypes
             */
            std::vector<std::vector<typename MpiConverterTools<TBwdB::FieldDimension>::Coordinate> > mBTypes;

         #else
            /**
             * @brief Storage for the forward datatypes
             */
            std::vector<MPI_Datatype>  mFTypes;

            /**
             * @brief Storage for the backward datatypes
             */
            std::vector<MPI_Datatype> mBTypes;
         #endif //defined GEOMHDISCC_MPIPACK_MANUAL
   };

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::MpiConverter()
   {
      // Check that all dimensions match
      Debug::StaticAssert< (TFwdA::FieldDimension == TBwdA::FieldDimension) >();
      Debug::StaticAssert< (TBwdA::FieldDimension == TFwdB::FieldDimension) >();
      Debug::StaticAssert< (TFwdB::FieldDimension == TBwdB::FieldDimension) >();

      // Check that the data type is the same
      Debug::StaticTypeAssert<typename TFwdA::PointType , typename TBwdB::PointType>();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::~MpiConverter()
   {
      // Cleanup the types memory
      this->cleanupTypes();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::convertFwd(const TFwdA &in, StoragePairProviderMacro<TFwdB, TBwdB>  &storage)
   {
      // Send the data
      this->sendFwd(in);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::convertBwd(const TBwdB &in, StoragePairProviderMacro<TFwdA, TBwdA>  &storage)
   {
      // Send the data
      this->sendBwd(in);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> TFwdA& MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::getFwd(StoragePairProviderMacro<TFwdA, TBwdA>  &storage)
   {
      // Get storage for output value 
      TFwdA &rOut = storage.provideFwd();

      // Receive converted data
      this->receiveFwd(rOut);

      return rOut;
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> TBwdB& MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::getBwd(StoragePairProviderMacro<TFwdB, TBwdB>  &storage)
   {
      // Get storage for output value 
      TBwdB &rOut = storage.provideBwd();

      // Receive converted data
      this->receiveBwd(rOut);

      return rOut;
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::setup()
   {
      // Initialize positions
      this->resetFwdPositions();
      this->resetBwdPositions();

      // setup the communication requests
      this->setupRequests();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::setupCommunication(const int packs, const TransformDirection::Id direction)
   {
      // Store the number of packs in active transfer
      this->mActiveSend = this->mPacks;
      this->mActiveReceive = this->mPacks;
      this->mActiveDirection = this->mDirection;

      // Store the number of packs in the next communication in direction
      this->mPacks = packs;
      this->mDirection = direction;
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::prepareForwardReceive()
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

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::initiateBackwardSend()
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

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::prepareBackwardReceive()
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

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::initiateForwardSend()
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

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::init(SharedResolution spRes, const Dimensions::Transform::Id fwdDim, TFwdA &fwdTmp, TBwdB &bwdTmp, const ArrayI& fwdPacks, const ArrayI& bwdPacks)
   {
      StageTimer stage;
      stage.start("creating MPI datatypes",1);

      // Store the possible pack sizes
      this->mForwardPacks = fwdPacks;
      this->mBackwardPacks = bwdPacks;

      // Initialise the active packs
      this->mActiveSend = 0;
      this->mActiveReceive = 0;
      this->mActiveDirection = TransformDirection::FORWARD;

      // Store Transform ID
      this->mTraId = fwdDim;
      
      // Create index converter
      this->mspIdxConv = SharedPtrMacro<TIdx>(new TIdx(spRes, this->mTraId));

      // initialise the data types
      this->initTypes(spRes, fwdTmp, bwdTmp);

      stage.done();
      stage.start("cleaning empty datatypes",1);

      // initialise the size and CPU lists
      this->initLists();

      stage.done();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::initTypes(SharedResolution spRes, TFwdA &fTmp, TBwdB &bTmp)
   {
      std::map<typename MpiConverterTools<TBwdB::FieldDimension>::Coordinate, typename MpiConverterTools<TBwdB::FieldDimension>::Coordinate> localFwdMap;
      MpiConverterTools<TFwdA::FieldDimension>::buildLocalFwdMap(localFwdMap, spRes, this->mTraId);
      std::map<typename MpiConverterTools<TBwdB::FieldDimension>::Coordinate, typename MpiConverterTools<TBwdB::FieldDimension>::Coordinate> localBwdMap;
      MpiConverterTools<TBwdB::FieldDimension>::buildLocalBwdMap(localBwdMap, spRes, this->mTraId, this->mspIdxConv);

      // Loop over group cpus
      this->mBTypes.reserve(FrameworkMacro::transformCpus(this->mTraId).size());
      this->mFTypes.reserve(FrameworkMacro::transformCpus(this->mTraId).size());
      for(int id = 0; id < FrameworkMacro::transformCpus(this->mTraId).size(); id++)
      {
      	 // Synchronize 
     	   FrameworkMacro::syncTransform(this->mTraId);

         // Create TBwdB datatypes
         #if defined GEOMHDISCC_MPIPACK_MANUAL
            this->mBTypes.push_back(std::vector<typename MpiConverterTools<TBwdB::FieldDimension>::Coordinate>());
            MpiConverterTools<TBwdB::FieldDimension>::buildBwdDatatype(this->mBTypes.back(), localBwdMap, spRes, this->mTraId, bTmp, id);
         #else
            MPI_Datatype type = MpiConverterTools<TBwdB::FieldDimension>::buildBwdDatatype(localBwdMap, spRes, this->mTraId, bTmp, id);
            this->mBTypes.push_back(type);
         #endif //defined GEOMHDISCC_MPIPACK_MANUAL

      	 // Synchronize 
     	   FrameworkMacro::syncTransform(this->mTraId);

         // Create TFwdA datatypes
         #if defined GEOMHDISCC_MPIPACK_MANUAL
            this->mFTypes.push_back(std::vector<typename MpiConverterTools<TFwdA::FieldDimension>::Coordinate>());
            MpiConverterTools<TFwdA::FieldDimension>::buildFwdDatatype(this->mFTypes.back(), localFwdMap, spRes, this->mTraId, fTmp, id);
         #else
            type = MpiConverterTools<TFwdA::FieldDimension>::buildFwdDatatype(localFwdMap, spRes, this->mTraId, fTmp, id);
            this->mFTypes.push_back(type);
         #endif //defined GEOMHDISCC_MPIPACK_MANUAL
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::syncFwdBuffer()
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

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::syncBwdBuffer()
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

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::syncFwdRecvBuffer()
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

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::syncBwdRecvBuffer()
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

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::syncFwdSendBuffer()
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

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::syncBwdSendBuffer()
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

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::sendFwd(const TFwdA& data)
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
         #if defined GEOMHDISCC_MPIPACK_MANUAL
            MpiConverterTools<TFwdA::FieldDimension>::template pack<typename TFwdA::PointType>(this->mspFBuffers->at(id), this->mspFBuffers->pos(id), data, this->mFTypes.at(id));
         #else
            int ierr = MPI_Pack(const_cast<typename TFwdA::PointType *>(data.data().data()), 1, this->mFTypes.at(id), this->mspFBuffers->at(id), this->sizeFPacket(id), &(this->mspFBuffers->pos(id)), FrameworkMacro::transformComm(this->mTraId));
            FrameworkMacro::check(ierr, 761);
         #endif //defined GEOMHDISCC_MPIPACK_MANUAL
      }

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::FWDSENDCONV);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::sendBwd(const TBwdB& data)
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
         #if defined GEOMHDISCC_MPIPACK_MANUAL
            MpiConverterTools<TBwdB::FieldDimension>::template pack<typename TBwdB::PointType>(this->mspBBuffers->at(id), this->mspBBuffers->pos(id), data, this->mBTypes.at(id));
         #else
            int ierr = MPI_Pack(const_cast<typename TBwdB::PointType *>(data.data().data()), 1, this->mBTypes.at(id), this->mspBBuffers->at(id), this->sizeBPacket(id), &(this->mspBBuffers->pos(id)), FrameworkMacro::transformComm(this->mTraId));
            FrameworkMacro::check(ierr, 762);
         #endif //defined GEOMHDISCC_MPIPACK_MANUAL
      }

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWDSENDCONV);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::receiveFwd(TFwdA &rData)
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
            #ifdef GEOMHDISCC_DEBUG
               MPI_Status stats[this->nFCpu()];
               int ierr = MPI_Waitsome(this->nFCpu(), this->pRecvFRequests(this->mPacks), &count, idx.data(), stats);
               FrameworkMacro::check(ierr, 771);
               DebuggerMacro_msg("Received FWD packs", 5);
            #else
               int ierr = MPI_Waitsome(this->nFCpu(), this->pRecvFRequests(this->mPacks), &count, idx.data(), MPI_STATUSES_IGNORE);
               FrameworkMacro::check(ierr, 771);
            #endif //GEOMHDISCC_DEBUG

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

               #if defined GEOMHDISCC_MPIPACK_MANUAL
                  MpiConverterTools<TFwdA::FieldDimension>::template unpack<typename TFwdA::PointType>(rData, this->mFTypes.at(pos), this->mspFBuffers->at(pos), this->mspFBuffers->pos(pos));
               #else
                  int ierr = MPI_Unpack(this->mspFBuffers->at(pos), this->sizeFPacket(pos), &(this->mspFBuffers->pos(pos)), rData.rData().data(), 1, this->mFTypes.at(pos), FrameworkMacro::transformComm(this->mTraId));
                  FrameworkMacro::check(ierr, 763);
               #endif //defined GEOMHDISCC_MPIPACK_MANUAL
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

            #if defined GEOMHDISCC_MPIPACK_MANUAL
               MpiConverterTools<TFwdA::FieldDimension>::template unpack<typename TFwdA::PointType>(rData, this->mFTypes.at(id), this->mspFBuffers->at(id), this->mspFBuffers->pos(id));
            #else
               int ierr = MPI_Unpack(this->mspFBuffers->at(id), this->sizeFPacket(id), &(this->mspFBuffers->pos(id)), rData.rData().data(), 1, this->mFTypes.at(id), FrameworkMacro::transformComm(this->mTraId));
               FrameworkMacro::check(ierr, 764);
            #endif //defined GEOMHDISCC_MPIPACK_MANUAL

         }

         // Stop detailed profiler
         DetailedProfilerMacro_stop(ProfilerMacro::FWDRECVCONV);
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::receiveBwd(TBwdB &rData)
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
            #ifdef GEOMHDISCC_DEBUG
               MPI_Status stats[this->nBCpu()];
               int ierr = MPI_Waitsome(this->nBCpu(), this->pRecvBRequests(this->mPacks), &count, idx.data(), stats);
               FrameworkMacro::check(ierr, 772);
               DebuggerMacro_msg("Received BWD packs", 5);
            #else 
               int ierr = MPI_Waitsome(this->nBCpu(), this->pRecvBRequests(this->mPacks), &count, idx.data(), MPI_STATUSES_IGNORE);
               FrameworkMacro::check(ierr, 772);
            #endif //GEOMHDISCC_DEBUG

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

               #if defined GEOMHDISCC_MPIPACK_MANUAL
                  MpiConverterTools<TBwdB::FieldDimension>::template unpack<typename TBwdB::PointType>(rData, this->mBTypes.at(pos), this->mspBBuffers->at(pos), this->mspBBuffers->pos(pos));
               #else
                  int ierr = MPI_Unpack(this->mspBBuffers->at(pos), this->sizeBPacket(pos), &(this->mspBBuffers->pos(pos)), rData.rData().data(), 1, this->mBTypes.at(pos), FrameworkMacro::transformComm(this->mTraId));
                  FrameworkMacro::check(ierr, 765);
               #endif //defined GEOMHDISCC_MPIPACK_MANUAL
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

            #if defined GEOMHDISCC_MPIPACK_MANUAL
               MpiConverterTools<TBwdB::FieldDimension>::template unpack<typename TBwdB::PointType>(rData, this->mBTypes.at(id), this->mspBBuffers->at(id), this->mspBBuffers->pos(id));
            #else
               int ierr = MPI_Unpack(this->mspBBuffers->at(id), this->sizeBPacket(id), &(this->mspBBuffers->pos(id)), rData.rData().data(), 1, this->mBTypes.at(id), FrameworkMacro::transformComm(this->mTraId));
               FrameworkMacro::check(ierr, 766);
            #endif //defined GEOMHDISCC_MPIPACK_MANUAL
         }

         // Stop detailed profiler
         DetailedProfilerMacro_stop(ProfilerMacro::BWDRECVCONV);
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::initLists()
   {
      int sze;
      std::vector<int> unusedF;
      std::vector<int> unusedB;

      // Loop over all CPUs
      for(int i = 0; i < FrameworkMacro::transformCpus(this->mTraId).size(); i++)
      {
         // Compute buffer sizes for F group
         #if defined GEOMHDISCC_MPIPACK_MANUAL
            sze = this->mFTypes.at(i).size();
         #else
            MPI_Pack_size(1, this->mFTypes.at(i), FrameworkMacro::transformComm(this->mTraId), &sze);
         #endif //defined GEOMHDISCC_MPIPACK_MANUAL
         if(sze != 0)
         {
            this->mFSizes.push_back(sze);
            this->mFCpuGroup.push_back(i);

         // Get a list of unused forward datatypes
         } else
         {
            unusedF.push_back(i);
         }

         // Compute buffer sizes for B group
         #if defined GEOMHDISCC_MPIPACK_MANUAL
            sze = this->mBTypes.at(i).size();
         #else
            MPI_Pack_size(1, this->mBTypes.at(i), FrameworkMacro::transformComm(this->mTraId), &sze);
         #endif //defined GEOMHDISCC_MPIPACK_MANUAL
         if(sze != 0)
         {
            this->mBSizes.push_back(sze);
            this->mBCpuGroup.push_back(i);

         // Get a list of unused backward datatypes
         } else
         {
            unusedB.push_back(i);
         }
      }

      std::vector<int>::reverse_iterator rit;
      // Erase the unused forward datatypes
      for(rit = unusedF.rbegin(); rit != unusedF.rend(); ++rit)
      {
         this->mFTypes.erase(this->mFTypes.begin() + *rit);
      }

      // Erase the unused backward datatypes
      for(rit = unusedB.rbegin(); rit != unusedB.rend(); ++rit)
      {
         this->mBTypes.erase(this->mBTypes.begin() + *rit);
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::cleanupTypes()
   {
      #if defined GEOMHDISCC_MPIPACK_MANUAL
         // No cleanup is required
         
      #else
         // Cleanup the F Types
         typename std::vector<MPI_Datatype>::iterator  it;

         for(it = this->mFTypes.begin(); it != this->mFTypes.end(); ++it)
         {
            MPI_Type_free(&(*it));
         }

         for(it = this->mBTypes.begin(); it != this->mBTypes.end(); ++it)
         {
            MPI_Type_free(&(*it));
         }
      #endif //defined GEOMHDISCC_MPIPACK_MANUAL
   }

#ifdef GEOMHDISCC_STORAGEPROFILE
   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::profileStorage() const
   {
      MHDFloat memTypes = 0.0;
      MHDFloat memComm = 0.0;

      // General communication storage
      memComm += 4.0*5.0;
      memComm += 4.0*(this->mForwardPacks.size() + this->mBackwardPacks.size());
      memComm += 4.0*(this->mFSizes.size() + this->nFCpu());
      memComm += 4.0*(this->mBSizes.size() + this->nBCpu());

      // Requests communication storage
      memComm += 0.0;

      // MPI datatypes storage
      memTypes += 0.0;

      StorageProfilerMacro_update(StorageProfilerMacro::MPI, memTypes + memComm);
      #ifdef GEOMHDISCC_STORAGEPROFILER_DETAILED
         StorageProfilerMacro_update(StorageProfilerMacro::MPITYPES, memTypes);

         StorageProfilerMacro_update(StorageProfilerMacro::MPICOMM, memComm);
      #endif // GEOMHDISCC_STORAGEPROFILER_DETAILED
   }
#endif // GEOMHDISCC_STORAGEPROFILE

}
}

#endif // MPICONVERTER_HPP
