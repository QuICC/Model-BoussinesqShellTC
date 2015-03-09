/**
 * @file MpiConverter.hpp
 * @brief Implementation of the FDSH MPI data converter 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef MPICONVERTER_HPP
#define MPICONVERTER_HPP

// Configuration includes
//
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
         virtual void setupCommunication(const int packs);

         /**
          * @brief Start persistent send for forward transform
          */
         virtual void initiateForwardSend();

         /**
          * @brief Post persistent receive for forward transform
          */
         virtual void prepareForwardReceive();

         /**
          * @brief Start communication for forward transform
          */
         virtual void initiateForwardCommunication();

         /**
          * @brief Start persistent send for backward transform
          */
         virtual void initiateBackwardSend();

         /**
          * @brief Post persistent receive for backward transform
          */
         virtual void prepareBackwardReceive();

         /**
          * @brief Start communication for backward transform
          */
         virtual void initiateBackwardCommunication();

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
          * @param fwdDim  Dimension index for forward transform
          * @param fTmp    TFwdA temporary
          * @param bTmp    TBwdB temporary
          */
         void initTypes(SharedResolution spRes, const Dimensions::Transform::Id fwdDim, TFwdA& fTmp, TBwdB& bTmp);

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
          * @brief Cleanup the data types
          */
         void cleanupTypes();

         /**
          * @brief Storage for the forward datatypes
          */
         std::vector<MPI_Datatype>  mFTypes;

         /**
          * @brief Storage for the backward datatypes
          */
         std::vector<MPI_Datatype> mBTypes;
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
      // initialise the send and receive positions
      this->initPositions();

      // setup the communication requests
      this->setupRequests();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::setupCommunication(const int packs)
   {
      // Store the number of packs in the next communication
      this->mPacks = packs;
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::prepareForwardReceive()
   {
      // Don't do anything if the number of packs is zero
      if(this->mPacks > 0)
      {
         // Make sure calls are posted at the right moment
         int flag;
         MPI_Testall(this->nFCpu(), this->pRecvFRequests(this->mPacks), &flag, MPI_STATUSES_IGNORE);
         if(!flag)
         {
            MPI_Abort(MPI_COMM_WORLD, 99);
         }

         // Prepost the receive calls
         MPI_Startall(this->nFCpu(), this->pRecvFRequests(this->mPacks));
         this->resetRecvPositions();
         this->mIsReceiving = true;
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::initiateBackwardSend()
   {
      // Don't do anything if the number of packs is zero
      if(this->mPacks > 0)
      {
         // Store the number of packs in active send
         this->mActiveBSendPacks = this->mPacks;

         // Make sure calls are posted at the right moment
         int flag;
         MPI_Testall(this->nBCpu(), this->pSendBRequests(this->mPacks), &flag, MPI_STATUSES_IGNORE);
         if(!flag)
         {
            MPI_Abort(MPI_COMM_WORLD, 99);
         }

         // Post non blocking send calls 
         MPI_Startall(this->nBCpu(), this->pSendBRequests(this->mPacks));
         this->resetSendPositions();
         this->mIsSending = true;
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::initiateForwardCommunication()
   {
      // Prepose forward receive
      this->prepareForwardReceive();

      // Start backward send
      this->initiateBackwardSend();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::prepareBackwardReceive()
   {
      // Don't do anything if the number of packs is zero
      if(this->mPacks > 0)
      {
         // Make sure calls are posted at the right moment
         int flag;
         MPI_Testall(this->nBCpu(), this->pRecvBRequests(this->mPacks), &flag, MPI_STATUSES_IGNORE);
         if(!flag)
         {
            MPI_Abort(MPI_COMM_WORLD, 99);
         }

         // Prepost the receive calls
         MPI_Startall(this->nBCpu(), this->pRecvBRequests(this->mPacks));
         this->resetRecvPositions();
         this->mIsReceiving = true;
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::initiateForwardSend()
   {
      // Don't do anything if the number of packs is zero
      if(this->mPacks > 0)
      {
         // Store the number of packs in active send
         this->mActiveFSendPacks = this->mPacks;

         // Make sure calls are posted at the right moment
         int flag;
         MPI_Testall(this->nFCpu(), this->pSendFRequests(this->mPacks), &flag, MPI_STATUSES_IGNORE);
         if(!flag)
         {
            MPI_Abort(MPI_COMM_WORLD, 99);
         }

         // Post non blocking send calls 
         MPI_Startall(this->nFCpu(), this->pSendFRequests(this->mPacks));
         this->resetSendPositions();
         this->mIsSending = true;
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::initiateBackwardCommunication()
   {
      // Prepose backward receive
      this->prepareBackwardReceive();

      // Start forward send
      this->initiateForwardSend();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::init(SharedResolution spRes, const Dimensions::Transform::Id fwdDim, TFwdA &fwdTmp, TBwdB &bwdTmp, const ArrayI& fwdPacks, const ArrayI& bwdPacks)
   {
      // Store the possible pack sizes
      this->mForwardPacks = fwdPacks;
      this->mBackwardPacks = bwdPacks;

      // initialise the previous active packs to a possible value
      if(this->mForwardPacks.size() > 0 && this->mBackwardPacks.size() > 0)
      {
         this->mActiveBSendPacks = this->mForwardPacks(0);
      } else
      {
         this->mActiveBSendPacks = 0;
      }
      if(this->mBackwardPacks.size() > 0 && this->mForwardPacks.size() > 0)
      {
         this->mActiveFSendPacks = this->mBackwardPacks(0);
      } else
      {
         this->mActiveFSendPacks = 0;
      }
      
      // Create index converter
      this->mspIdxConv = SharedPtrMacro<TIdx>(new TIdx(spRes, fwdDim));

      // initialise the data types
      this->initTypes(spRes, fwdDim, fwdTmp, bwdTmp);

      // initialise the size and CPU lists
      this->initLists();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::initTypes(SharedResolution spRes, const Dimensions::Transform::Id fwdDim, TFwdA &fTmp, TBwdB &bTmp)
   {
      // Loop over all the cpus
      for(int id = 0; id < FrameworkMacro::nCpu(); id++)
      {
         // Create TBwdB datatypes
         MPI_Datatype type = MpiConverterTools<TBwdB::FieldDimension>::buildBwdDatatype(spRes, fwdDim, bTmp, id, this->mspIdxConv);
         this->mBTypes.push_back(type);

         // Create TFwdA datatypes
         type = MpiConverterTools<TFwdA::FieldDimension>::buildFwdDatatype(spRes, fwdDim, fTmp, id, this->mspIdxConv);
         this->mFTypes.push_back(type);
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::sendFwd(const TFwdA& data)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::FWDSENDWAIT);

      // Check if communication interface is busy
      if(this->mIsSending)
      {
         int flag;
         // Make sure previous communication has finished
         if(this->mActiveBSendPacks > 0)
         {
            // (depending on MPI implementation the double test (test+wait) is required for the expected result)
            MPI_Testall(this->nFCpu(), this->pSendBRequests(this->mActiveBSendPacks), &flag, MPI_STATUSES_IGNORE);

            // If not all are ready yet wait for completion
            if(! flag)
            {
               MPI_Waitall(this->nFCpu(), this->pSendBRequests(this->mActiveBSendPacks), MPI_STATUSES_IGNORE);
            }
         }

         // Make sure communication to be used are all finished
         // (depending on MPI implementation the double test (test+wait) is required for the expected result)
         MPI_Testall(this->nFCpu(), this->pSendFRequests(this->mPacks), &flag, MPI_STATUSES_IGNORE);

         // If not all are ready yet wait for completion
         if(! flag)
         {
            MPI_Waitall(this->nFCpu(), this->pSendFRequests(this->mPacks), MPI_STATUSES_IGNORE);
         }

         // Reset communication status
         this->mIsSending = false;
      }

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::FWDSENDWAIT);

      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::FWDSENDCONV);

      // Pack data into send buffer
      for(int id = 0; id < this->nFCpu(); ++id)
      {
         MPI_Pack(const_cast<typename TFwdA::PointType *>(data.data().data()), 1, this->mFTypes.at(id), this->mspFBuffers->at(id), this->sizeFPacket(id), &(this->mSendPositions.at(id)), MPI_COMM_WORLD);
      }

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::FWDSENDCONV);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::sendBwd(const TBwdB& data)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWDSENDWAIT);

      // Check if communication interface is busy
      if(this->mIsSending)
      {
         int flag;
         // Make sure previous communication has finished
         if(this->mActiveFSendPacks > 0)
         {
            // (depending on MPI implementation the double test (test+wait) is required for the expected result)
            MPI_Testall(this->nBCpu(), this->pSendFRequests(this->mActiveFSendPacks), &flag, MPI_STATUSES_IGNORE);
            // If not all are ready yet wait for completion
            if(! flag)
            {
               MPI_Waitall(this->nBCpu(), this->pSendFRequests(this->mActiveFSendPacks), MPI_STATUSES_IGNORE);
            }
         }

         // Make sure communication to be used are all finished
         // (depending on MPI implementation the double test (test+wait) is required for the expected result)
         MPI_Testall(this->nBCpu(), this->pSendBRequests(this->mPacks), &flag, MPI_STATUSES_IGNORE);
         // If not all are ready yet wait for completion
         if(! flag)
         {
            MPI_Waitall(this->nBCpu(), this->pSendBRequests(this->mPacks), MPI_STATUSES_IGNORE);
         }

         // Reset communication status
         this->mIsSending = false;
      }

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWDSENDWAIT);

      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWDSENDCONV);

      // Pack data into send buffer
      for(int id = 0; id < this->nBCpu(); ++id)
      {
         MPI_Pack(const_cast<typename TBwdB::PointType *>(data.data().data()), 1, this->mBTypes.at(id), this->mspBBuffers->at(id), this->sizeBPacket(id), &(this->mSendPositions.at(id)), MPI_COMM_WORLD);
      }

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWDSENDCONV);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::receiveFwd(TFwdA &rData)
   {
      // Check if communication interface is busy
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
            MPI_Waitsome(this->nFCpu(), this->pRecvFRequests(this->mPacks), &count, idx.data(), MPI_STATUSES_IGNORE);

            // Stop detailed profiler
            DetailedProfilerMacro_stop(ProfilerMacro::FWDRECVWAIT);

            // Start detailed profiler
            DetailedProfilerMacro_start(ProfilerMacro::FWDRECVCONV);

            // Unpack already received data from receive buffer
            for(int id = 0; id < count; ++id)
            {
               MPI_Unpack(this->mspFBuffers->at(idx(id)), this->sizeFPacket(idx(id)), &(this->mRecvPositions.at(idx(id))), rData.rData().data(), 1, this->mFTypes.at(idx(id)), MPI_COMM_WORLD);
            }

            // Stop detailed profiler
            DetailedProfilerMacro_stop(ProfilerMacro::FWDRECVCONV);

            // Update the number of missing receives
            keepWaiting -= count;
         }

         // Reset communication status
         this->mIsReceiving = false;
      } else
      {
         // Start detailed profiler
         DetailedProfilerMacro_start(ProfilerMacro::FWDRECVCONV);

         // Unpack data from receive buffer
         for(int id = 0; id < this->nFCpu(); ++id)
         {
            MPI_Unpack(this->mspFBuffers->at(id), this->sizeFPacket(id), &(this->mRecvPositions.at(id)), rData.rData().data(), 1, this->mFTypes.at(id), MPI_COMM_WORLD);
         }

         // Stop detailed profiler
         DetailedProfilerMacro_stop(ProfilerMacro::FWDRECVCONV);
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::receiveBwd(TBwdB &rData)
   {
      // Check if communication interface is busy
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
            MPI_Waitsome(this->nBCpu(), this->pRecvBRequests(this->mPacks), &count, idx.data(), MPI_STATUSES_IGNORE);

            // Stop detailed profiler
            DetailedProfilerMacro_stop(ProfilerMacro::BWDRECVWAIT);

            // Start detailed profiler
            DetailedProfilerMacro_start(ProfilerMacro::BWDRECVCONV);

            // Unpack already received data from receive buffer
            for(int id = 0; id < count; ++id)
            {
               MPI_Unpack(this->mspBBuffers->at(idx(id)), this->sizeBPacket(idx(id)), &(this->mRecvPositions.at(idx(id))), rData.rData().data(), 1, this->mBTypes.at(idx(id)), MPI_COMM_WORLD);
            }

            // Stop detailed profiler
            DetailedProfilerMacro_stop(ProfilerMacro::BWDRECVCONV);

            // Update the number of missing receives
            keepWaiting -= count;
         }

         // Reset communication status
         this->mIsReceiving = false;
      } else
      {
         // Start detailed profiler
         DetailedProfilerMacro_start(ProfilerMacro::BWDRECVCONV);

         // Unpack data from receive buffer
         for(int id = 0; id < this->nBCpu(); ++id)
         {
            MPI_Unpack(this->mspBBuffers->at(id), this->sizeBPacket(id), &(this->mRecvPositions.at(id)), rData.rData().data(), 1, this->mBTypes.at(id), MPI_COMM_WORLD);
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
      for(int id = 0; id < FrameworkMacro::nCpu(); ++id)
      {
         // Compute buffer sizes for F group
         MPI_Pack_size(1, this->mFTypes.at(id), MPI_COMM_WORLD, &sze);
         if(sze != 0)
         {
            this->mFSizes.push_back(sze);
            this->mFCpuGroup.push_back(id);

         // Get a list of unused forward datatypes
         } else
         {
            unusedF.push_back(id);
         }

         // Compute buffer sizes for B group
         MPI_Pack_size(1, this->mBTypes.at(id), MPI_COMM_WORLD, &sze);
         if(sze != 0)
         {
            this->mBSizes.push_back(sze);
            this->mBCpuGroup.push_back(id);

         // Get a list of unused backward datatypes
         } else
         {
            unusedB.push_back(id);
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
   }

#ifdef GEOMHDISCC_STORAGEPROFILE
   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::profileStorage() const
   {
      MHDFloat memTypes = 0.0;
      MHDFloat memComm = 0.0;

      // General communication storage
      memComm += 4.0*5.0;
      memComm += 4.0*(this->mForwardPacks.size() + this->mBackwardPacks.size());
      memComm += 4.0*(this->mFSizes.size() + this->mFCpuGroup.size());
      memComm += 4.0*(this->mBSizes.size() + this->mBCpuGroup.size());
      memComm += 4.0*(this->mSendPositions.size() + this->mRecvPositions.size());

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
