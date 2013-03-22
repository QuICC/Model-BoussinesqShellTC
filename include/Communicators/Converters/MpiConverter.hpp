/** \file MpiConverter.hpp
 *  \brief Implementation of the FDSH MPI data converter.
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
    * \brief Implementation of the FDSH MPI data converter.
    */
   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> class MpiConverter: public MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>
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
          * @brief Start communication for forward transform
          */
         virtual void initiateForwardCommunication();

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

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB>::MpiConverter()
   {
      // Check that all dimensions match
      Debug::StaticAssert< (TFwdA::FieldDimension == TBwdA::FieldDimension) >();
      Debug::StaticAssert< (TBwdA::FieldDimension == TFwdB::FieldDimension) >();
      Debug::StaticAssert< (TFwdB::FieldDimension == TBwdB::FieldDimension) >();

      // Check that the data type is the same
      Debug::StaticTypeAssert<typename TFwdA::PointType , typename TBwdB::PointType>();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB>::~MpiConverter()
   {
      // Cleanup the types memory
      this->cleanupTypes();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB>::convertFwd(const TFwdA &in, StoragePairProviderMacro<TFwdB, TBwdB>  &storage)
   {
      DetailedProfilerMacro_start(ProfilerMacro::FWDCONVSEND);

      // Send the data
      this->sendFwd(in);

      DetailedProfilerMacro_stop(ProfilerMacro::FWDCONVSEND);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB>::convertBwd(const TBwdB &in, StoragePairProviderMacro<TFwdA, TBwdA>  &storage)
   {
      DetailedProfilerMacro_start(ProfilerMacro::BWDCONVSEND);

      // Send the data
      this->sendBwd(in);

      DetailedProfilerMacro_stop(ProfilerMacro::BWDCONVSEND);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> TFwdA& MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB>::getFwd(StoragePairProviderMacro<TFwdA, TBwdA>  &storage)
   {
      DetailedProfilerMacro_start(ProfilerMacro::BWDCONVRECV);

      // Get storage for output value 
      TFwdA &rOut = storage.provideFwd();

      // Receive converted data
      this->receiveFwd(rOut);

      DetailedProfilerMacro_stop(ProfilerMacro::BWDCONVRECV);

      return rOut;
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> TBwdB& MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB>::getBwd(StoragePairProviderMacro<TFwdB, TBwdB>  &storage)
   {
      DetailedProfilerMacro_start(ProfilerMacro::FWDCONVRECV);

      // Get storage for output value 
      TBwdB &rOut = storage.provideBwd();

      // Receive converted data
      this->receiveBwd(rOut);

      DetailedProfilerMacro_stop(ProfilerMacro::FWDCONVRECV);

      return rOut;
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB>::setup()
   {
      // initialise the send and receive positions
      this->initPositions();

      // setup the communication requests
      this->setupRequests();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB>::setupCommunication(const int packs)
   {
      // Store the number of packs in the next communication
      this->mPacks = packs;
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB>::initiateForwardCommunication()
   {
      // Synchronize CPUs
      FrameworkMacro::synchronize();

      // Store the number of packs in active send
      this->mActiveBSendPacks = this->mPacks;

      // Prepost the receive calls
      MPI_Startall(this->nFCpu(), this->pRecvFRequests(this->mPacks));
      this->resetRecvPositions();
      this->mIsReceiving = true;

      // Post non blocking send calls 
      MPI_Startall(this->nBCpu(), this->pSendBRequests(this->mPacks));
      this->resetSendPositions();
      this->mIsSending = true;
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB>::initiateBackwardCommunication()
   {
      // Synchronize CPUs
      FrameworkMacro::synchronize();

      // Store the number of packs in active send
      this->mActiveFSendPacks = this->mPacks;

      // Prepost the receive calls
      MPI_Startall(this->nBCpu(), this->pRecvBRequests(this->mPacks));
      this->resetRecvPositions();
      this->mIsReceiving = true;

      // Post non blocking send calls 
      MPI_Startall(this->nFCpu(), this->pSendFRequests(this->mPacks));
      this->resetSendPositions();
      this->mIsSending = true;
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB>::init(SharedResolution spRes, const Dimensions::Transform::Id fwdDim, TFwdA &fwdTmp, TBwdB &bwdTmp, const ArrayI& fwdPacks, const ArrayI& bwdPacks)
   {
      // Store the possible pack sizes
      this->mForwardPacks = fwdPacks;
      this->mBackwardPacks = bwdPacks;

      // initialise the previous active packs to a possible value
      this->mActiveFSendPacks = this->mBackwardPacks(0);
      this->mActiveBSendPacks = this->mForwardPacks(0);

      // initialise the data types
      this->initTypes(spRes, fwdDim, fwdTmp, bwdTmp);

      // initialise the size and CPU lists
      this->initLists();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB>::initTypes(SharedResolution spRes, const Dimensions::Transform::Id fwdDim, TFwdA &fTmp, TBwdB &bTmp)
   {
      // Loop over all the cpus
      for(int id = 0; id < FrameworkMacro::nCpu(); id++)
      {
         // Create TBwdB datatypes
         MPI_Datatype type = MpiConverterTools<TBwdB::FieldDimension>::buildBwdDatatype(spRes, fwdDim, bTmp, id);
         this->mBTypes.push_back(type);

         // Create TFwdA datatypes
         type = MpiConverterTools<TFwdA::FieldDimension>::buildFwdDatatype(spRes, fwdDim, fTmp, id);
         this->mFTypes.push_back(type);
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB>::sendFwd(const TFwdA& data)
   {
      // Check if communication interface is busy
      if(this->mIsSending)
      {
         int flag;
         // Make sure previous communication has finished
         // (depending on MPI implementation the double test (test+wait) is required for the expected result)
         MPI_Testall(this->nFCpu(), this->pSendBRequests(this->mActiveBSendPacks), &flag, MPI_STATUSES_IGNORE);

         // If not all are ready yet wait for completion
         if(! flag)
         {
            MPI_Waitall(this->nFCpu(), this->pSendBRequests(this->mActiveBSendPacks), MPI_STATUSES_IGNORE);
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

      // Pack data into send buffer
      for(int id = 0; id < this->nFCpu(); ++id)
      {
         MPI_Pack(const_cast<typename TFwdA::PointType *>(data.data().data()), 1, this->mFTypes.at(id), this->mspFBuffers->at(id), this->sizeFPacket(id), &(this->mSendPositions.at(id)), MPI_COMM_WORLD);
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB>::sendBwd(const TBwdB& data)
   {
      // Check if communication interface is busy
      if(this->mIsSending)
      {
         int flag;
         // Make sure previous communication has finished
         // (depending on MPI implementation the double test (test+wait) is required for the expected result)
         MPI_Testall(this->nBCpu(), this->pSendFRequests(this->mActiveFSendPacks), &flag, MPI_STATUSES_IGNORE);
         // If not all are ready yet wait for completion
         if(! flag)
         {
            MPI_Waitall(this->nBCpu(), this->pSendFRequests(this->mActiveFSendPacks), MPI_STATUSES_IGNORE);
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

      // Pack data into send buffer
      for(int id = 0; id < this->nBCpu(); ++id)
      {
         MPI_Pack(const_cast<typename TBwdB::PointType *>(data.data().data()), 1, this->mBTypes.at(id), this->mspBBuffers->at(id), this->sizeBPacket(id), &(this->mSendPositions.at(id)), MPI_COMM_WORLD);
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB>::receiveFwd(TFwdA &rData)
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
            // Wait for some of the requests to finish
            MPI_Waitsome(this->nFCpu(), this->pRecvFRequests(this->mPacks), &count, idx.data(), MPI_STATUSES_IGNORE);

            // Unpack already received data from receive buffer
            for(int id = 0; id < count; ++id)
            {
               MPI_Unpack(this->mspFBuffers->at(idx(id)), this->sizeFPacket(idx(id)), &(this->mRecvPositions.at(idx(id))), rData.rData().data(), 1, this->mFTypes.at(idx(id)), MPI_COMM_WORLD);
            }

            // Update the number of missing receives
            keepWaiting -= count;
         }

         // Reset communication status
         this->mIsReceiving = false;
      } else
      {
         // Unpack data from receive buffer
         for(int id = 0; id < this->nFCpu(); ++id)
         {
            MPI_Unpack(this->mspFBuffers->at(id), this->sizeFPacket(id), &(this->mRecvPositions.at(id)), rData.rData().data(), 1, this->mFTypes.at(id), MPI_COMM_WORLD);
         }
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB>::receiveBwd(TBwdB &rData)
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
            // Wait for some of the requests to finish
            MPI_Waitsome(this->nBCpu(), this->pRecvBRequests(this->mPacks), &count, idx.data(), MPI_STATUSES_IGNORE);

            // Unpack already received data from receive buffer
            for(int id = 0; id < count; ++id)
            {
               MPI_Unpack(this->mspBBuffers->at(idx(id)), this->sizeBPacket(idx(id)), &(this->mRecvPositions.at(idx(id))), rData.rData().data(), 1, this->mBTypes.at(idx(id)), MPI_COMM_WORLD);
            }

            // Update the number of missing receives
            keepWaiting -= count;
         }

         // Reset communication status
         this->mIsReceiving = false;
      } else
      {
         // Unpack data from receive buffer
         for(int id = 0; id < this->nBCpu(); ++id)
         {
            MPI_Unpack(this->mspBBuffers->at(id), this->sizeBPacket(id), &(this->mRecvPositions.at(id)), rData.rData().data(), 1, this->mBTypes.at(id), MPI_COMM_WORLD);
         }
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB>::initLists()
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

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB>::cleanupTypes()
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
   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB>::profileStorage() const
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
