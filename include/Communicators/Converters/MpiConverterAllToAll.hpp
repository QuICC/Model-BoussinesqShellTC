/**
 * @file MpiConverterAllToAll.hpp
 * @brief Implementation of the MPI data converter with AllToAll data exchange
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef MPICONVERTERALLTOALL_HPP
#define MPICONVERTERALLTOALL_HPP

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

namespace GeoMHDiSCC {

namespace Parallel {

   /**
    * @brief Implementation of the MPI data converter with AllToAll data exchange.
    */
   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> class MpiConverterAllToAll: public MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>
   {
      public:
         /**
          * @brief Constructor
          */
         MpiConverterAllToAll();

         /**
          * @brief Destructor
          */
         virtual ~MpiConverterAllToAll();

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
          * @brief Set forward communication sizes
          */
         void setFwdCommSizes();

         /**
          * @brief Set backward communication sizes
          */
         void setBwdCommSizes();

         /**
          * @brief Setup the communication data
          */
         void setupCommData();

         /**
          * @brief Forward communication sizes
          */
         ArrayI mFCommSizes;

         /**
          * @brief backward communication sizes
          */
         ArrayI mBCommSizes;

         /**
          * @brief MPI datatypes for forward communication
          */
         std::vector<MPI_Datatype>  mFCommTypes;

         /**
          * @brief MPI datatypes for backward communication
          */
         std::vector<MPI_Datatype>  mBCommTypes;
   };

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> MpiConverterAllToAll<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::MpiConverterAllToAll()
   {
      // Check that all dimensions match
      Debug::StaticAssert< (TFwdA::FieldDimension == TBwdA::FieldDimension) >();
      Debug::StaticAssert< (TBwdA::FieldDimension == TFwdB::FieldDimension) >();
      Debug::StaticAssert< (TFwdB::FieldDimension == TBwdB::FieldDimension) >();

      // Check that the data type is the same
      Debug::StaticTypeAssert<typename TFwdA::PointType , typename TBwdB::PointType>();

      // Keep empty communications
      this->mNeedEmptyComm = true;
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> MpiConverterAllToAll<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::~MpiConverterAllToAll()
   {
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterAllToAll<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::setup()
   {
      // Initialize positions
      this->resetFwdPositions();
      this->resetBwdPositions();

      // setup the communication data
      this->setupCommData();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterAllToAll<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::setupCommunication(const int packs, const TransformDirection::Id direction)
   {
      // Store the number of packs in the next communication in direction
      this->mPacks = packs;
      this->mDirection = direction;
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterAllToAll<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::prepareForwardReceive()
   {
      // Don't do anything if the number of packs is zero
      if(this->mPacks > 0)
      {
         this->resetBwdPositions();
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterAllToAll<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::initiateBackwardSend()
   {
      // Don't do anything if the number of packs is zero
      if(this->mPacks > 0)
      {
         // Start detailed profiler
         DetailedProfilerMacro_start(ProfilerMacro::BWDRECVWAIT);

         this->resetFwdPositions();
         this->setFwdCommSizes();
         this->setBwdCommSizes();

         // All-to-all data exchange
         #if defined QUICC_MPIPACK_MANUAL
            MPI_Alltoallv(this->mspBBuffers->data(), this->mBCommSizes.data(), this->mspBBuffers->zero(), this->mBCommTypes[0], this->mspFBuffers->data(), this->mFCommSizes.data(), this->mspFBuffers->zero(), this->mFCommTypes[0], FrameworkMacro::transformComm(this->mTraId)); 
         #else
            MPI_Alltoallw(this->mspBBuffers->data(), this->mBCommSizes.data(), this->mspBBuffers->zero(), &this->mBCommTypes[0], this->mspFBuffers->data(), this->mFCommSizes.data(), this->mspFBuffers->zero(), &this->mFCommTypes[0], FrameworkMacro::transformComm(this->mTraId)); 
         #endif //defined QUICC_MPIPACK_MANUAL

         // Stop detailed profiler
         DetailedProfilerMacro_stop(ProfilerMacro::BWDRECVWAIT);
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterAllToAll<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::prepareBackwardReceive()
   {
      // Don't do anything if the number of packs is zero
      if(this->mPacks > 0)
      {
         this->resetFwdPositions();
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterAllToAll<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::initiateForwardSend()
   {
      // Don't do anything if the number of packs is zero
      if(this->mPacks > 0)
      {
         // Start detailed profiler
         DetailedProfilerMacro_start(ProfilerMacro::FWDRECVWAIT);

         this->resetBwdPositions();
         this->setFwdCommSizes();
         this->setBwdCommSizes();

         // All-to-all data exchange
         #if defined QUICC_MPIPACK_MANUAL
            MPI_Alltoallv(this->mspFBuffers->data(), this->mFCommSizes.data(), this->mspFBuffers->zero(), this->mFCommTypes[0], this->mspBBuffers->data(), this->mBCommSizes.data(), this->mspBBuffers->zero(), this->mBCommTypes[0], FrameworkMacro::transformComm(this->mTraId)); 
         #else
            MPI_Alltoallw(this->mspFBuffers->data(), this->mFCommSizes.data(), this->mspFBuffers->zero(), &this->mFCommTypes[0], this->mspBBuffers->data(), this->mBCommSizes.data(), this->mspBBuffers->zero(), &this->mBCommTypes[0], FrameworkMacro::transformComm(this->mTraId)); 
         #endif //defined QUICC_MPIPACK_MANUAL

         // Stop detailed profiler
         DetailedProfilerMacro_stop(ProfilerMacro::FWDRECVWAIT);
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterAllToAll<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::sendFwd(const TFwdA& data)
   {
      // Should be in backward direction
      assert(this->mDirection == TransformDirection::BACKWARD);

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

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterAllToAll<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::sendBwd(const TBwdB& data)
   {
      // Should be in backward direction
      assert(this->mDirection == TransformDirection::FORWARD);

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

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterAllToAll<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::receiveFwd(TFwdA &rData)
   {
      // Should be in backward direction
      assert(this->mDirection == TransformDirection::FORWARD);

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

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterAllToAll<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::receiveBwd(TBwdB &rData)
   {
      // Should be in backward direction
      assert(this->mDirection == TransformDirection::BACKWARD);

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

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterAllToAll<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::setFwdCommSizes()
   {
      // All-to-all data FWD sizes
      for(size_t i = 0; i < this->mFSizes.size(); i++)
      {
         this->mFCommSizes(i) = this->mPacks*this->mFSizes.at(i);
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterAllToAll<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::setBwdCommSizes()
   {
      // All-to-all data FWD sizes
      for(size_t i = 0; i < this->mBSizes.size(); i++)
      {
         this->mBCommSizes(i) = this->mPacks*this->mBSizes.at(i);
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterAllToAll<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::setupCommData()
   {
      #if defined QUICC_MPIPACK_MANUAL
         this->mFCommTypes.push_back(MpiTypes::template type<typename TFwdA::PointType>());

         this->mBCommTypes.push_back(MpiTypes::template type<typename TBwdB::PointType>());
      #else
         for(size_t i = 0; i < this->mFSizes.size(); i++)
         {
            this->mFCommTypes.push_back(MPI_PACKED);
         }

         for(size_t i = 0; i < this->mBSizes.size(); i++)
         {
            this->mBCommTypes.push_back(MPI_PACKED);
         }
      #endif //defined QUICC_MPIPACK_MANUAL

      this->mFCommSizes.resize(this->mFSizes.size());
      this->mBCommSizes.resize(this->mBSizes.size());
   }

#ifdef QUICC_STORAGEPROFILE
   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterAllToAll<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::profileStorage() const
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
      #ifdef QUICC_STORAGEPROFILER_DETAILED
         StorageProfilerMacro_update(StorageProfilerMacro::MPITYPES, memTypes);

         StorageProfilerMacro_update(StorageProfilerMacro::MPICOMM, memComm);
      #endif // QUICC_STORAGEPROFILER_DETAILED
   }
#endif // QUICC_STORAGEPROFILE

}
}

#endif // MPICONVERTERALLTOALL_HPP
