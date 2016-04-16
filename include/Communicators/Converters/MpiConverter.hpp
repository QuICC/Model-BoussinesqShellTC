/**
 * @file MpiConverter.hpp
 * @brief Implementation of the MPI data converter 
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
#if defined GEOMHDISCC_MPICOMM_SENDRECV
   #include "Communicators/Converters/MpiConverterSendRecv.hpp"
   #define MpiConverterCommMacro MpiConverterSendRecv
#elif defined GEOMHDISCC_MPICOMM_ALLTOALL
   #include "Communicators/Converters/MpiConverterAllToAll.hpp"
   #define MpiConverterCommMacro MpiConverterAllToAll
#endif //defined GEOMHDISCC_MPICOMM_SENDRECV
#include "Communicators/Converters/MpiConverterTools.hpp"
#include "StorageProviders/StoragePairProviderMacro.h"
#include "Resolutions/Resolution.hpp"
#include "Timers/StageTimer.hpp"

namespace GeoMHDiSCC {

namespace Parallel {

   /**
    * @brief Implementation of the MPI data converter.
    */
   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> class MpiConverter: public MpiConverterCommMacro<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>
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

      #ifdef GEOMHDISCC_STORAGEPROFILE
         /**
         * @brief Do storage profiling
         */
         virtual void profileStorage() const;
      #endif // GEOMHDISCC_STORAGEPROFILE
         
      protected:

      private:
         /**
          * @brief Initialise the datatypes for packing
          *
          * @param spRes   Shared Resolution
          * @param fTmp    TFwdA temporary
          * @param bTmp    TBwdB temporary
          */
         void initPackTypes(SharedResolution spRes, TFwdA& fTmp, TBwdB& bTmp);

         /**
          * @brief Initialise the sizes and CPU lists
          */
         void initLists();
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

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::init(SharedResolution spRes, const Dimensions::Transform::Id fwdDim, TFwdA &fwdTmp, TBwdB &bwdTmp, const ArrayI& fwdPacks, const ArrayI& bwdPacks)
   {
      StageTimer stage;
      stage.start("creating MPI packing datatypes",1);

      // Store the possible pack sizes
      this->mForwardPacks = fwdPacks;
      this->mBackwardPacks = bwdPacks;

      // Store Transform ID
      this->mTraId = fwdDim;
      
      // Create index converter
      this->mspIdxConv = SharedPtrMacro<TIdx>(new TIdx(spRes, this->mTraId));

      // initialise the data types
      this->initPackTypes(spRes, fwdTmp, bwdTmp);

      stage.done();
      stage.start("Creating CPU and size lists",1);

      // initialise the size and CPU lists
      this->initLists();

      stage.done();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::initPackTypes(SharedResolution spRes, TFwdA &fTmp, TBwdB &bTmp)
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
         if(sze != 0 || this->mNeedEmptyComm)
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
         if(sze != 0 || this->mNeedEmptyComm)
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
