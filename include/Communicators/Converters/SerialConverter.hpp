/**
 * @file SerialConverter.hpp
 * @brief Implementation of the serial data converter 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SERIALCONVERTER_HPP
#define SERIALCONVERTER_HPP

// Debug includes
//
#include "Profiler/ProfilerMacro.h"
#include "StaticAsserts/StaticAssert.hpp"

// Configuration includes
//

// System includes
//
#include <cassert>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Enums/Dimensions.hpp"
#include "Enums/DimensionTools.hpp"
#include "Enums/TransformDirection.hpp"
#include "StorageProviders/StoragePairProviderMacro.h"
#include "Communicators/Converters/SerialConverterBase.hpp"

namespace GeoMHDiSCC {

namespace Parallel {

   /**
    * @brief Implementation of the serial data converter.
    */
   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> class SerialConverter : public SerialConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>
   {
      public:
         /**
          * @brief Constructor
          */
         SerialConverter();

         /**
          * @brief Destructor
          */
         virtual ~SerialConverter();

         /**
          * @brief Initialise the converter
          *
          * @param spRes   Shared Resolution
          * @param id      Dimension index for forward transform
          */
         void init(SharedResolution spRes, const Dimensions::Transform::Id id);

         /**
          * @brief Setup the converter
          */
         virtual void setup();

         /**
          * @brief Get the converted data from TBwdA to TFwdB conversion
          */
         virtual TFwdA& getFwd(StoragePairProviderMacro<TFwdA, TBwdA>& storage);

         /**
          * @brief Get the converted data from TFwdB to TBwdA conversion
          */
         virtual TBwdB& getBwd(StoragePairProviderMacro<TFwdB, TBwdB>& storage);

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

      private:
   };

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> SerialConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::SerialConverter()
   {
      // Check that all dimensions match
      Debug::StaticAssert< (TFwdA::FieldDimension == TBwdA::FieldDimension) >();
      Debug::StaticAssert< (TBwdA::FieldDimension == TFwdB::FieldDimension) >();
      Debug::StaticAssert< (TFwdB::FieldDimension == TBwdB::FieldDimension) >();

      // Check that the data type is the same
      Debug::StaticTypeAssert<typename TFwdA::PointType , typename TBwdB::PointType>();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> SerialConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::~SerialConverter()
   {
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> TFwdA& SerialConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::getFwd(StoragePairProviderMacro<TFwdA, TBwdA>  &storage)
   {
      // Recover storage from provider
      TFwdA &rOut = storage.recoverFwd();

      return rOut;
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> TBwdB& SerialConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::getBwd(StoragePairProviderMacro<TFwdB, TBwdB>  &storage)
   {
      // Recover storage from provider
      TBwdB &rOut = storage.recoverBwd();

      return rOut;
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void SerialConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::setup()
   {
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void SerialConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::setupCommunication(const int packs, const TransformDirection::Id direction)
   {
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void SerialConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::initiateForwardSend()
   {
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void SerialConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::prepareForwardReceive()
   {
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void SerialConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::initiateBackwardSend()
   {
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void SerialConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::prepareBackwardReceive()
   {
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void SerialConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::init(SharedResolution spRes, const Dimensions::Transform::Id id)
   {
      // Store the shared pointer to the transform resolution
      this->mspTRes = spRes->cpu()->dim(id);
      
      // Create index converter
      this->mspIdxConv = SharedPtrMacro<TIdx>(new TIdx(spRes, id));
   }

#ifdef QUICC_STORAGEPROFILE
   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void SerialConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::profileStorage() const
   {
   }
#endif // QUICC_STORAGEPROFILE

}
}

#endif // SERIALCONVERTER_HPP
