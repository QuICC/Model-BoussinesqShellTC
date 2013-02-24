/** \file SerialConverter.hpp
 *  \brief Implementation of the serial data converter.
 */

#ifndef SERIALCONVERTER_HPP
#define SERIALCONVERTER_HPP

// Configuration includes
//
#include "Profiler/ProfilerMacro.h"
#include "StaticAssert/StaticAssert.hpp"

// System includes
//
#include <assert.h>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/Enums/Dimensions.hpp"
#include "StorageProviders/StoragePairProvider.hpp"
#include "Communicators/Converters/SerialConverterBase.hpp"

namespace GeoMHDiSCC {

   /**
    * \brief Implementation of the serial data converter.
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
          * @brief Get the converted data from TBwdA to TFwdB conversion
          */
         virtual TFwdA& getFwd(StoragePairProvider<TFwdA, TBwdA>& storage);

         /**
          * @brief Get the converted data from TFwdB to TBwdA conversion
          */
         virtual TBwdB& getBwd(StoragePairProvider<TFwdB, TBwdB>& storage);

         /**
          * @brief Setup the converter
          */
         virtual void setup();

         /**
          * @brief Setup upcoming communication
          *
          * @param packs Number of packets in communication packing
          */
         void setupCommunication(const int packs);

         /**
          * @brief Start communication for forward transform
          */
         void initiateForwardCommunication();

         /**
          * @brief Start communication for backward transform
          */
         void initiateBackwardCommunication();

         /**
          * @brief Initialise the converter
          *
          * @param spRes      Shared Resolution
          * @param fwdDim     Dimension index for forward transform
          */
         void init(SharedResolution spRes, const int fwdDim);

      #ifdef GEOMHDISCC_STORAGEPROFILE
         /**
         * @brief Do storage profiling
         */
         virtual void profileStorage() const;
      #endif // GEOMHDISCC_STORAGEPROFILE
         
      protected:

      private:
   };

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> SerialConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::SerialConverter()
   {
      // Check that all dimensions match
      StaticAssert< (TFwdA::FieldDimension == TBwdA::FieldDimension) >();
      StaticAssert< (TBwdA::FieldDimension == TFwdB::FieldDimension) >();
      StaticAssert< (TFwdB::FieldDimension == TBwdB::FieldDimension) >();

      // Check that the data type is the same
      StaticTypeAssert<typename TFwdA::CoefficientType , typename TBwdB::CoefficientType>();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> SerialConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::~SerialConverter()
   {
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void SerialConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::setup()
   {
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void SerialConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::setupCommunication(const int packs)
   {
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void SerialConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::initiateForwardCommunication()
   {
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void SerialConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::initiateBackwardCommunication()
   {
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> TFwdA& SerialConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::getFwd(StoragePairProvider<TFwdA, TBwdA>  &storage)
   {
      // Recover storage from provider
      TFwdA &rOut = storage.recoverFwd();

      return rOut;
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> TBwdB& SerialConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::getBwd(StoragePairProvider<TFwdB, TBwdB>  &storage)
   {
      // Recover storage from provider
      TBwdB &rOut = storage.recoverBwd();

      return rOut;
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void SerialConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::init(SharedResolution spRes, const int fwdDim)
   {
      // Store the shared pointer to the transform resolution
      this->mspTRes = spRes->cpu()->dim(fwdDim);
   }

#ifdef GEOMHDISCC_STORAGEPROFILE
   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void SerialConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::profileStorage() const
   {
   }
#endif // GEOMHDISCC_STORAGEPROFILE

}

#endif // SERIALCONVERTER_HPP
