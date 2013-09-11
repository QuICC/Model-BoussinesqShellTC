/**
 * @file IConverter.hpp
 * @brief Implementation of the interface for a data converter 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

#ifndef ICONVERTER_HPP
#define ICONVERTER_HPP

// Configuration includes
//

// System includes
//
#include <tr1/tuple>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "StorageProviders/StoragePairProviderMacro.h"

namespace GeoMHDiSCC {

namespace Parallel {

   /**
    * @brief Implementation of the interface for a data converter.
    */
   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> class IConverter
   {
      public:
         /**
          * @brief Constructor
          */
         IConverter();

         /**
          * @brief Destructor
          */
         virtual ~IConverter();

         /**
          * @brief Set up the converter
          */
         virtual void setup() = 0;

         /**
          * @brief Convert data from TFwdA to TBwdB
          */
         virtual void convertFwd(const TFwdA &in, StoragePairProviderMacro<TFwdB, TBwdB> &storage) = 0;

         /**
          * @brief Convert data from TBwdB to TFwdA
          */
         virtual void convertBwd(const TBwdB &in, StoragePairProviderMacro<TFwdA, TBwdA> &storage) = 0;

         /**
          * @brief Get the converted data from TBwdB to TFwdA conversion
          */
         virtual TFwdA& getFwd(StoragePairProviderMacro<TFwdA, TBwdA> &storage) = 0;

         /**
          * @brief Get the converted data from TFwdA to TBwdB conversion
          */
         virtual TBwdB& getBwd(StoragePairProviderMacro<TFwdB, TBwdB> &storage) = 0;

         /**
          * @brief Setup upcoming communication
          *
          * @param packs Number of packets in communication packing
          */
         virtual void setupCommunication(const int packs) = 0;

         /**
          * @brief Initiate communication for forward transform
          */
         virtual void initiateForwardCommunication() = 0;

         /**
          * @brief Initiate communication for backward transform
          */
         virtual void initiateBackwardCommunication() = 0;

      #ifdef GEOMHDISCC_STORAGEPROFILE
         /**
         * @brief Do storage profiling
         */
         virtual void profileStorage() const = 0;
      #endif // GEOMHDISCC_STORAGEPROFILE
         
      protected:

      private:
   };

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> IConverter<TFwdA, TBwdA, TFwdB, TBwdB>::IConverter()
   {
      // Check that all dimensions match
      Debug::StaticAssert< (TFwdA::FieldDimension == TBwdA::FieldDimension) >();
      Debug::StaticAssert< (TBwdA::FieldDimension == TFwdB::FieldDimension) >();
      Debug::StaticAssert< (TFwdB::FieldDimension == TBwdB::FieldDimension) >();

      // Check that the data type is the same
      Debug::StaticTypeAssert<typename TFwdA::PointType , typename TBwdB::PointType>();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> IConverter<TFwdA, TBwdA, TFwdB, TBwdB>::~IConverter()
   {
   }

}
}

#endif // ICONVERTER_HPP
