/** \file IConverter.hpp
 *  \brief Implementation of the interface for a data converter.
 *
 *  \mhdBug Needs test
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
#include "StorageProviders/StoragePairProvider.hpp"

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
         virtual void convertFwd(const TFwdA &in, StoragePairProvider<TFwdB, TBwdB> &storage) = 0;

         /**
          * @brief Convert data from TBwdB to TFwdA
          */
         virtual void convertBwd(const TBwdB &in, StoragePairProvider<TFwdA, TBwdA> &storage) = 0;

         /**
          * @brief Get the converted data from TBwdB to TFwdA conversion
          */
         virtual TFwdA& getFwd(StoragePairProvider<TFwdA, TBwdA> &storage) = 0;

         /**
          * @brief Get the converted data from TFwdA to TBwdB conversion
          */
         virtual TBwdB& getBwd(StoragePairProvider<TFwdB, TBwdB> &storage) = 0;

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

         /**
          * @brief Set the communication buffers pointers
          *
          * @brief fBuffers Forward communication buffers
          * @brief bBuffers Backward communication buffers
          */
         void setBuffers(std::vector<char *> &fBuffers, std::vector<char *> &bBuffers);

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
          * @brief Get global MPI rank of CPU from forward CPU group
          *
          * @param id CPU group id
          */
         int fCpu(const int id) const;

         /**
          * @brief Get global MPI rank of CPU from backward CPU group
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
          * @brief List of CPU ranks involved in the forward conversion
          */
         std::vector<int>  mFCpuGroup;

         /**
          * @brief List of CPU ranks involved in the backward conversion
          */
         std::vector<int>  mBCpuGroup;

         /**
          * @brief List of the forward buffer sizes
          */
         std::vector<int>  mFSizes;

         /**
          * @brief List of the backward buffer sizes
          */
         std::vector<int>  mBSizes;

         /**
          * @brief Forward communication buffer pointer
          */
         std::vector<char *> *mpFBuffers;

         /**
          * @brief Backward communication buffer pointer
          */
         std::vector<char *> *mpBBuffers;

         /**
          * @brief Communication packs counter
          */
         int mPacks;

         /**
          * @brief Possible forward transform packs
          */
         ArrayI   mForwardPacks;

         /**
          * @brief Possible backward transform packs
          */
         ArrayI   mBackwardPacks;

      private:
   };

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> inline const std::vector<int>& IConverter<TFwdA, TBwdA, TFwdB, TBwdB>::fwdSizes() const
   {
      return this->mFSizes;
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> inline const std::vector<int>& IConverter<TFwdA, TBwdA, TFwdB, TBwdB>::bwdSizes() const
   {
      return this->mBSizes;
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> inline int IConverter<TFwdA, TBwdA, TFwdB, TBwdB>::nFCpu() const
   {
      return this->mFCpuGroup.size();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> inline int IConverter<TFwdA, TBwdA, TFwdB, TBwdB>::nBCpu() const
   {
      return this->mBCpuGroup.size();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> inline int IConverter<TFwdA, TBwdA, TFwdB, TBwdB>::fCpu(const int id) const
   {
      return this->mFCpuGroup.at(id);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> inline int IConverter<TFwdA, TBwdA, TFwdB, TBwdB>::bCpu(const int id) const
   {
      return this->mBCpuGroup.at(id);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> inline int IConverter<TFwdA, TBwdA, TFwdB, TBwdB>::sizeFPacket(const int id) const
   {
      return this->mPacks*this->mFSizes.at(id);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> inline int IConverter<TFwdA, TBwdA, TFwdB, TBwdB>::sizeBPacket(const int id) const
   {
      return this->mPacks*this->mBSizes.at(id);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> IConverter<TFwdA, TBwdA, TFwdB, TBwdB>::IConverter()
      : mIsSending(false), mIsReceiving(false), mPacks(0)
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

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> void IConverter<TFwdA, TBwdA, TFwdB, TBwdB>::setBuffers(std::vector<char *> &fBuffers,std::vector<char *> &bBuffers)
   {
      // Set TForward buffers pointer
      this->mpFBuffers = &fBuffers;

      // Set TBackward buffers pointer
      this->mpBBuffers = &bBuffers;
   }

}
}

#endif // ICONVERTER_HPP
