/** \file ConverterBase.hpp
 *  \brief Implementation of the base for data converter.
 */

#ifndef CONVERTERBASE_HPP
#define CONVERTERBASE_HPP

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

namespace EPMPhoenix {

   /**
    * @brief Implementation of the base for data converter.
    */
   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> class ConverterBase
   {
      public:
         /**
          * @brief Constructor
          */
         ConverterBase();

         /**
          * @brief Destructor
          */
         virtual ~ConverterBase();

         /**
          * @brief Convert data from TFwdA to TBwdB
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
          * @brief Get the converted data from TBwdA to TFwdB conversion
          */
         virtual TFwdA& getFwd(StoragePairProvider<TFwdA, TBwdA> &storage) = 0;

         /**
          * @brief Get the converted data from TFwdB to TBwdA conversion
          */
         virtual TBwdB& getBwd(StoragePairProvider<TFwdB, TBwdB> &storage) = 0;

         /**
          * @brief Setup upcoming communication
          *
          * @param packs Number of packets in communication packing
          */
         virtual void setupCommunication(const int packs) = 0;

         /**
          * @brief Start communication for forward transform
          */
         virtual void initiateForwardCommunication() = 0;

         /**
          * @brief Start communication for backward transform
          */
         virtual void initiateBackwardCommunication() = 0;

      #ifdef EPMPHOENIX_STORAGEPROFILE
         /**
         * @brief Do storage profiling
         */
         virtual void profileStorage() const = 0;
      #endif // EPMPHOENIX_STORAGEPROFILE

         /**
          * @brief Set the communication buffers pointers
          *
          * @brief fBuffers TForward communication buffers
          * @brief bBuffers TBackward communication buffers
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

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> inline const std::vector<int>& ConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::fwdSizes() const
   {
      return this->mFSizes;
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> inline const std::vector<int>& ConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::bwdSizes() const
   {
      return this->mBSizes;
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> inline int ConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::nFCpu() const
   {
      return this->mFCpuGroup.size();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> inline int ConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::nBCpu() const
   {
      return this->mBCpuGroup.size();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> inline int ConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::fCpu(const int id) const
   {
      return this->mFCpuGroup.at(id);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> inline int ConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::bCpu(const int id) const
   {
      return this->mBCpuGroup.at(id);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> inline int ConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::sizeFPacket(const int id) const
   {
      return this->mPacks*this->mFSizes.at(id);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> inline int ConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::sizeBPacket(const int id) const
   {
      return this->mPacks*this->mBSizes.at(id);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> ConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::ConverterBase()
      : mIsSending(false), mIsReceiving(false), mPacks(0)
   {
      // Check that all dimensions match
      StaticAssert< (TFwdA::FieldDimension == TBwdA::FieldDimension) >();
      StaticAssert< (TBwdA::FieldDimension == TFwdB::FieldDimension) >();
      StaticAssert< (TFwdB::FieldDimension == TBwdB::FieldDimension) >();

      // Check that the data type is the same
      StaticTypeAssert<typename TFwdA::CoefficientType , typename TBwdB::CoefficientType>();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> ConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::~ConverterBase()
   {
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> void ConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::setBuffers(std::vector<char *> &fBuffers,std::vector<char *> &bBuffers)
   {
      // Set TForward buffers pointer
      this->mpFBuffers = &fBuffers;

      // Set TBackward buffers pointer
      this->mpBBuffers = &bBuffers;
   }

}

#endif // CONVERTERBASE_HPP
