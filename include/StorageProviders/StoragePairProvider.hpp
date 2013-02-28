/** \file StoragePairProvider.hpp
 *  \brief Templated implementation of a data pair storage provider.
 *
 *  \mhdBug Needs test
 */

#ifndef STORAGEPAIRPROVIDER_HPP
#define STORAGEPAIRPROVIDER_HPP

// System includes
//
#include <assert.h>
#include <queue>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Exceptions/Exception.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief Templated implementation of a data pair storage provider.
    *
    * \tparam TForward Input type of a forward transform
    * \tparam TBackward Input type of a backward transform
    */
   template <typename TForward, typename TBackward> class StoragePairProvider
   {
      public:
         /**
          * @brief Constructor
          */
         StoragePairProvider();

         /**
          * @brief Destructor
          */
         virtual ~StoragePairProvider();

         /**
          * @brief Initialise the provider with one unit of each
          *
          * @param setupFwd Setup for TForward
          * @param setupBwd Setup for TBackward
          */
         void init(const typename TForward::SetupType& setupFwd, const typename TBackward::SetupType& setupBwd);

         /**
          * @brief Resize the storage
          *
          * @param nFTmp Number of required forward storage units
          * @param nBTmp Number of required backward storage units
          */
         void resize(const int nFTmp, const int nBTmp);

         /**
          * @brief Provide temporary storage for forward transform
          */
         TForward&  provideFwd();

         /**
          * @brief Provide temporary storage for backward transform
          */
         TBackward&  provideBwd();

         /**
          * @brief Hold storage so that it can be recovered later
          *
          * @param tmp Storage to hold (put in used queue)
          */
         void  holdFwd(TForward &tmp);

         /**
          * @brief Hold storage so that it can be recovered later
          *
          * @param tmp Storage to free (put in used queue)
          */
         void  holdBwd(TBackward &tmp);

         /**
          * @brief Recover the unfreed temporary forward transform storage used previously
          */
         TForward&  recoverFwd();

         /**
          * @brief Recover the unfreed temporary backward transform storage used previously
          */
         TBackward&  recoverBwd();

         /**
          * @brief Free tempory storage after use and put back into queue
          *
          * @param tmp Storage to free (put back in queue)
          */
         void  freeFwd(TForward &tmp);

         /**
          * @brief Free tempory storage after use and put back into queue
          *
          * @param tmp Storage to free (put back in queue)
          */
         void  freeBwd(TBackward &tmp);

         /**
          * @brief Number of TForward temporaries
          */
         int nFTmp() const;

         /**
          * @brief Number of TBackward temporaries
          */
         int nBTmp() const;

         /**
          * @brief Direct access to the temporaries
          *
          * @param i Index of the temporary
          */
         TForward& rFTmp(const int i);

         /**
          * @brief Direct access to the temporaries
          *
          * @param i Index of the temporary
          */
         TBackward& rBTmp(const int i);

         /**
          * @brief Direct assess to the vector of TForward of temporaries
          */
         std::vector<TForward>&  rFTmps();

         /**
          * @brief Direct assess to the vector of TBackward of temporaries
          */
         std::vector<TBackward>&  rBTmps();

     #ifdef GEOMHDISCC_STORAGEPROFILE
         /**
         * @brief Get the memory requirements
         */
         MHDFloat requiredStorage() const;
     #endif // GEOMHDISCC_STORAGEPROFILE
         
      protected:

         /**
          * @brief Number of forward tempory storage scalars
          */
         int mNFTmp;

         /**
          * @brief Number of backward tempory storage scalars
          */
         int mNBTmp;

         /**
          * @brief Available temporary storage queue for TForward data
          */
         std::queue<TForward *> mFQueue;

         /**
          * @brief Available temporary storage queue for TBackward data
          */
         std::queue<TBackward *> mBQueue;

         /**
          * @brief In-use temporary storage queue for TForward data
          */
         std::queue<TForward *> mUsedFQueue;

         /**
          * @brief In-use temporary storage queue for TBackward data
          */
         std::queue<TBackward *> mUsedBQueue;

         /**
          * @brief Vector of TForward storage
          */
         std::vector<TForward>  mFTmp;

         /**
          * @brief Vector of TBackward storage
          */
         std::vector<TBackward>  mBTmp;

         /**
          * @brief TForward data pointer
          */
         TForward   *mpFTmp;

         /**
          * @brief TBackward data pointer
          */
         TBackward   *mpBTmp;

      private:
         /**
          * @brief Initialise temporary storage
          *
          * @param setupFwd Setup for TForward
          * @param setupBwd Setup for TBackward
          */
         void initStorage(const typename TForward::SetupType& setupFwd, const typename TBackward::SetupType& setupBwd);

         /**
          * @brief Initialise temporary storage queues (available storage)
          */
         void initQueues();
   };

   template <typename TForward, typename TBackward> inline int StoragePairProvider<TForward, TBackward>::nFTmp() const
   {
      return this->mNFTmp;
   }

   template <typename TForward, typename TBackward> inline int StoragePairProvider<TForward, TBackward>::nBTmp() const
   {
      return this->mNBTmp;
   }

   template <typename TForward, typename TBackward> inline TForward& StoragePairProvider<TForward, TBackward>::rFTmp(const int i)
   {
      return this->mFTmp.at(i);
   }

   template <typename TForward, typename TBackward> inline TBackward& StoragePairProvider<TForward, TBackward>::rBTmp(const int i)
   {
      return this->mBTmp.at(i);
   }

   template <typename TForward, typename TBackward> inline std::vector<TForward>& StoragePairProvider<TForward, TBackward>::rFTmps()
   {
      return this->mFTmp;
   }

   template <typename TForward, typename TBackward> inline std::vector<TBackward>& StoragePairProvider<TForward, TBackward>::rBTmps()
   {
      return this->mBTmp;
   }

   template <typename TForward, typename TBackward> inline void StoragePairProvider<TForward, TBackward>::freeFwd(TForward &tmp)
   {
      this->mFQueue.push(&tmp);
   }

   template <typename TForward, typename TBackward> inline void StoragePairProvider<TForward, TBackward>::freeBwd(TBackward &tmp)
   {
      this->mBQueue.push(&tmp);
   }

   template <typename TForward, typename TBackward> inline void StoragePairProvider<TForward, TBackward>::holdFwd(TForward &tmp)
   {
      this->mUsedFQueue.push(&tmp);
   }

   template <typename TForward, typename TBackward> inline void StoragePairProvider<TForward, TBackward>::holdBwd(TBackward &tmp)
   {
      this->mUsedBQueue.push(&tmp);
   }

   template <typename TForward, typename TBackward> TForward&  StoragePairProvider<TForward, TBackward>::recoverFwd()
   {
      // Add in an assert for non empty queue
      assert(this->mUsedFQueue.size());

      this->mpFTmp = this->mUsedFQueue.front();
      this->mUsedFQueue.pop();

      return *this->mpFTmp;
   }

   template <typename TForward, typename TBackward> TBackward&  StoragePairProvider<TForward, TBackward>::recoverBwd()
   {
      // Add in an assert for non empty queue
      assert(this->mUsedBQueue.size());

      this->mpBTmp = this->mUsedBQueue.front();
      this->mUsedBQueue.pop();

      return *this->mpBTmp;
   }

   template <typename TForward, typename TBackward> StoragePairProvider<TForward, TBackward>::StoragePairProvider()
      : mNFTmp(0), mNBTmp(0)
   {
   }

   template <typename TForward, typename TBackward> void StoragePairProvider<TForward, TBackward>::init(const typename TForward::SetupType& setupFwd, const typename TBackward::SetupType& setupBwd)
   {
      // Set initial number of TForward storages
      this->mNFTmp = 2;

      // Set initial number of TBackward storages
      this->mNBTmp = 2;

      // Init the temporary storage
      this->initStorage(setupFwd, setupBwd);

      // Init the temporary storage queues
      this->initQueues();
   }

   template <typename TForward, typename TBackward> void StoragePairProvider<TForward, TBackward>::resize(const int nFTmp, const int nBTmp)
   {
      // Set number of TForward storages
      this->mNFTmp = nFTmp;

      // Set number of TBackward storages
      this->mNBTmp = nBTmp;

      // Make sure storage has been initialised
      if(this->mFTmp.size() == 0 || this->mBTmp.size())
      {
         throw Exception("StoragePairProvider::resize", "Can't resize uninitialised memory storage");
      } else
      {
         // Resize TForward storage
         this->mFTmp.resize(this->mNFTmp, this->mFTmp.at(0));

         // Resize TBackward storage
         this->mBTmp.resize(this->mNBTmp, this->mBTmp.at(0));
      }

      // Update queue
      this->initQueues();
   }

   template <typename TForward, typename TBackward> void StoragePairProvider<TForward, TBackward>::initStorage(const typename TForward::SetupType& setupFwd, const typename TBackward::SetupType& setupBwd)
   {
      // Initialise the TForward storage data
      for(int i=0; i < this->mNFTmp; ++i)
      {
         this->mFTmp.push_back(TForward(setupFwd));
      }

      // Initialise the TBackward storage data
      for(int i=0; i < this->mNBTmp; ++i)
      {
         this->mBTmp.push_back(TBackward(setupBwd));
      }
   }

   template <typename TForward, typename TBackward> void StoragePairProvider<TForward, TBackward>::initQueues()
   {
      // Make sure queue is empty
      if(!this->mFQueue.empty())
      {
         while(!this->mFQueue.empty()) this->mFQueue.pop();
      }

      // Initialise the TForward storage queue
      for(int i=0; i < this->mNFTmp; ++i)
      {
         this->mFQueue.push(&(this->mFTmp.at(i)));
      }

      // Make sure queue is empty
      if(!this->mBQueue.empty())
      {
         while(!this->mBQueue.empty()) this->mBQueue.pop();
      }

      // Initialise the TBackward storage queue
      for(int i=0; i < this->mNBTmp; ++i)
      {
         this->mBQueue.push(&(this->mBTmp.at(i)));
      }
   }

   template <typename TForward, typename TBackward> TForward&  StoragePairProvider<TForward, TBackward>::provideFwd()
   {
      // Add in an assert for non empty queue
      assert(this->mFQueue.size());

      this->mpFTmp = this->mFQueue.front();
      this->mFQueue.pop();

      return *this->mpFTmp;
   }

   template <typename TForward, typename TBackward> TBackward&  StoragePairProvider<TForward, TBackward>::provideBwd()
   {
      // Add in an assert for non empty queue
      assert(this->mBQueue.size());

      this->mpBTmp = this->mBQueue.front();
      this->mBQueue.pop();

      return *this->mpBTmp;
   }

   #ifdef GEOMHDISCC_STORAGEPROFILE
   template <typename TForward, typename TBackward> MHDFloat  StoragePairProvider<TForward, TBackward>::requiredStorage() const
   {
      MHDFloat mem = 0.0;

      if(this->nFTmp() > 0)
      {
         mem += this->nFTmp()*this->mFTmp.at(0).requiredStorage();
      }

      if(this->nBTmp() > 0)
      {
         mem += this->nBTmp()*this->mBTmp.at(0).requiredStorage();
      }

      return mem;
   }
   #endif // GEOMHDISCC_STORAGEPROFILE

}

#endif // STORAGEPAIRPROVIDER_HPP
