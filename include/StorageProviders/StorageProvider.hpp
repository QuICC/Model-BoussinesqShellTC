/** 
 * @file StorageProvider.hpp
 * @brief Templated implementation of a single data storage provider.
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef STORAGEPROVIDER_HPP
#define STORAGEPROVIDER_HPP

// System includes
//
#include <cassert>
#include <queue>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"

namespace QuICC {

   /**
    * @brief Templated implementation of a single data storage provider.
    *
    * \tparam TData Data type
    */
   template <typename TData> class StorageProvider
   {
      public:
         /**
          * @brief Constructor
          */
         StorageProvider();

         /**
          * @brief Destructor
          */
         ~StorageProvider();

         /**
          * @brief Initialise the provider
          *
          * @param spSetup Shared setup for TData
          */
         void init(typename TData::SharedSetupType spSetup);

         /**
          * @brief Resize the storage
          *
          * @param nTmp Number of required storage units
          */
         void resize(const int nTmp);

         /**
          * @brief Provide temporary storage
          */
         TData&  provide();

         /**
          * @brief Hold storage so that it can be recovered later
          *
          * @param tmp Storage to hold (put in used queue)
          */
         void  hold(TData &tmp);

         /**
          * @brief Recover the unfreed temporary storage
          */
         TData&  recover();

         /**
          * @brief Free tempory storage after use and put back into queue
          *
          * @param tmp Storage to free (put back in queue)
          */
         void  free(TData &tmp);

         /**
          * @brief Number of temporaries
          */
         int nTmp() const;

         /**
          * @brief Direct access to the temporaries
          *
          * @param i Index of the temporary
          */
         TData& rTmp(const int i);

         /**
          * @brief Direct assess to the vector of temporaries
          */
         std::vector<TData>&  rTmps();

     #ifdef QUICC_STORAGEPROFILE
         /**
         * @brief Get the memory requirements
         */
         MHDFloat requiredStorage() const;
     #endif // QUICC_STORAGEPROFILE
         
      protected:

         /**
          * @brief Number of forward tempory storage scalars
          */
         int mNTmp;

         /**
          * @brief Available temporary storage queue
          */
         std::queue<TData *> mQueue;

         /**
          * @brief In-use temporary storage queue
          */
         std::queue<TData *> mUsedQueue;

         /**
          * @brief Vector of TData storage
          */
         std::vector<TData>  mTmp;

         /**
          * @brief TData data pointer
          */
         TData   *mpTmp;

      private:
         /**
          * @brief Initialise temporary storage
          *
          * @param spSetup Shared setup for TData
          */
         void initStorage(typename TData::SharedSetupType spSetup);

         /**
          * @brief Initialise temporary storage queues (available storage)
          */
         void initQueues();
   };

   template <typename TData> inline int StorageProvider<TData>::nTmp() const
   {
      return this->mNTmp;
   }

   template <typename TData> inline TData& StorageProvider<TData>::rTmp(const int i)
   {
      return this->mTmp.at(i);
   }

   template <typename TData> inline std::vector<TData>& StorageProvider<TData>::rTmps()
   {
      return this->mTmp;
   }

   template <typename TData> inline void StorageProvider<TData>::free(TData &tmp)
   {
      this->mQueue.push(&tmp);
   }

   template <typename TData> inline void StorageProvider<TData>::hold(TData &tmp)
   {
      this->mUsedQueue.push(&tmp);
   }

   template <typename TData> TData&  StorageProvider<TData>::recover()
   {
      // Add in an assert for non empty queue
      assert(this->mUsedQueue.size());

      this->mpTmp = this->mUsedQueue.front();
      this->mUsedQueue.pop();

      return *this->mpTmp;
   }

   template <typename TData> StorageProvider<TData>::StorageProvider()
      : mNTmp(0)
   {
   }

   template <typename TData> StorageProvider<TData>::~StorageProvider()
   {
   }

   template <typename TData> void StorageProvider<TData>::init(typename TData::SharedSetupType spSetup)
   {
      // Set initial number of TData storages
      this->mNTmp = 1;

      // Init the temporary storage
      this->initStorage(spSetup);

      // Init the temporary storage queues
      this->initQueues();
   }

   template <typename TData> void StorageProvider<TData>::resize(const int nTmp)
   {
      // Set number of TData storages
      this->mNTmp = nTmp;

      // Make sure storage has been initialised
      if(this->mTmp.size() == 0)
      {
         throw Exception("StorageProvider::resize", "Can't resize uninitialised memory storage");
      } else
      {
         // Resize TData storage
         this->mTmp.resize(this->mNTmp, this->mTmp.at(0));
      }

      // Init the temporary storage queues
      this->initQueues();
   }

   template <typename TData> void StorageProvider<TData>::initStorage(typename TData::SharedSetupType spSetup)
   {
      // Initialise the TData storage data
      for(int i=0; i < this->mNTmp; ++i)
      {
         this->mTmp.push_back(TData(spSetup));
      }
   }

   template <typename TData> void StorageProvider<TData>::initQueues()
   {
      // Make sure queue is empty
      if(!this->mQueue.empty())
      {
         while(!this->mQueue.empty()) this->mQueue.pop();
      }

      // Initialise the TData storage queue
      for(int i=0; i < this->mNTmp; ++i)
      {
         this->mQueue.push(&(this->mTmp.at(i)));
      }
   }

   template <typename TData> TData&  StorageProvider<TData>::provide()
   {
      // Add in an assert for non empty queue
      assert(this->mBQueue.size());

      this->mpTmp = this->mQueue.front();
      this->mQueue.pop();

      return *this->mpTmp;
   }

   #ifdef QUICC_STORAGEPROFILE
   template <typename TData> MHDFloat  StorageProvider<TData>::required() const
   {
      MHDFloat mem = 0.0;

      if(this->nTmp() > 0)
      {
         mem += this->nTmp()*this->mTmp.at(0).requiredStorage();
      }

      return mem;
   }
   #endif // QUICC_STORAGEPROFILE

}

#endif // STORAGEPROVIDER_HPP
