/** 
 * @file DynamicPairProvider.hpp
 * @brief Templated implementation of a data pair storage provider adapting its size dynamically.
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef DYNAMICPAIRPROVIDER_HPP
#define DYNAMICPAIRPROVIDER_HPP

// System includes
//
#include <cassert>
#include <queue>
#include <list>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Exceptions/Exception.hpp"

namespace QuICC {

   /**
    * @brief Templated implementation of a data pair storage provider adapting its size dynamically.
    *
    * \tparam TForward Input type of a forward transform
    * \tparam TBackward Input type of a backward transform
    */
   template <typename TForward, typename TBackward> class DynamicPairProvider
   {
      public:
         /**
          * @brief Constructor
          */
         DynamicPairProvider();

         /**
          * @brief Destructor
          */
         ~DynamicPairProvider();

         /**
          * @brief Initialise the provider with one unit of each
          *
          * \mhdBug Should ultimatively be removed or at least simplified
          *
          * @param spSetupFwd Setup for TForward
          * @param spSetupBwd Setup for TBackward
          */
         void init(typename TForward::SharedSetupType spSetupFwd, typename TBackward::SharedSetupType spSetupBwd);

         /**
          * @brief Resize the storage
          *
          * \mhdBug Should ultimatively be removed or at least simplified
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
          * @brief Hold storage with specific ID so that it can be recovered later
          *
          * @param tmp Storage to hold (put in used map)
          * @param id  ID used as key in map
          */
         void  holdFwd(TForward &tmp, const int id);

         /**
          * @brief Hold storage with specific ID so that it can be recovered later
          *
          * @param tmp Storage to free (put in used map)
          * @param id  ID used as key in map
          */
         void  holdBwd(TBackward &tmp, const int id);

         /**
          * @brief Recover the unfreed temporary forward transform storage used previously
          */
         TForward&  recoverFwd();

         /**
          * @brief Recover the unfreed temporary backward transform storage used previously
          */
         TBackward&  recoverBwd();

         /**
          * @brief Recover the unfreed temporary forward transform storage with specific ID used previously
          */
         TForward&  recoverFwd(const int id);

         /**
          * @brief Recover the unfreed temporary backward transform storage with specifig ID used previously
          */
         TBackward&  recoverBwd(const int id);

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
          * @brief Direct assess to the vector of TForward of temporaries
          *
          * \mhdBug Should ultimatively be removed
          */
         std::list<TForward>&  rFTmps();

         /**
          * @brief Direct assess to the vector of TBackward of temporaries
          *
          * \mhdBug Should ultimatively be removed
          */
         std::list<TBackward>&  rBTmps();

     #ifdef QUICC_STORAGEPROFILE
         /**
         * @brief Get the memory requirements
         */
         MHDFloat requiredStorage() const;
     #endif // QUICC_STORAGEPROFILE
         
      protected:

      private:
         /**
          * @brief Clear all the pointer queues
          */
         void clearQueues();

         /**
          * @brief Add forward storage
          */
         void addFStorage();

         /**
          * @brief Add backward storage
          */
         void addBStorage();

         /**
          * @brief Shared forward setup
          */
         typename TForward::SharedSetupType mspFSetup;

         /**
          * @brief Shared backward setup
          */
         typename TBackward::SharedSetupType mspBSetup;

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
          * @brief In-use temporary storage map for TForward data
          */
         std::map<int,TForward *> mUsedFMap;

         /**
          * @brief In-use temporary storage map for TBackward data
          */
         std::map<int,TBackward *> mUsedBMap;

         /**
          * @brief Vector of TForward storage
          */
         std::list<TForward>  mFTmp;

         /**
          * @brief Vector of TBackward storage
          */
         std::list<TBackward>  mBTmp;

         /**
          * @brief TForward data pointer
          */
         TForward   *mpFTmp;

         /**
          * @brief TBackward data pointer
          */
         TBackward   *mpBTmp;
   };

   template <typename TForward, typename TBackward> inline std::list<TForward>& DynamicPairProvider<TForward, TBackward>::rFTmps()
   {
      return this->mFTmp;
   }

   template <typename TForward, typename TBackward> inline std::list<TBackward>& DynamicPairProvider<TForward, TBackward>::rBTmps()
   {
      return this->mBTmp;
   }

   template <typename TForward, typename TBackward> inline void DynamicPairProvider<TForward, TBackward>::freeFwd(TForward &tmp)
   {
      this->mFQueue.push(&tmp);
   }

   template <typename TForward, typename TBackward> inline void DynamicPairProvider<TForward, TBackward>::freeBwd(TBackward &tmp)
   {
      this->mBQueue.push(&tmp);
   }

   template <typename TForward, typename TBackward> inline void DynamicPairProvider<TForward, TBackward>::holdFwd(TForward &tmp)
   {
      this->mUsedFQueue.push(&tmp);
   }

   template <typename TForward, typename TBackward> inline void DynamicPairProvider<TForward, TBackward>::holdBwd(TBackward &tmp)
   {
      this->mUsedBQueue.push(&tmp);
   }

   template <typename TForward, typename TBackward> inline void DynamicPairProvider<TForward, TBackward>::holdFwd(TForward &tmp, const int id)
   {
      this->mUsedFMap.insert(std::make_pair(id, &tmp));
   }

   template <typename TForward, typename TBackward> inline void DynamicPairProvider<TForward, TBackward>::holdBwd(TBackward &tmp, const int id)
   {
      this->mUsedBMap.insert(std::make_pair(id, &tmp));
   }

   template <typename TForward, typename TBackward> TForward&  DynamicPairProvider<TForward, TBackward>::provideFwd()
   {
      // Check if storage is avaiable
      if(this->mFQueue.empty())
      {
         this->addFStorage();
      }

      // Assert for non empty storage
      assert(!this->mFQueue.empty());

      // Get pointer from forward queue
      this->mpFTmp = this->mFQueue.front();

      // Remove pointer from forward queue
      this->mFQueue.pop();

      return *this->mpFTmp;
   }

   template <typename TForward, typename TBackward> TBackward&  DynamicPairProvider<TForward, TBackward>::provideBwd()
   {
      // Check if storage is avaiable
      if(this->mBQueue.empty())
      {
         this->addBStorage();
      }

      // Assert for non empty storage
      assert(!this->mBQueue.empty());

      // Get pointer from backward queue
      this->mpBTmp = this->mBQueue.front();

      // Remove pointer from backward queue
      this->mBQueue.pop();

      return *this->mpBTmp;
   }

   template <typename TForward, typename TBackward> TForward&  DynamicPairProvider<TForward, TBackward>::recoverFwd()
   {
      // Add in an assert for non empty queue
      assert(this->mUsedFQueue.size() > 0);

      // Get pointer from used queue
      this->mpFTmp = this->mUsedFQueue.front();

      // Remove pointer from used queue
      this->mUsedFQueue.pop();

      return *this->mpFTmp;
   }

   template <typename TForward, typename TBackward> TBackward&  DynamicPairProvider<TForward, TBackward>::recoverBwd()
   {
      // Add in an assert for non empty queue
      assert(this->mUsedBQueue.size() > 0);

      // Get pointer from used queue
      this->mpBTmp = this->mUsedBQueue.front();

      // Remove pointer from used queue
      this->mUsedBQueue.pop();

      return *this->mpBTmp;
   }

   template <typename TForward, typename TBackward> TForward&  DynamicPairProvider<TForward, TBackward>::recoverFwd(const int id)
   {
      // Add in an assert for non empty map
      assert(this->mUsedFMap.count(id) > 0);

      // Get pointer from used map
      this->mpFTmp = this->mUsedFMap.find(id)->second;

      // Remove pointer from used map
      this->mUsedFMap.erase(id);

      return *this->mpFTmp;
   }

   template <typename TForward, typename TBackward> TBackward&  DynamicPairProvider<TForward, TBackward>::recoverBwd(const int id)
   {
      // Add in an assert for non empty map
      assert(this->mUsedBMap.count(id) > 0);

      // Get pointer from used map
      this->mpBTmp = this->mUsedBMap.find(id)->second;

      // Remove pointer from used map
      this->mUsedBMap.erase(id);

      return *this->mpBTmp;
   }

   template <typename TForward, typename TBackward> DynamicPairProvider<TForward, TBackward>::DynamicPairProvider()
   {
   }

   template <typename TForward, typename TBackward> DynamicPairProvider<TForward, TBackward>::~DynamicPairProvider()
   {
   }

   template <typename TForward, typename TBackward> void DynamicPairProvider<TForward, TBackward>::init(typename TForward::SharedSetupType spSetupFwd, typename TBackward::SharedSetupType spSetupBwd)
   {
      // Set the forward setup
      this->mspFSetup = spSetupFwd;

      // Set the backward setup
      this->mspBSetup = spSetupBwd;

      // Init the temporary storage queues
      this->clearQueues();

      // Add first forward storage
      this->addFStorage();

      // Add first backward storage
      this->addBStorage();
   }

   template <typename TForward, typename TBackward> void DynamicPairProvider<TForward, TBackward>::resize(const int nFTmp, const int nBTmp)
   {
      //
      // This is a simple placeholder to generalize interface for static storage providers
      //
   }

   template <typename TForward, typename TBackward> void DynamicPairProvider<TForward, TBackward>::addFStorage()
   {
      // Avoid growing to big (error in transform most likely)
      if(this->mFTmp.size() > 20)
      {
         throw Exception("Forward transform storage providers became too large!");
      }

      // Create forward storage scalar
      this->mFTmp.push_back(TForward(this->mspFSetup));

      // Add storage to forward queue
      this->mFQueue.push(&(this->mFTmp.back()));
   }

   template <typename TForward, typename TBackward> void DynamicPairProvider<TForward, TBackward>::addBStorage()
   {
      // Avoid growing to big (error in transform most likely)
      if(this->mBTmp.size() > 20)
      {
         throw Exception("Backward transform storage providers became too large!");
      }

      // Create backward storage scalar
      TBackward   tmp(this->mspBSetup);
      this->mBTmp.push_back(tmp);

      // Add storage to backward queue
      this->mBQueue.push(&(this->mBTmp.back()));
   }

   template <typename TForward, typename TBackward> void DynamicPairProvider<TForward, TBackward>::clearQueues()
   {
      // Make sure forward queue is empty
      if(!this->mFQueue.empty())
      {
         while(!this->mFQueue.empty()) this->mFQueue.pop();
      }

      // Make sure forward used queue is empty
      if(!this->mUsedFQueue.empty())
      {
         while(!this->mUsedFQueue.empty()) this->mUsedFQueue.pop();
      }

      // Make sure forward used map is empty
      this->mUsedFMap.clear();

      // Make sure backward queue is empty
      if(!this->mBQueue.empty())
      {
         while(!this->mBQueue.empty()) this->mBQueue.pop();
      }

      // Make sure backward used queue is empty
      if(!this->mUsedBQueue.empty())
      {
         while(!this->mUsedBQueue.empty()) this->mUsedBQueue.pop();
      }

      // Make sure backward used map is empty
      this->mUsedBMap.clear();
   }

   #ifdef QUICC_STORAGEPROFILE
   template <typename TForward, typename TBackward> MHDFloat  DynamicPairProvider<TForward, TBackward>::requiredStorage() const
   {
      MHDFloat mem = 0.0;

      if(this->mFTmp.size() > 0)
      {
         mem += this->mFTmp.size()*this->mFTmp.front().requiredStorage();
      }

      if(this->mBTmp.size() > 0)
      {
         mem += this->mBTmp.size()*this->mBTmp.front().requiredStorage();
      }

      return mem;
   }
   #endif // QUICC_STORAGEPROFILE

}

#endif // DYNAMICPAIRPROVIDER_HPP
