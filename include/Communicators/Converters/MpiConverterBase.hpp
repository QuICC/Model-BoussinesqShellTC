/** \file MpiConverterBase.hpp
 *  \brief Templated implementation of the base of a MPI data converter.
 */

#ifndef MPICONVERTERBASE_HPP
#define MPICONVERTERBASE_HPP

// Configuration includes
//
#include "Framework/FrameworkMacro.h"
#include "StorageProfiler/StorageProfilerMacro.h"
#include "Profiler/ProfilerMacro.h"

// System includes
//
#include <assert.h>
#include <set>
#include <map>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MpiTypes.hpp"
#include "Resolutions/Resolution.hpp"
#include "Communicators/Converters/IConverter.hpp"

namespace GeoMHDiSCC {

namespace Parallel {

   /**
    * \brief Templated implementation of the base of a MPI data converter.
    */
   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> class MpiConverterBase: public IConverter<TFwdA,TBwdA,TFwdB,TBwdB>
   {
      public:
         /**
          * @brief Constructor
          */
         MpiConverterBase();

         /**
          * @brief Destructor
          */
         virtual ~MpiConverterBase();
         
      protected:
         /**
          * @brief Build an MPI Datatype
          *
          * @param spRes   Shared Resolution
          * @param fwdDim  Dimension index for forward transform
          * @param data    Input data
          * @param type    Created MPI data type
          * @param cpuId   ID of the CPU
          */
         void buildFwdDatatype(SharedResolution spRes, const int fwdDim, TFwdA &data, MPI_Datatype &type, const int cpuId);

         /**
          * @brief Build an MPI Datatype
          *
          * @param spRes   Shared Resolution
          * @param fwdDim  Dimension index for forward transform
          * @param data    Input data
          * @param type    Created MPI data type
          * @param cpuId   ID of the CPU
          */
         void buildBwdDatatype(SharedResolution spRes, const int fwdDim, TBwdB &data, MPI_Datatype &type, const int cpuId);

         /**
          * @brief Extract shared indexes
          *
          * @param sharedMap     Storage for the shared index map
          * @param localIdxMap   Local node key to indexes map 
          * @param remoteKeys    Remote node index keys
          */
         void extractShared(std::map<std::tr1::tuple<int,int,int>,std::tr1::tuple<int,int,int> >& sharedMap, const std::map<std::tr1::tuple<int,int,int>,std::tr1::tuple<int,int,int> >& localIdxMap, const std::set<std::tr1::tuple<int,int,int> >& remoteKeys);

         /**
          * @brief Create type
          *
          * @param data       The concerned data
          * @param sharedMap  Shared index map
          * @param type       MPI datatype storage
          */
         template <typename TData> void buildType(TData &data, const std::map<std::tr1::tuple<int,int,int>,std::tr1::tuple<int,int,int> >& sharedMap, MPI_Datatype &type);

         /**
          * @brief Reset Receive positions
          */
         void resetRecvPositions();

         /**
          * @brief Reset Send positions
          */
         void resetSendPositions();

         /**
          * @brief Get a pointer to the receive backward requests
          *
          * @param size Pack size of the requested request
          */
         MPI_Request * pRecvBRequests(const int size);

         /**
          * @brief Get a pointer to the receive forward requests
          *
          * @param size Pack size of the requested request
          */
         MPI_Request * pRecvFRequests(const int size);

         /**
          * @brief Get a pointer to the send backward requests
          *
          * @param size Pack size of the requested request
          */
         MPI_Request * pSendBRequests(const int size);

         /**
          * @brief Get a pointer to the send forward requests
          *
          * @param size Pack size of the requested request
          */
         MPI_Request * pSendFRequests(const int size);

         /**
          * @brief Get ring recv source for id
          *
          * @param id   ID of the CPU
          * @param ref  ID of the reference CPU
          * @param size Size of the CPU group
          */
         int  recvSrc(const int id, const int ref, const int size) const;

         /**
          * @brief Get ring send destination for id
          *
          * @param id   ID of the CPU
          * @param ref  ID of the reference CPU
          * @param size Size of the CPU group
          */
         int  sendDest(const int id, const int ref, const int size) const;

         /**
          * @brief Initialise the positions
          */
         void initPositions();

         /**
          * @brief Setup the MPI communication requests requests
          */
         void setupRequests();

         /**
          * @brief Cleanup the MPI communication requests
          */
         void cleanupRequests();

         /**
          * @brief The number of packs in the "previous/active" forward send
          */
         int   mActiveFSendPacks;

         /**
          * @brief The number of packs in the "previous/active" backward send
          */
         int   mActiveBSendPacks;

         /**
          * @brief Storage for the receive position pointers
          */
         std::vector<int>  mRecvPositions;

         /**
          * @brief Storage for the send position pointers
          */
         std::vector<int>  mSendPositions;

         /**
          * @brief Storage for the non blocking communication requests: Recv F
          */
         std::map<int, std::vector<MPI_Request> >  mRecvFRequests;

         /**
          * @brief Storage for the non blocking communication requests: Recv B
          */
         std::map<int, std::vector<MPI_Request> >  mRecvBRequests;

         /**
          * @brief Storage for the non blocking communication requests: Send F
          */
         std::map<int, std::vector<MPI_Request> >  mSendFRequests;

         /**
          * @brief Storage for the non blocking communication requests: Send B
          */
         std::map<int, std::vector<MPI_Request> >  mSendBRequests;

      private:
   };

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> inline MPI_Request * MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::pRecvBRequests(const int size)
   {
      return &(this->mRecvBRequests[size].front());
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> inline MPI_Request * MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::pRecvFRequests(const int size)
   {
      return &(this->mRecvFRequests[size].front());
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> inline MPI_Request * MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::pSendBRequests(const int size)
   {
      return &(this->mSendBRequests[size].front());
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> inline MPI_Request * MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::pSendFRequests(const int size)
   {
      return &(this->mSendFRequests[size].front());
   }
      
   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::MpiConverterBase()
      : IConverter<TFwdA, TBwdA, TFwdB, TBwdB>(), mActiveFSendPacks(0), mActiveBSendPacks(0)
   {
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::~MpiConverterBase()
   {
      // Cleanup the requests memory
      this->cleanupRequests();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> void MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::extractShared(std::map<std::tr1::tuple<int,int,int>,std::tr1::tuple<int,int,int> >& sharedMap, const std::map<std::tr1::tuple<int,int,int>,std::tr1::tuple<int,int,int> >& localIdxMap, const std::set<std::tr1::tuple<int,int,int> >& remoteKeys)
   {
      // List of local index keys 
      std::set<std::tr1::tuple<int,int,int> >  localKeys;
      // Create map iterator
      std::map<std::tr1::tuple<int,int,int>,std::tr1::tuple<int,int,int> >::const_iterator mapIt;
      
      // Extract the set of local keys
      for(mapIt = localIdxMap.begin(); mapIt != localIdxMap.end(); ++mapIt)
      {
         localKeys.insert(mapIt->first);
      }
      
      // Storage for the shared keys
      std::set<std::tr1::tuple<int,int,int> > sharedKeys;

      // Create the list of common indexes
      std::set_intersection(localKeys.begin(), localKeys.end(), remoteKeys.begin(), remoteKeys.end(), std::inserter(sharedKeys, sharedKeys.begin()));

      // Clear the shared map
      sharedMap.clear();

      // Fill shared map
      std::set<std::tr1::tuple<int,int,int> >::iterator sit;
      for(sit = sharedKeys.begin(); sit != sharedKeys.end(); sit++)
      {
         sharedMap.insert(std::make_pair(*sit, localIdxMap.find(*sit)->second));
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> template <typename TData> void MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::buildType(TData &data, const std::map<std::tr1::tuple<int,int,int>,std::tr1::tuple<int,int,int> >& sharedMap, MPI_Datatype &type)
   {
      // Get the number of elements
      int nElements = sharedMap.size();
      // Prepare data required for MPI datatype
      MPI_Aint    displ[nElements];
      int         blocks[nElements];

      // Prepare data required for MPI datatype
      MPI_Aint    base;
      MPI_Aint    element;

      // Set the base of the datatype to 0 (thus using the memory "origin", this is the safest approach that I could make to work)
      base = 0;

      // Create iterator
      std::map<std::tr1::tuple<int,int,int>,std::tr1::tuple<int,int,int> >::const_iterator it;

      // Create MPI displacement list
      int tot = 0;
      for(it = sharedMap.begin(); it != sharedMap.end(); ++it)
      {
         // Get address of stored coordinates
         MPI_Get_address(&data.rPoint(std::tr1::get<0>(it->second), std::tr1::get<1>(it->second), std::tr1::get<2>(it->second)), &element);

         // Fill datatype information
         displ[tot] = element - base;
         blocks[tot] = 1;

         // Increment datatype size
         tot++;
      }

      // Create MPI datatype
      MPI_Type_create_hindexed(nElements, blocks, displ, MpiTypes::type<typename TFwdA::CoefficientType>(), &type);
      // Commit MPI datatype
      MPI_Type_commit(&type);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> void MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::buildFwdDatatype(SharedResolution spRes, const int fwdDim, TFwdA &data, MPI_Datatype &type, const int cpuId)
   {
      // Create  map of the local indexes to unique keys
      std::map<std::tr1::tuple<int,int,int>,std::tr1::tuple<int,int,int> >  localIdxMap;
      // List of remote keys
      std::set<std::tr1::tuple<int,int,int> >  remoteKeys;

      // Storage for the index
      std::tr1::tuple<int,int,int> coord;
      
      // Storage for the simulation wide indexes
      int i_, j_, k_;

      // Storage for the generated key
      std::tr1::tuple<int,int,int> key;

      // Setup in the 3D case
      if(TFwdB::FieldDimension == Dimensions::THREED)
      {
         // Create the list of local indexes
         for(int k=0; k < spRes->cpu(FrameworkMacro::id())->dim(fwdDim)->dim3D(); ++k)
         {
            // Extract "physical" index of third dimension
            k_ = spRes->cpu(FrameworkMacro::id())->dim(fwdDim)->idx3D(k);

            for(int j=0; j < spRes->cpu(FrameworkMacro::id())->dim(fwdDim)->dim2D(k); ++j)
            {
               // Extract "physical" index of second dimension
               j_ = spRes->cpu(FrameworkMacro::id())->dim(fwdDim)->idx2D(j,k);

               for(int i=0; i < spRes->cpu(FrameworkMacro::id())->dim(fwdDim)->dimFwd(j,k); ++i)
               {
                  // Extract "physical" index of first dimension
                  i_ = spRes->cpu(FrameworkMacro::id())->dim(fwdDim)->idxFwd(i,j,k);

                  // Combine array indexes into coordinate tuple
                  coord = std::tr1::make_tuple(i, j, k);

                  // Create key as (1D, 3D, 2D)
                  key = std::tr1::make_tuple(i_, j_, k_);

                  // add key->coordinate to map
                  localIdxMap.insert(std::pair<std::tr1::tuple<int,int,int>,std::tr1::tuple<int,int,int> >(key, coord));
               }
            }
         }

         // Create the list of remote indexes
         for(int k=0; k < spRes->cpu(cpuId)->dim(fwdDim+1)->dim3D(); ++k)
         {
            // Extract "physical" index of third dimension
            k_ = spRes->cpu(cpuId)->dim(fwdDim+1)->idx3D(k);

            for(int j=0; j < spRes->cpu(cpuId)->dim(fwdDim+1)->dim2D(k); ++j)
            {
               // Extract "physical" index of second dimension
               j_ = spRes->cpu(cpuId)->dim(fwdDim+1)->idx2D(j,k);

               for(int i=0; i < spRes->cpu(cpuId)->dim(fwdDim+1)->dimBwd(j,k); ++i)
               {
                  // Extract "physical" index of first dimension
                  i_ = spRes->cpu(cpuId)->dim(fwdDim+1)->idxBwd(i,j,k);

                  // Create key as (1D, 3D, 2D)
                  key = std::tr1::make_tuple(j_, k_, i_);

                  // Add key to remote set
                  remoteKeys.insert(key);
               }
            }
         }

      // Setup in the 2D case
      } else if(TFwdB::FieldDimension == Dimensions::TWOD)
      {
         // Create the list of local indexes
         for(int j=0; j < spRes->cpu(FrameworkMacro::id())->dim(fwdDim)->dim2D(); ++j)
         {
            // Extract "physical" index of second dimension
            j_ = spRes->cpu(FrameworkMacro::id())->dim(fwdDim)->idx2D(j);

            for(int i=0; i < spRes->cpu(FrameworkMacro::id())->dim(fwdDim)->dimFwd(j); ++i)
            {
               // Extract "physical" index of first dimension
               i_ = spRes->cpu(FrameworkMacro::id())->dim(fwdDim)->idxFwd(i,j);

               // Combine array indexes into coordinate tuple
               coord = std::tr1::make_tuple(i, j, 0);

               // Create key as (1D, 2D, 0)
               key = std::tr1::make_tuple(i_, j_, 0);

               // add key->coordinate to map
               localIdxMap.insert(std::pair<std::tr1::tuple<int,int,int>,std::tr1::tuple<int,int,int> >(key, coord));
            }
         }

         // Create the list of remote indexes
         for(int j=0; j < spRes->cpu(cpuId)->dim(fwdDim)->dim2D(); ++j)
         {
            // Extract "physical" index of second dimension
            j_ = spRes->cpu(cpuId)->dim(fwdDim)->idx2D(j);

            for(int i=0; i < spRes->cpu(cpuId)->dim(fwdDim+1)->dimBwd(j); ++i)
            {
               // Extract "physical" index of first dimension
               i_ = spRes->cpu(cpuId)->dim(fwdDim+1)->idxBwd(i,j);

               // Create key as (1D, 2D, 0)
               key = std::tr1::make_tuple(j_, i_, 0);

               // Add key to remote set
               remoteKeys.insert(key);
            }
         }

      // Setup in the 1D case
      } else if(TFwdB::FieldDimension == Dimensions::ONED)
      {
         // Create the list of local indexes
         for(int i=0; i < spRes->cpu(FrameworkMacro::id())->dim(fwdDim)->dimFwd(); ++i)
         {
            // Extract "physical" index of first dimension
            i_ = spRes->cpu(FrameworkMacro::id())->dim(fwdDim)->idxFwd(i);

            // Combine array indexes into coordinate tuple
            coord = std::tr1::make_tuple(i, 0, 0);

            // Create key as (1D, 0, 0)
            key = std::tr1::make_tuple(i_, 0, 0);

            // add key->coordinate to map
            localIdxMap.insert(std::pair<std::tr1::tuple<int,int,int>,std::tr1::tuple<int,int,int> >(key, coord));
         }

         // Create the list of remote indexes
         for(int i=0; i < spRes->cpu(cpuId)->dim(fwdDim+1)->dimBwd(); ++i)
         {
            // Extract "physical" index of first dimension
            i_ = spRes->cpu(cpuId)->dim(fwdDim+1)->idxBwd(i);

            // Create key as (1D, 0, 0)
            key = std::tr1::make_tuple(i_, 0, 0);

            // Add key to remote set
            remoteKeys.insert(key);
         }
      }

      // Extract map of shared indexes (stored as keys)
      std::map<std::tr1::tuple<int,int,int>,std::tr1::tuple<int,int,int> >  sharedMap;
      this->extractShared(sharedMap, localIdxMap, remoteKeys);

      // Create the datatype
      this->buildType(data, sharedMap, type);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> void MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::buildBwdDatatype(SharedResolution spRes, const int fwdDim, TBwdB &data, MPI_Datatype &type, const int cpuId)
   {
      // Create  map of the local indexes to unique keys
      std::map<std::tr1::tuple<int,int,int>,std::tr1::tuple<int,int,int> >  localIdxMap;
      // List of remote keys
      std::set<std::tr1::tuple<int,int,int> >  remoteKeys;

      // Storage for the index
      std::tr1::tuple<int,int,int> coord;
      
      // Storage for the simulation wide indexes
      int i_, j_, k_;

      // Storage for the generated key
      std::tr1::tuple<int,int,int> key;

      // Setup in the 3D case
      if(TFwdB::FieldDimension == Dimensions::THREED)
      {
         // Create the list of local indexes
         for(int k=0; k < spRes->cpu(FrameworkMacro::id())->dim(fwdDim+1)->dim3D(); ++k)
         {
            // Extract "physical" index of third dimension
            k_ = spRes->cpu(FrameworkMacro::id())->dim(fwdDim+1)->idx3D(k);

            for(int j=0; j < spRes->cpu(FrameworkMacro::id())->dim(fwdDim+1)->dim2D(k); ++j)
            {
               // Extract "physical" index of second dimension
               j_ = spRes->cpu(FrameworkMacro::id())->dim(fwdDim+1)->idx2D(j,k);

               for(int i=0; i < spRes->cpu(FrameworkMacro::id())->dim(fwdDim+1)->dimBwd(j,k); ++i)
               {
                  // Extract "physical" index of first dimension
                  i_ = spRes->cpu(FrameworkMacro::id())->dim(fwdDim+1)->idxBwd(i,j,k);

                  // Combine array indexes into coordinate tuple
                  coord = std::tr1::make_tuple(i, j, k);

                  // Create key as (1D, 3D, 2D)
                  key = std::tr1::make_tuple(j_, k_, i_);

                  // add key->coordinate to map
                  localIdxMap.insert(std::pair<std::tr1::tuple<int,int,int>,std::tr1::tuple<int,int,int> >(key, coord));
               }
            }
         }

         // Create the list of remote indexes
         for(int k=0; k < spRes->cpu(cpuId)->dim(fwdDim)->dim3D(); ++k)
         {
            // Extract "physical" index of third dimension
            k_ = spRes->cpu(cpuId)->dim(fwdDim)->idx3D(k);

            for(int j=0; j < spRes->cpu(cpuId)->dim(fwdDim)->dim2D(k); ++j)
            {
               // Extract "physical" index of second dimension
               j_ = spRes->cpu(cpuId)->dim(fwdDim)->idx2D(j,k);

               for(int i=0; i < spRes->cpu(cpuId)->dim(fwdDim)->dimFwd(j,k); ++i)
               {
                  // Extract "physical" index of first dimension
                  i_ = spRes->cpu(cpuId)->dim(fwdDim)->idxFwd(i,j,k);

                  // Create key as (1D, 3D, 2D)
                  key = std::tr1::make_tuple(i_, j_, k_);

                  // Add key to remote set
                  remoteKeys.insert(key);
               }
            }
         }

         // Setup in the 2D case
      } else if(TFwdB::FieldDimension == Dimensions::TWOD)
      {
         // Create the list of local indexes
         for(int j=0; j < spRes->cpu(FrameworkMacro::id())->dim(fwdDim+1)->dim2D(); ++j)
         {
            // Extract "physical" index of second dimension
            j_ = spRes->cpu(FrameworkMacro::id())->dim(fwdDim+1)->idx2D(j);

            for(int i=0; i < spRes->cpu(FrameworkMacro::id())->dim(fwdDim+1)->dimBwd(j); ++i)
            {
               // Extract "physical" index of first dimension
               i_ = spRes->cpu(FrameworkMacro::id())->dim(fwdDim+1)->idxBwd(i,j);

               // Combine array indexes into coordinate tuple
               coord = std::tr1::make_tuple(i, j, 0);

               // Create key as (1D, 2D, 0)
               key = std::tr1::make_tuple(j_, i_, 0);

               // add key->coordinate to map
               localIdxMap.insert(std::pair<std::tr1::tuple<int,int,int>,std::tr1::tuple<int,int,int> >(key, coord));
            }
         }

         // Create the list of remote indexes
         for(int j=0; j < spRes->cpu(cpuId)->dim(fwdDim)->dim2D(); ++j)
         {
            // Extract "physical" index of second dimension
            j_ = spRes->cpu(cpuId)->dim(fwdDim)->idx2D(j);

            for(int i=0; i < spRes->cpu(cpuId)->dim(fwdDim)->dimFwd(j); ++i)
            {
               // Extract "physical" index of first dimension
               i_ = spRes->cpu(cpuId)->dim(fwdDim)->idxFwd(i,j);

               // Create key as (1D, 2D, 0)
               key = std::tr1::make_tuple(i_, j_, 0);

               // Add key to remote set
               remoteKeys.insert(key);
            }
         }

      // Setup in the 1D case
      } else if(TFwdB::FieldDimension == Dimensions::ONED)
      {
         // Create the list of local indexes
         for(int i=0; i < spRes->cpu(FrameworkMacro::id())->dim(fwdDim+1)->dimBwd(); ++i)
         {
            // Extract "physical" index of first dimension
            i_ = spRes->cpu(FrameworkMacro::id())->dim(fwdDim+1)->idxBwd(i);

            // Combine array indexes into coordinate tuple
            coord = std::tr1::make_tuple(i, 0, 0);

            // Create key as (1D, 0, 0)
            key = std::tr1::make_tuple(i_, 0, 0);

            // add key->coordinate to map
            localIdxMap.insert(std::pair<std::tr1::tuple<int,int,int>,std::tr1::tuple<int,int,int> >(key, coord));
         }

         // Create the list of remote indexes
         for(int i=0; i < spRes->cpu(cpuId)->dim(fwdDim)->dimFwd(); ++i)
         {
            // Extract "physical" index of first dimension
            i_ = spRes->cpu(cpuId)->dim(fwdDim+1)->idxFwd(i);

            // Create key as (1D, 0, 0)
            key = std::tr1::make_tuple(i_, 0, 0);

            // Add key to remote set
            remoteKeys.insert(key);
         }
      }

      // Extract map of shared indexes (stored as keys)
      std::map<std::tr1::tuple<int,int,int>,std::tr1::tuple<int,int,int> >  sharedMap;
      this->extractShared(sharedMap, localIdxMap, remoteKeys);

      // Create the datatype
      this->buildType(data, sharedMap, type);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> void MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::initPositions()
   {
      // Get maximum position size
      int maxSize = std::max(this->nFCpu(), this->nBCpu());

      // Initialise the position values
      for(int i = 0; i < maxSize; ++i)
      {
         this->mRecvPositions.push_back(0);

         this->mSendPositions.push_back(0);
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> void MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::resetRecvPositions()
   {
      // Create position iterator
      std::vector<int>::iterator it;

      // Reset all positions to zero
      for(it = this->mRecvPositions.begin(); it != this->mRecvPositions.end(); ++it)
      {
         (*it) = 0;
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> void MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::resetSendPositions()
   {
      // Create position iterator
      std::vector<int>::iterator it;

      // Reset all positions to zero
      for(it = this->mSendPositions.begin(); it != this->mSendPositions.end(); ++it)
      {
         (*it) = 0;
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> int MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::sendDest(const int id, const int ref, const int size) const
   {
      // Create send ring
      return ((id + 1 + ref) % size);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> int MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::recvSrc(const int id, const int ref, const int size) const
   {
      // Create recv ring
      return ((size - 1 - id + ref) % size);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> void MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::setupRequests()
   {
      // Storage for global location flags
      int dest;
      int src;
      int tag;
      // Storage for CPU group location flags
      int grpMe;
      int grpDest;
      int grpSrc;

      // Storage for the number of packs
      int packs;

      // Initialise forward transform requests
      for(int k = 0; k < this->mForwardPacks.size(); ++k)
      {
         // Get the pack size
         packs = this->mForwardPacks(k);

         // Initialise receive forward with empty requests
         this->mRecvFRequests.insert(std::make_pair<int, std::vector<MPI_Request> >(packs, std::vector<MPI_Request>()));
         for(int id = 0; id < this->nFCpu(); ++id)
         {
            this->mRecvFRequests[packs].push_back(MPI_REQUEST_NULL);
         }

         // Initialise send backward with empty requests
         this->mSendBRequests.insert(std::make_pair<int, std::vector<MPI_Request> >(packs, std::vector<MPI_Request>()));
         for(int id = 0; id < this->nBCpu(); ++id)
         {
            this->mSendBRequests[packs].push_back(MPI_REQUEST_NULL);
         }

         // Create receive forward requests
         for(int id = 0; id < this->nFCpu(); ++id)
         {
            // Get CPU group index of local node
            grpMe = (*std::find(this->mFCpuGroup.begin(), this->mFCpuGroup.end(), FrameworkMacro::id()));

            // Get source index in CPU group
            grpSrc = this->recvSrc(id, grpMe, this->nFCpu());

            // Get global MPI source rank
            src = this->fCpu(grpSrc);

            // Set MPI tag
            tag = src;

            //Safety asserts
            assert(grpSrc < this->mpFBuffers->size());
            assert(grpSrc < this->mFSizes.size());
            assert(grpSrc < this->mRecvFRequests[packs].size());

            // initialise the Recv request
            MPI_Recv_init(this->mpFBuffers->at(grpSrc), packs*this->mFSizes.at(grpSrc), MPI_PACKED, src, tag, MPI_COMM_WORLD, &(this->mRecvFRequests[packs].at(grpSrc)));
         }

         // Create send backward requests
         for(int id = 0; id < this->nBCpu(); ++id)
         {
            // Set MPI tag
            tag = FrameworkMacro::id();
            // Get CPU group index of local node
            grpMe = (*std::find(this->mBCpuGroup.begin(), this->mBCpuGroup.end(), tag));
            // Get destination index in CPU group
            grpDest = this->sendDest(id, grpMe, this->nBCpu());
            // Get global MPI destination rank
            dest = this->bCpu(grpDest);

            //Safety asserts
            assert(grpDest < this->mpBBuffers->size());
            assert(grpDest < this->mBSizes.size());
            assert(grpDest < this->mSendBRequests[packs].size());

            // initialise the Send request
            MPI_Send_init(this->mpBBuffers->at(grpDest), packs*this->mBSizes.at(grpDest), MPI_PACKED, dest, tag, MPI_COMM_WORLD, &(this->mSendBRequests[packs].at(grpDest)));
         }
      }

      // Initialise backward transform requests
      for(int k = 0; k < this->mBackwardPacks.size(); ++k)
      {
         // Get the pack size
         packs = this->mBackwardPacks(k);

         // Initialise receive backward with empty requests
         this->mRecvBRequests.insert(std::make_pair<int, std::vector<MPI_Request> >(packs, std::vector<MPI_Request>()));
         for(int id = 0; id < this->nBCpu(); ++id)
         {
            this->mRecvBRequests[packs].push_back(MPI_REQUEST_NULL);
         }

         // Initialise send forward with empty requests
         this->mSendFRequests.insert(std::make_pair<int, std::vector<MPI_Request> >(packs, std::vector<MPI_Request>()));
         for(int id = 0; id < this->nFCpu(); ++id)
         {
            this->mSendFRequests[packs].push_back(MPI_REQUEST_NULL);
         }

         // Create receive backward requests
         for(int id = 0; id < this->nBCpu(); ++id)
         {
            // Get CPU group index of local node
            grpMe = (*std::find(this->mBCpuGroup.begin(), this->mBCpuGroup.end(), FrameworkMacro::id()));
            // Get source index in CPU group
            grpSrc = this->recvSrc(id, grpMe, this->nBCpu());
            // Get global MPI source rank
            src = this->bCpu(grpSrc);
            // Set MPI tag
            tag = src;
            // initialise the Recv request
            MPI_Recv_init(this->mpBBuffers->at(grpSrc), packs*this->mBSizes.at(grpSrc), MPI_PACKED, src, tag, MPI_COMM_WORLD, &(this->mRecvBRequests[packs].at(grpSrc)));
         }

         // Create send forward requests
         for(int id = 0; id < this->nFCpu(); ++id)
         {
            // Set MPI tag
            tag = FrameworkMacro::id();
            // Get CPU group index of local node
            grpMe = (*std::find(this->mFCpuGroup.begin(), this->mFCpuGroup.end(), tag));
            // Get destination index in CPU group
            grpDest = this->sendDest(id, grpMe, this->nFCpu());
            // Get global MPI destination rank
            dest = this->fCpu(grpDest);
            // initialise the Send request
            MPI_Send_init(this->mpFBuffers->at(grpDest), packs*this->mFSizes.at(grpDest), MPI_PACKED, dest, tag, MPI_COMM_WORLD, &(this->mSendFRequests[packs].at(grpDest)));
         }
      }
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB> void MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB>::cleanupRequests()
   {
      // Create iterator
      std::map<int, std::vector<MPI_Request> >::iterator it;

      // Free requests from Recv B
      for(it = this->mRecvBRequests.begin(); it != this->mRecvBRequests.end(); it++)
      {
         for(unsigned int i = 0; i < (*it).second.size(); ++i)
         {
            MPI_Request_free(&((*it).second.at(i)));
         }
      }

      // Free requests from Recv F
      for(it = this->mRecvFRequests.begin(); it != this->mRecvFRequests.end(); it++)
      {
         for(unsigned int i = 0; i < (*it).second.size(); ++i)
         {
            MPI_Request_free(&((*it).second.at(i)));
         }
      }

      // Free requests from Send B
      for(it = this->mSendBRequests.begin(); it != this->mSendBRequests.end(); it++)
      {
         for(unsigned int i = 0; i < (*it).second.size(); ++i)
         {
            MPI_Request_free(&((*it).second.at(i)));
         }
      }

      // Free requests from Send F
      for(it = this->mSendFRequests.begin(); it != this->mSendFRequests.end(); it++)
      {
         for(unsigned int i = 0; i < (*it).second.size(); ++i)
         {
            MPI_Request_free(&((*it).second.at(i)));
         }
      }
   }

}
}

#endif // MPICONVERTERBASE_HPP
