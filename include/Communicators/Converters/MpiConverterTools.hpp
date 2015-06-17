/**
 * @file MpiConverterTools.hpp
 * @brief Implementation of some tools used by the MPI converter 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef MPICONVERTERTOOLS_HPP
#define MPICONVERTERTOOLS_HPP

// Debug includes
//

// Configuration includes
//
#include "Framework/FrameworkMacro.h"

// System includes
//
#include <cassert>
#include <set>
#include <map>
#include <tr1/tuple>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MpiTypes.hpp"
#include "Enums/Dimensions.hpp"
#include "Enums/DimensionTools.hpp"
#include "Resolutions/Resolution.hpp"

namespace GeoMHDiSCC {

namespace Parallel {

   /**
    * @brief Implementation of some tools used by the MPI converter.
    */
   template <Dimensions::Type TID> class MpiConverterTools;

   /**
    * @brief Specialisation of the tools used by the MPI converter for 3D simulation.
    */
   template <> class MpiConverterTools<Dimensions::THREED>
   {
      public:
         /// Typedef for a three indexes specifing a point
         typedef std::tr1::tuple<int,int,int>   Coordinate;

         static void buildLocalFwdMap(std::map<Coordinate,Coordinate>& rLocalIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim);

         template <typename TIdx> static void buildLocalBwdMap(std::map<Coordinate,Coordinate>& rLocalIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, SharedPtrMacro<TIdx> spIdxConv);

         /**
          * @brief Build a forward MPI Datatype
          *
          * @param spRes   Shared Resolution
          * @param fwdDim  Dimension index for forward transform
          * @param data    Input data
          * @param cpuId   ID of the CPU
          */
         template <typename TData> static MPI_Datatype buildFwdDatatype(const std::map<Coordinate,Coordinate>& localIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, typename Datatypes::FlatScalarField<TData,Dimensions::THREED> &data, const int cpuId);

         /**
          * @brief Build a backward MPI Datatype
          *
          * @param spRes   Shared Resolution
          * @param fwdDim  Dimension index for forward transform
          * @param data    Input data
          * @param cpuId   ID of the CPU
          */
         template <typename TData> static MPI_Datatype buildBwdDatatype(const std::map<Coordinate,Coordinate>& localIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, typename Datatypes::FlatScalarField<TData,Dimensions::THREED> &data, const int cpuId);

         /**
          * @brief Extract shared indexes
          *
          * @param sharedMap     Storage for the shared index map
          * @param localIdxMap   Local node key to indexes map 
          * @param remoteKeys    Remote node index keys
          */
         static void extractShared(std::map<Coordinate,Coordinate>& sharedMap, const std::map<Coordinate,Coordinate>& localIdxMap, const std::set<Coordinate>& remoteKeys);

         /**
          * @brief Create type
          *
          * @param data       The concerned data
          * @param sharedMap  Shared index map
          * @param type       MPI datatype storage
          */
         template <typename TData> static MPI_Datatype buildType(typename Datatypes::FlatScalarField<TData,Dimensions::THREED> &data, const std::map<Coordinate,Coordinate>& sharedMap);
         
      protected:

      private:
         /**
          * @brief Constructor
          */
         MpiConverterTools();

         /**
          * @brief Destructor
          */
         ~MpiConverterTools();

   };

   template <typename TData> MPI_Datatype MpiConverterTools<Dimensions::THREED>::buildType(typename Datatypes::FlatScalarField<TData,Dimensions::THREED> &data, const std::map<MpiConverterTools<Dimensions::THREED>::Coordinate,MpiConverterTools<Dimensions::THREED>::Coordinate>& sharedMap)
   {
      // Get the number of elements
      int nElements = sharedMap.size();

      // Prepare data required for MPI datatype
      MPI_Aint    displ[nElements];
      int         blocks[nElements];

      // Prepare data required for MPI datatype
      MPI_Aint    base;
      MPI_Aint    element;

      // Set the base of the datatype to the start of the data matrix
      MPI_Get_address(data.rData().data(),&element);
      base = element;

      // Create map iterator
      std::map<Coordinate,Coordinate>::const_iterator it;

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
      MPI_Datatype   type;
      MPI_Type_create_hindexed(nElements, blocks, displ, MpiTypes::type<TData>(), &type);

      // Commit MPI datatype
      MPI_Type_commit(&type);

      return type;
   }

   inline void MpiConverterTools<Dimensions::THREED>::buildLocalFwdMap(std::map<Coordinate,Coordinate>&  rLocalIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim)
   {
      // Make sure key map is empty
      rLocalIdxMap.clear();

      // Storage for the generated key
      Coordinate key;
      // Storage for the coordinate indexes
      Coordinate coord;
      
      // Storage for the simulation wide indexes
      int i_, j_, k_;

      //
      // Create the list of local indexes
      //

      // Loop over slow data dimension
      for(int k = 0; k < spRes->cpu()->dim(fwdDim)->dim<Dimensions::Data::DAT3D>(); ++k)
      {
         // Extract "physical" index of slow data dimension
         k_ = spRes->cpu()->dim(fwdDim)->idx<Dimensions::Data::DAT3D>(k);

         // Loop over middle data dimension
         for(int j = 0; j < spRes->cpu()->dim(fwdDim)->dim<Dimensions::Data::DAT2D>(k); ++j)
         {
            // Extract "physical" index of middle data dimension
            j_ = spRes->cpu()->dim(fwdDim)->idx<Dimensions::Data::DAT2D>(j,k);

            // Loop over forward data dimension
            for(int i = 0; i < spRes->cpu()->dim(fwdDim)->dim<Dimensions::Data::DATF1D>(k); ++i)
            {
               // Extract "physical" index of forward data dimension
               i_ = spRes->cpu()->dim(fwdDim)->idx<Dimensions::Data::DATF1D>(i,k);

               // Combine array indexes into coordinate tuple
               coord = std::tr1::make_tuple(i, j, k);

               // Create key as (1D, 2D, 3D)
               key = spRes->counter()->makeKey(fwdDim, i_, j_, k_);

               // add key->coordinate to map
               rLocalIdxMap.insert(std::pair<Coordinate,Coordinate>(key, coord));
            }
         }
      }
   }

   template <typename TIdx> void MpiConverterTools<Dimensions::THREED>::buildLocalBwdMap(std::map<Coordinate,Coordinate>& rLocalIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, SharedPtrMacro<TIdx> spIdxConv)
   {
      // Make sure key map is empty
      rLocalIdxMap.clear();

      // Storage for the index
      Coordinate coord;
      
      // Storage for the simulation wide indexes
      int i_, j_, k_;

      // Storage for the generated key
      Coordinate key;

      //
      // Create the list of local indexes
      //

      // Loop over slow data dimension
      for(int k = 0; k < spRes->cpu()->dim(Dimensions::jump(fwdDim,1))->dim<Dimensions::Data::DAT3D>(); ++k)
      {
         // Extract "physical" index of slow data dimension
         k_ = spRes->cpu()->dim(Dimensions::jump(fwdDim,1))->idx<Dimensions::Data::DAT3D>(k);

         // Loop over middle data dimension
         for(int j = 0; j < spRes->cpu()->dim(Dimensions::jump(fwdDim,1))->dim<Dimensions::Data::DAT2D>(k); ++j)
         {
            // Extract "physical" index of middle data dimension
            j_ = spRes->cpu()->dim(Dimensions::jump(fwdDim,1))->idx<Dimensions::Data::DAT2D>(j,k);

            // Loop over backward data dimension
            for(int i = 0; i < spRes->cpu()->dim(Dimensions::jump(fwdDim,1))->dim<Dimensions::Data::DATB1D>(k); ++i)
            {
               // Extract "physical" index of backward data dimension
               i_ = spRes->cpu()->dim(Dimensions::jump(fwdDim,1))->idx<Dimensions::Data::DATB1D>(i,k);

               // Combine array indexes into coordinate tuple
               coord = std::tr1::make_tuple(i + spIdxConv->centralPadding(i_, k), j, k);

               // Create key as (2D, 3D, 1D)
               key = spRes->counter()->makeKey(Dimensions::jump(fwdDim,1), i_, j_, k_);

               // add key->coordinate to map
               rLocalIdxMap.insert(std::pair<Coordinate,Coordinate>(key, coord));
            }
         }
      }
   }


   template <typename TData> MPI_Datatype MpiConverterTools<Dimensions::THREED>::buildFwdDatatype(const std::map<Coordinate,Coordinate>&  localIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, typename Datatypes::FlatScalarField<TData,Dimensions::THREED> &data, const int cpuId)
   {
      // List of remote keys
      std::set<Coordinate>  remoteKeys;

      // Storage for the simulation wide indexes
      int i_, j_, k_;

      // Storage for the generated key
      Coordinate key;

      //
      // Create the list of remote indexes in next transform
      //

      // Remote is also local CPU
      MatrixI  matRemote;
      if(cpuId == FrameworkMacro::transformId(fwdDim))
      {
         // Loop over slow data dimension
         for(int k=0; k < spRes->cpu()->dim(Dimensions::jump(fwdDim,1))->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            // Extract "physical" index of slow data dimension
            k_ = spRes->cpu()->dim(Dimensions::jump(fwdDim,1))->idx<Dimensions::Data::DAT3D>(k);

            // Loop over middle data dimension
            for(int j=0; j < spRes->cpu()->dim(Dimensions::jump(fwdDim,1))->dim<Dimensions::Data::DAT2D>(k); ++j)
            {
               // Extract "physical" index of middle data dimension
               j_ = spRes->cpu()->dim(Dimensions::jump(fwdDim,1))->idx<Dimensions::Data::DAT2D>(j,k);

               // Loop over backward data dimension
               for(int i=0; i < spRes->cpu()->dim(Dimensions::jump(fwdDim,1))->dim<Dimensions::Data::DATB1D>(k); ++i)
               {
                  // Extract "physical" index of backward data dimension
                  i_ = spRes->cpu()->dim(Dimensions::jump(fwdDim,1))->idx<Dimensions::Data::DATB1D>(i,k);

                  // Create key as (2D, 3D, 1D) indexes (i.e. data gets transposed during communication)
                  key = spRes->counter()->makeKey(Dimensions::jump(fwdDim,1), i_, j_, k_);

                  // Add key to remote set
                  remoteKeys.insert(key);
               }
            }
         }

         // Convert remote keys to matrix to send througg MPI
         matRemote.resize(3, remoteKeys.size());
         int i = 0; 
         for(std::set<Coordinate>::iterator it = remoteKeys.begin(); it != remoteKeys.end(); ++it)
         {
            matRemote(0,i) = std::tr1::get<0>(*it);
            matRemote(1,i) = std::tr1::get<1>(*it);
            matRemote(2,i) = std::tr1::get<2>(*it);
            i++;
         }
   
         // Broadcast size
         i = remoteKeys.size();
         MPI_Bcast(&i, 1, MPI_INT, cpuId, FrameworkMacro::transformComm(fwdDim));

         // Broadcast data
         MPI_Bcast(matRemote.data(), matRemote.cols()*matRemote.rows(), MPI_INT, cpuId, FrameworkMacro::transformComm(fwdDim)); 

      // Remote CPU needs to generate list 
      } else
      {
         // Get size
         int nCoords;
         MPI_Bcast(&nCoords, 1, MPI_INT, cpuId, FrameworkMacro::transformComm(fwdDim));

         // Get remote keys as matrix
         matRemote.resize(3, nCoords);
         MPI_Bcast(matRemote.data(), matRemote.cols()*matRemote.rows(), MPI_INT, cpuId, FrameworkMacro::transformComm(fwdDim)); 

         // Convert matrix to remoteKeys set
         for(int i = 0; i < nCoords; i++)
         {
            key = std::tr1::make_tuple(matRemote(0,i), matRemote(1,i), matRemote(2,i));
            remoteKeys.insert(key);
         }
      }

      // Synchronize
      FrameworkMacro::syncTransform(fwdDim);

      // Extract map of shared indexes (stored as keys)
      std::map<Coordinate,Coordinate>  sharedMap;
      MpiConverterTools<Dimensions::THREED>::extractShared(sharedMap, localIdxMap, remoteKeys);

      // Create the datatype
      MPI_Datatype type = MpiConverterTools<Dimensions::THREED>::buildType(data, sharedMap);

      return type;
   }

   template <typename TData> MPI_Datatype MpiConverterTools<Dimensions::THREED>::buildBwdDatatype(const std::map<Coordinate,Coordinate>&  localIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, typename Datatypes::FlatScalarField<TData,Dimensions::THREED> &data, const int cpuId)
   {
      // List of remote keys
      std::set<Coordinate>  remoteKeys;
      
      // Storage for the simulation wide indexes
      int i_, j_, k_;

      // Storage for the generated key
      Coordinate key;

      //
      // Create the list of remote indexes
      //

      // Remote is also local CPU
      MatrixI  matRemote;
      if(cpuId == FrameworkMacro::transformId(fwdDim))
      {
         // Loop over slow data dimension
         for(int k=0; k < spRes->cpu()->dim(fwdDim)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            // Extract "physical" index of slow data dimension
            k_ = spRes->cpu()->dim(fwdDim)->idx<Dimensions::Data::DAT3D>(k);

            // Loop over middle data dimension
            for(int j=0; j < spRes->cpu()->dim(fwdDim)->dim<Dimensions::Data::DAT2D>(k); ++j)
            {
               // Extract "physical" index of middle data dimension
               j_ = spRes->cpu()->dim(fwdDim)->idx<Dimensions::Data::DAT2D>(j,k);

               // Loop over forward data dimension
               for(int i=0; i < spRes->cpu()->dim(fwdDim)->dim<Dimensions::Data::DATF1D>(k); ++i)
               {
                  // Extract "physical" index of forward data dimension
                  i_ = spRes->cpu()->dim(fwdDim)->idx<Dimensions::Data::DATF1D>(i,k);

                  // Create key as (1D, 2D, 3D)
                  key = spRes->counter()->makeKey(fwdDim, i_, j_, k_);

                  // Add key to remote set
                  remoteKeys.insert(key);
               }
            }
         }

         // Convert remote keys to matrix to send through MPI
         matRemote.resize(3, remoteKeys.size());
         int i = 0; 
         for(std::set<Coordinate>::iterator it = remoteKeys.begin(); it != remoteKeys.end(); ++it)
         {
            matRemote(0,i) = std::tr1::get<0>(*it);
            matRemote(1,i) = std::tr1::get<1>(*it);
            matRemote(2,i) = std::tr1::get<2>(*it);
            i++;
         }

         // Broadcast size
         i = remoteKeys.size();
         MPI_Bcast(&i, 1, MPI_INT, cpuId, FrameworkMacro::transformComm(fwdDim));

         // Broadcast data
         MPI_Bcast(matRemote.data(), matRemote.rows()*matRemote.cols(), MPI_INT, cpuId, FrameworkMacro::transformComm(fwdDim)); 

      // Remote CPU needs to generate list 
      } else
      {
         // Get size
         int nCoords;
         MPI_Bcast(&nCoords, 1, MPI_INT, cpuId, FrameworkMacro::transformComm(fwdDim));

         // Get remot ekeys as matrix
         matRemote.resize(3, nCoords);
         MPI_Bcast(matRemote.data(), matRemote.cols()*matRemote.rows(), MPI_INT, cpuId, FrameworkMacro::transformComm(fwdDim)); 

         // Convert matrix to remoteKeys set
         for(int i = 0; i < nCoords; i++)
         {
            key = std::tr1::make_tuple(matRemote(0,i), matRemote(1,i), matRemote(2,i));
            remoteKeys.insert(key);
         }
      }

      // Synchronize
      FrameworkMacro::syncTransform(fwdDim);

      // Extract map of shared indexes (stored as keys)
      std::map<Coordinate,Coordinate>  sharedMap;
      MpiConverterTools<Dimensions::THREED>::extractShared(sharedMap, localIdxMap, remoteKeys);

      // Create the datatype
      MPI_Datatype type = MpiConverterTools<Dimensions::THREED>::buildType(data, sharedMap);

      return type;
   }

   /**
    * @brief Specialisation of the tools used by the MPI converter for 2D simulation.
    */
   template <> class MpiConverterTools<Dimensions::TWOD>
   {
      public:
         /// Typedef for a three indexes specifing a point
         typedef std::pair<int,int>   Coordinate;

         /**
          * @brief Build a forward MPI Datatype
          *
          * @param spRes   Shared Resolution
          * @param fwdDim  Dimension index for forward transform
          * @param data    Input data
          * @param type    Created MPI data type
          * @param cpuId   ID of the CPU
          */
         template <typename TData> static MPI_Datatype buildFwdDatatype(SharedResolution spRes, const Dimensions::Transform::Id fwdDim, typename Datatypes::FlatScalarField<TData,Dimensions::TWOD> &data, const int cpuId);

         /**
          * @brief Build a backward MPI Datatype
          *
          * @param spRes   Shared Resolution
          * @param fwdDim  Dimension index for forward transform
          * @param data    Input data
          * @param type    Created MPI data type
          * @param cpuId   ID of the CPU
          */
         template <typename TData> static MPI_Datatype buildBwdDatatype(SharedResolution spRes, const Dimensions::Transform::Id fwdDim, typename Datatypes::FlatScalarField<TData,Dimensions::TWOD> &data, const int cpuId);

         /**
          * @brief Extract shared indexes
          *
          * @param sharedMap     Storage for the shared index map
          * @param localIdxMap   Local node key to indexes map 
          * @param remoteKeys    Remote node index keys
          */
         static void extractShared(std::map<Coordinate,Coordinate>& sharedMap, const std::map<Coordinate,Coordinate>& localIdxMap, const std::set<Coordinate>& remoteKeys);

         /**
          * @brief Create type
          *
          * @param data       The concerned data
          * @param sharedMap  Shared index map
          * @param type       MPI datatype storage
          */
         template <typename TData> static MPI_Datatype buildType(typename Datatypes::FlatScalarField<TData,Dimensions::TWOD> &data, const std::map<Coordinate,Coordinate>& sharedMap);
         
      protected:

      private:
         /**
          * @brief Constructor
          */
         MpiConverterTools();

         /**
          * @brief Destructor
          */
         ~MpiConverterTools();

   };

   template <typename TData> MPI_Datatype MpiConverterTools<Dimensions::TWOD>::buildType(typename Datatypes::FlatScalarField<TData,Dimensions::TWOD> &data, const std::map<MpiConverterTools<Dimensions::TWOD>::Coordinate,MpiConverterTools<Dimensions::TWOD>::Coordinate>& sharedMap)
   {
      // Get the number of elements
      int nElements = sharedMap.size();

      // Prepare data required for MPI datatype
      MPI_Aint    displ[nElements];
      int         blocks[nElements];

      // Prepare data required for MPI datatype
      MPI_Aint    base;
      MPI_Aint    element;

      // Set the base of the datatype to the start of the data matrix
      MPI_Get_address(data.rData().data(),&element);
      base = element;

      // Create iterator
      std::map<Coordinate,Coordinate>::const_iterator it;

      // Create MPI displacement list
      int tot = 0;
      for(it = sharedMap.begin(); it != sharedMap.end(); ++it)
      {
         // Get address of stored coordinates
         MPI_Get_address(&data.rPoint(it->second.first, it->second.second), &element);

         // Fill datatype information
         displ[tot] = element - base;
         blocks[tot] = 1;

         // Increment datatype size
         tot++;
      }

      // Create MPI datatype
      MPI_Datatype  type;
      MPI_Type_create_hindexed(nElements, blocks, displ, MpiTypes::type<TData>(), &type);

      // Commit MPI datatype
      MPI_Type_commit(&type);

      return type;
   }

   template <typename TData> MPI_Datatype MpiConverterTools<Dimensions::TWOD>::buildFwdDatatype(SharedResolution spRes, const Dimensions::Transform::Id fwdDim, typename Datatypes::FlatScalarField<TData,Dimensions::TWOD> &data, const int cpuId)
   {
      // Create  map of the local indexes to unique keys
      std::map<Coordinate,Coordinate>  localIdxMap;
      // List of remote keys
      std::set<Coordinate>  remoteKeys;

      // Storage for the index
      Coordinate coord;
      
      // Storage for the simulation wide indexes
      int i_, j_;

      // Storage for the generated key
      Coordinate key;

      //
      // Create the list of local indexes
      //

      // Loop over slow data dimension
      for(int j=0; j < spRes->cpu(FrameworkMacro::id())->dim(fwdDim)->dim<Dimensions::Data::DAT2D>(); ++j)
      {
         // Extract "physical" index of slow data dimension
         j_ = spRes->cpu(FrameworkMacro::id())->dim(fwdDim)->idx<Dimensions::Data::DAT2D>(j);

         // Loop over forward data dimension
         for(int i=0; i < spRes->cpu(FrameworkMacro::id())->dim(fwdDim)->dim<Dimensions::Data::DATF1D>(j); ++i)
         {
            // Extract "physical" index of forward data dimension
            i_ = spRes->cpu(FrameworkMacro::id())->dim(fwdDim)->idx<Dimensions::Data::DATF1D>(i,j);

            // Combine array indexes into coordinate
            coord = std::make_pair(i, j);

            // Create key as (1D, 2D)
            key = std::make_pair(i_, j_);

            // add key->coordinate to map
            localIdxMap.insert(std::pair<Coordinate,Coordinate>(key, coord));
         }
      }

      //
      // Create the list of remote indexes
      //

      // Loop over slow data dimension
      for(int j=0; j < spRes->cpu(cpuId)->dim(fwdDim)->dim<Dimensions::Data::DAT2D>(); ++j)
      {
         // Extract "physical" index of slow data dimension
         j_ = spRes->cpu(cpuId)->dim(fwdDim)->idx<Dimensions::Data::DAT2D>(j);

         // Loop over backward data dimension
         for(int i=0; i < spRes->cpu(cpuId)->dim(Dimensions::jump(fwdDim,1))->dim<Dimensions::Data::DATB1D>(j); ++i)
         {
            // Extract "physical" index of backward data dimension
            i_ = spRes->cpu(cpuId)->dim(Dimensions::jump(fwdDim,1))->idx<Dimensions::Data::DATB1D>(i,j);

            // Create key as (2D, 1D)
            key = std::make_pair(j_, i_);

            // Add key to remote set
            remoteKeys.insert(key);
         }
      }

      // Extract map of shared indexes (stored as keys)
      std::map<Coordinate,Coordinate>  sharedMap;
      MpiConverterTools<Dimensions::TWOD>::extractShared(sharedMap, localIdxMap, remoteKeys);

      // Create the datatype
      MPI_Datatype type = MpiConverterTools<Dimensions::TWOD>::buildType(data, sharedMap);

      return type;
   }

   template <typename TData> MPI_Datatype MpiConverterTools<Dimensions::TWOD>::buildBwdDatatype(SharedResolution spRes, const Dimensions::Transform::Id fwdDim, typename Datatypes::FlatScalarField<TData,Dimensions::TWOD> &data, const int cpuId)
   {
      // Create  map of the local indexes to unique keys
      std::map<Coordinate,Coordinate>  localIdxMap;

      // List of remote keys
      std::set<Coordinate>  remoteKeys;

      // Storage for the index
      Coordinate coord;
      
      // Storage for the simulation wide indexes
      int i_, j_;

      // Storage for the generated key
      Coordinate key;

      //
      // Create the list of local indexes
      //

      // Loop over slow data dimension
      for(int j=0; j < spRes->cpu(FrameworkMacro::id())->dim(Dimensions::jump(fwdDim,1))->dim<Dimensions::Data::DAT2D>(); ++j)
      {
         // Extract "physical" index of slow data dimension
         j_ = spRes->cpu(FrameworkMacro::id())->dim(Dimensions::jump(fwdDim,1))->idx<Dimensions::Data::DAT2D>(j);

         // Loop over backward data dimension
         for(int i=0; i < spRes->cpu(FrameworkMacro::id())->dim(Dimensions::jump(fwdDim,1))->dim<Dimensions::Data::DATB1D>(j); ++i)
         {
            // Extract "physical" index of backward data dimension
            i_ = spRes->cpu(FrameworkMacro::id())->dim(Dimensions::jump(fwdDim,1))->idx<Dimensions::Data::DATB1D>(i,j);

            // Combine array indexes into coordinate pair
            coord = std::make_pair(i, j);

            // Create key as (2D, 1D)
            key = std::make_pair(j_, i_);

            // add key->coordinate to map
            localIdxMap.insert(std::pair<Coordinate,Coordinate>(key, coord));
         }
      }

      //
      // Create the list of remote indexes
      //

      // Loop over slow data dimension
      for(int j=0; j < spRes->cpu(cpuId)->dim(fwdDim)->dim<Dimensions::Data::DAT2D>(); ++j)
      {
         // Extract "physical" index of slow data dimension
         j_ = spRes->cpu(cpuId)->dim(fwdDim)->idx<Dimensions::Data::DAT2D>(j);

         // Loop ver forward data dimension
         for(int i=0; i < spRes->cpu(cpuId)->dim(fwdDim)->dim<Dimensions::Data::DATF1D>(j); ++i)
         {
            // Extract "physical" index of forward data dimension
            i_ = spRes->cpu(cpuId)->dim(fwdDim)->idx<Dimensions::Data::DATF1D>(i,j);

            // Create key as (1D, 2D)
            key = std::make_pair(i_, j_);

            // Add key to remote set
            remoteKeys.insert(key);
         }
      }

      // Extract map of shared indexes (stored as keys)
      std::map<Coordinate,Coordinate>  sharedMap;
      MpiConverterTools<Dimensions::TWOD>::extractShared(sharedMap, localIdxMap, remoteKeys);

      // Create the datatype
      MPI_Datatype type = MpiConverterTools<Dimensions::TWOD>::buildType(data, sharedMap, type);

      return type;
   }

   /**
    * @brief Specialisation of the tools used by the MPI converter for 1D simulation.
    */
   template <> class MpiConverterTools<Dimensions::ONED>
   {
      public:
         /// Typedef for a three indexes specifing a point
         typedef int Coordinate;

         /**
          * @brief Build a forward MPI Datatype
          *
          * @param spRes   Shared Resolution
          * @param fwdDim  Dimension index for forward transform
          * @param data    Input data
          * @param type    Created MPI data type
          * @param cpuId   ID of the CPU
          */
         template <typename TData> static MPI_Datatype buildFwdDatatype(SharedResolution spRes, const Dimensions::Transform::Id fwdDim, typename Datatypes::FlatScalarField<TData,Dimensions::ONED> &data, const int cpuId);

         /**
          * @brief Build a backward MPI Datatype
          *
          * @param spRes   Shared Resolution
          * @param fwdDim  Dimension index for forward transform
          * @param data    Input data
          * @param type    Created MPI data type
          * @param cpuId   ID of the CPU
          */
         template <typename TData> static MPI_Datatype buildBwdDatatype(SharedResolution spRes, const Dimensions::Transform::Id fwdDim, typename Datatypes::FlatScalarField<TData,Dimensions::ONED> &data, const int cpuId);

         /**
          * @brief Extract shared indexes
          *
          * @param sharedMap     Storage for the shared index map
          * @param localIdxMap   Local node key to indexes map 
          * @param remoteKeys    Remote node index keys
          */
         static void extractShared(std::map<Coordinate,Coordinate>& sharedMap, const std::map<Coordinate,Coordinate>& localIdxMap, const std::set<Coordinate>& remoteKeys);

         /**
          * @brief Create type
          *
          * @param data       The concerned data
          * @param sharedMap  Shared index map
          * @param type       MPI datatype storage
          */
         template <typename TData> static MPI_Datatype buildType(typename Datatypes::FlatScalarField<TData,Dimensions::ONED> &data, const std::map<Coordinate,Coordinate>& sharedMap);
         
      protected:

      private:
         /**
          * @brief Constructor
          */
         MpiConverterTools();

         /**
          * @brief Destructor
          */
         ~MpiConverterTools();

   };

   template <typename TData> MPI_Datatype MpiConverterTools<Dimensions::ONED>::buildType(typename Datatypes::FlatScalarField<TData,Dimensions::ONED> &data, const std::map<MpiConverterTools<Dimensions::ONED>::Coordinate,MpiConverterTools<Dimensions::ONED>::Coordinate>& sharedMap)
   {
      // Get the number of elements
      int nElements = sharedMap.size();

      // Prepare data required for MPI datatype
      MPI_Aint    displ[nElements];
      int         blocks[nElements];

      // Prepare data required for MPI datatype
      MPI_Aint    base;
      MPI_Aint    element;

      // Set the base of the datatype to the start of the data matrix
      MPI_Get_address(data.rData().data(),&element);
      base = element;

      // Create iterator
      std::map<Coordinate,Coordinate>::const_iterator it;

      // Create MPI displacement list
      int tot = 0;
      for(it = sharedMap.begin(); it != sharedMap.end(); ++it)
      {
         // Get address of stored coordinates
         MPI_Get_address(&data.rPoint(it->second), &element);

         // Fill datatype information
         displ[tot] = element - base;
         blocks[tot] = 1;

         // Increment datatype size
         tot++;
      }

      // Create MPI datatype
      MPI_Datatype type;
      MPI_Type_create_hindexed(nElements, blocks, displ, MpiTypes::type<TData>(), &type);

      // Commit MPI datatype
      MPI_Type_commit(&type);

      return type;
   }

   template <typename TData> MPI_Datatype MpiConverterTools<Dimensions::ONED>::buildFwdDatatype(SharedResolution spRes, const Dimensions::Transform::Id fwdDim, typename Datatypes::FlatScalarField<TData,Dimensions::ONED> &data, const int cpuId)
   {
      // Create  map of the local indexes to unique keys
      std::map<Coordinate,Coordinate>  localIdxMap;

      // List of remote keys
      std::set<Coordinate>  remoteKeys;

      // Storage for the index
      Coordinate coord;

      // Storage for the simulation wide indexes
      int i_;

      // Storage for the generated key
      Coordinate key;

      //
      // Create the list of local indexes
      //

      // Loop over forward data dimension
      for(int i=0; i < spRes->cpu(FrameworkMacro::id())->dim(fwdDim)->dim<Dimensions::Data::DATF1D>(); ++i)
      {
         // Extract "physical" index of forward data dimension
         i_ = spRes->cpu(FrameworkMacro::id())->dim(fwdDim)->idx<Dimensions::Data::DATF1D>(i);

         // Set coordindate
         coord = i;

         // Set key
         key = i_;

         // add key->coordinate to map
         localIdxMap.insert(std::pair<Coordinate,Coordinate>(key, coord));
      }

      //
      // Create the list of remote indexes
      //

      // Loop over backward data dimension
      for(int i=0; i < spRes->cpu(cpuId)->dim(Dimensions::jump(fwdDim,1))->dim<Dimensions::Data::DATB1D>(); ++i)
      {
         // Extract "physical" index of backward data dimension
         i_ = spRes->cpu(cpuId)->dim(Dimensions::jump(fwdDim,1))->idx<Dimensions::Data::DATB1D>(i);

         // Set key
         key = i_;

         // Add key to remote set
         remoteKeys.insert(key);
      }

      // Extract map of shared indexes (stored as keys)
      std::map<Coordinate,Coordinate>  sharedMap;
      MpiConverterTools<Dimensions::ONED>::extractShared(sharedMap, localIdxMap, remoteKeys);

      // Create the datatype
      MPI_Datatype type = MpiConverterTools<Dimensions::ONED>::buildType(data, sharedMap, type);

      return type;
   }

   template <typename TData> MPI_Datatype MpiConverterTools<Dimensions::ONED>::buildBwdDatatype(SharedResolution spRes, const Dimensions::Transform::Id fwdDim, typename Datatypes::FlatScalarField<TData,Dimensions::ONED> &data, const int cpuId)
   {
      // Create  map of the local indexes to unique keys
      std::map<Coordinate,Coordinate>  localIdxMap;
      // List of remote keys
      std::set<Coordinate>  remoteKeys;

      // Storage for the index
      Coordinate coord;
      
      // Storage for the simulation wide indexes
      int i_;

      // Storage for the generated key
      Coordinate key;

      //
      // Create the list of local indexes
      //

      // Loop ver backward data dimension
      for(int i=0; i < spRes->cpu(FrameworkMacro::id())->dim(Dimensions::jump(fwdDim,1))->dim<Dimensions::Data::DATB1D>(); ++i)
      {
         // Extract "physical" index of backward data dimension
         i_ = spRes->cpu(FrameworkMacro::id())->dim(Dimensions::jump(fwdDim,1))->idx<Dimensions::Data::DATB1D>(i);

         // Set coordinate 
         coord = i;

         // Set key
         key = i_;

         // add key->coordinate to map
         localIdxMap.insert(std::pair<Coordinate,Coordinate>(key, coord));
      }

      //
      // Create the list of remote indexes
      //

      // Loop over forward data dimension
      for(int i=0; i < spRes->cpu(cpuId)->dim(fwdDim)->dim<Dimensions::Data::DATF1D>(); ++i)
      {
         // Extract "physical" index of forward data dimension
         i_ = spRes->cpu(cpuId)->dim(Dimensions::jump(fwdDim,1))->idx<Dimensions::Data::DATF1D>(i);

         // Set key
         key = i_;

         // Add key to remote set
         remoteKeys.insert(key);
      }

      // Extract map of shared indexes (stored as keys)
      std::map<Coordinate,Coordinate>  sharedMap;
      MpiConverterTools<Dimensions::ONED>::extractShared(sharedMap, localIdxMap, remoteKeys);

      // Create the datatype
      MPI_Datatype type = MpiConverterTools<Dimensions::ONED>::buildType(data, sharedMap, type);

      return type;
   }

}
}

#endif // MPICONVERTERTOOLS_HPP
