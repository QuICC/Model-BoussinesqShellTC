/** \file MpiConverterTools.hpp
 *  \brief Implementation of some tools used by the MPI converter.
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

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MpiTypes.hpp"
#include "Enums/Dimensions.hpp"
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

         /**
          * @brief Build a forward MPI Datatype
          *
          * @param spRes   Shared Resolution
          * @param fwdDim  Dimension index for forward transform
          * @param data    Input data
          * @param type    Created MPI data type
          * @param cpuId   ID of the CPU
          */
         template <typename TFwd> void buildFwdDatatype(SharedResolution spRes, const Dimensions::Transform::Id fwdDim, TFwd &data, MPI_Datatype &type, const int cpuId);

         /**
          * @brief Build a backward MPI Datatype
          *
          * @param spRes   Shared Resolution
          * @param fwdDim  Dimension index for forward transform
          * @param data    Input data
          * @param type    Created MPI data type
          * @param cpuId   ID of the CPU
          */
         template <typename TBwd> void buildBwdDatatype(SharedResolution spRes, const Dimensions::Transform::Id fwdDim, TBwd &data, MPI_Datatype &type, const int cpuId);

         /**
          * @brief Extract shared indexes
          *
          * @param sharedMap     Storage for the shared index map
          * @param localIdxMap   Local node key to indexes map 
          * @param remoteKeys    Remote node index keys
          */
         void extractShared(std::map<Coordinate,Coordinate>& sharedMap, const std::map<Coordinate,Coordinate>& localIdxMap, const std::set<Coordinate>& remoteKeys);

         /**
          * @brief Create type
          *
          * @param data       The concerned data
          * @param sharedMap  Shared index map
          * @param type       MPI datatype storage
          */
         template <typename TData> void buildType(TData &data, const std::map<Coordinate,Coordinate>& sharedMap, MPI_Datatype &type);
         
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

   template <> void MpiConverterTools<Dimensions::THREED>::extractShared(std::map<MpiConverterToos::<Dimensions::THREED>::Coordinate,MpiConverterToos::<Dimensions::THREED>::Coordinate>& sharedMap, const std::map<MpiConverterToos::<Dimensions::THREED>::Coordinate,MpiConverterToos::<Dimensions::THREED>::Coordinate>& localIdxMap, const std::set<MpiConverterToos::<Dimensions::THREED>::Coordinate>& remoteKeys)
   {
      // List of local index keys 
      std::set<Coordinate>  localKeys;

      // Create map iterator
      std::map<Coordinate,Coordinate>::const_iterator mapIt;
      
      // Extract the set of local keys
      for(mapIt = localIdxMap.begin(); mapIt != localIdxMap.end(); ++mapIt)
      {
         localKeys.insert(mapIt->first);
      }
      
      // Storage for the shared keys
      std::set<Coordinate> sharedKeys;

      // Create the list of common indexes
      std::set_intersection(localKeys.begin(), localKeys.end(), remoteKeys.begin(), remoteKeys.end(), std::inserter(sharedKeys, sharedKeys.begin()));

      // Clear the shared map
      sharedMap.clear();

      // Fill shared map
      std::set<Coordinate>::iterator sit;
      for(sit = sharedKeys.begin(); sit != sharedKeys.end(); sit++)
      {
         sharedMap.insert(std::make_pair(*sit, localIdxMap.find(*sit)->second));
      }
   }

   template <> template <typename TData> void MpiConverterTools<Dimensions::THREED>::buildType(TData &data, const std::map<MpiConverterTools<Dimensions::THREED>::Coordinate,MpiConverterTools<Dimensions::THREED>::Coordinate>& sharedMap, MPI_Datatype &type)
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
      MPI_Type_create_hindexed(nElements, blocks, displ, MpiTypes::type<typename TFwdA::PointType>(), &type);

      // Commit MPI datatype
      MPI_Type_commit(&type);
   }

   template <> template <typename TFwd> void MpiConverterTools<Dimensions::THREED>::buildFwdDatatype(SharedResolution spRes, const Dimensions::Transform::Id fwdDim, TFwd &data, MPI_Datatype &type, const int cpuId)
   {
      // Create map of the local indexes to unique keys
      std::map<Coordinate,Coordinate>  localIdxMap;

      // List of remote keys
      std::set<Coordinate>  remoteKeys;

      // Storage for the coordinate indexes
      Coordinate coord;
      
      // Storage for the simulation wide indexes
      int i_, j_, k_;

      // Storage for the generated key
      Coordinate key;

      //
      // Create the list of local indexes
      //

      // Loop over slow data dimension
      for(int k=0; k < spRes->cpu(FrameworkMacro::id())->dim(fwdDim)->dim<Dimensions::Data::DAT3D>(); ++k)
      {
         // Extract "physical" index of slow data dimension
         k_ = spRes->cpu(FrameworkMacro::id())->dim(fwdDim)->idx<Dimensions::Data::DAT3D>(k);

         // Loop over middle data dimension
         for(int j=0; j < spRes->cpu(FrameworkMacro::id())->dim(fwdDim)->dim<Dimensions::Data::DAT2D>(k); ++j)
         {
            // Extract "physical" index of middle data dimension
            j_ = spRes->cpu(FrameworkMacro::id())->dim(fwdDim)->idx<Dimensions::Data::DAT2D>(j,k);

            // Loop over forward data dimension
            for(int i=0; i < spRes->cpu(FrameworkMacro::id())->dim(fwdDim)->dim<Dimensions::Data::DATF1D>(j,k); ++i)
            {
               // Extract "physical" index of forward data dimension
               i_ = spRes->cpu(FrameworkMacro::id())->dim(fwdDim)->idx<Dimensions::Data::DATF1D>(i,j,k);

               // Combine array indexes into coordinate tuple
               coord = std::tr1::make_tuple(i, j, k);

               // Create key as (1D, 2D, 3D)
               key = std::tr1::make_tuple(i_, j_, k_);

               // add key->coordinate to map
               localIdxMap.insert(std::pair<Coordiante,Coordinate>(key, coord));
            }
         }
      }

      //
      // Create the list of remote indexes in next transform
      //

      // Loop over slow data dimension
      for(int k=0; k < spRes->cpu(cpuId)->dim(Dimensions::Transform::jump<fwdDim,1>::id)->dim<Dimensions::Data::DAT3D>(); ++k)
      {
         // Extract "physical" index of slow data dimension
         k_ = spRes->cpu(cpuId)->dim(Dimensions::Transform::jump<fwdDim,1>::id)->idx<Dimensions::Data::DAT3D>(k);

         // Loop over middle data dimension
         for(int j=0; j < spRes->cpu(cpuId)->dim(Dimensions::Transform::jump<fwdDim,1>::id)->dim<Dimensions::Data::DAT2D>(k); ++j)
         {
            // Extract "physical" index of middle data dimension
            j_ = spRes->cpu(cpuId)->dim(Dimensions::Transform::jump<fwdDim,1>::id)->idx<Dimensions::Data::DAT2D>(j,k);

            // Loop over backward data dimension
            for(int i=0; i < spRes->cpu(cpuId)->dim(Dimensions::Transform::jump<fwdDim,1>::id)->dim<Dimensions::Data::DATB1D>(j,k); ++i)
            {
               // Extract "physical" index of backward data dimension
               i_ = spRes->cpu(cpuId)->dim(Dimensions::Transform::jump<fwdDim,1>::id)->idx<Dimensions::Data::DATB1D>(i,j,k);

               // Create key as (2D, 3D, 1D) indexes (i.e. data gets transposed during communication)
               key = std::tr1::make_tuple(j_, k_, i_);

               // Add key to remote set
               remoteKeys.insert(key);
            }
         }
      }

      // Extract map of shared indexes (stored as keys)
      std::map<Coordinate,Coordinate>  sharedMap;
      MpiConverterTools<Dimensions::THREED>::extractShared(sharedMap, localIdxMap, remoteKeys);

      // Create the datatype
      MpiConverterTools<Dimensions::THREED>::buildType(data, sharedMap, type);
   }

   template <> template <typename TBwd void MpiConverterTools<Dimensions::THREED>::buildBwdDatatype(SharedResolution spRes, const Dimensions::Transform::Id fwdDim, TBwd &data, MPI_Datatype &type, const int cpuId)
   {
      // Create  map of the local indexes to unique keys
      std::map<Coordiante,Coordinate>  localIdxMap;

      // List of remote keys
      std::set<Coordinate>  remoteKeys;

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
      for(int k=0; k < spRes->cpu(FrameworkMacro::id())->dim(Dimensions::Transform::jump<fwdDim,1>::id)->dim<Dimensions::Data:DAT3D>(); ++k)
      {
         // Extract "physical" index of slow data dimension
         k_ = spRes->cpu(FrameworkMacro::id())->dim(Dimensions::Transform::jump<fwdDim,1>::id)->idx<Dimensions::Data:DAT3D>(k);

         // Loop over middle data dimension
         for(int j=0; j < spRes->cpu(FrameworkMacro::id())->dim(Dimensions::Transform::jump<fwdDim,1>::id)->dim<Dimensions::Data:DAT2D>(k); ++j)
         {
            // Extract "physical" index of middle data dimension
            j_ = spRes->cpu(FrameworkMacro::id())->dim(Dimensions::Transform::jump<fwdDim,1>::id)->idx<Dimensions::Data:DAT2D>(j,k);

            // Loop over backward data dimension
            for(int i=0; i < spRes->cpu(FrameworkMacro::id())->dim(Dimensions::Transform::jump<fwdDim,1>::id)->dim<Dimensions::Data:DATB1D>(j,k); ++i)
            {
               // Extract "physical" index of backward data dimension
               i_ = spRes->cpu(FrameworkMacro::id())->dim(Dimensions::Transform::jump<fwdDim,1>::id)->idx<Dimensions::Data:DATB1D>(i,j,k);

               // Combine array indexes into coordinate tuple
               coord = std::tr1::make_tuple(i, j, k);

               // Create key as (2D, 3D, 1D)
               key = std::tr1::make_tuple(j_, k_, i_);

               // add key->coordinate to map
               localIdxMap.insert(std::pair<Coordinate,Coordinate>(key, coord));
            }
         }
      }

      //
      // Create the list of remote indexes
      //

      // Loop over slow data dimension
      for(int k=0; k < spRes->cpu(cpuId)->dim(fwdDim)->dim<Dimensions::Data:DAT3D>(); ++k)
      {
         // Extract "physical" index of slow data dimension
         k_ = spRes->cpu(cpuId)->dim(fwdDim)->idx<Dimensions::Data:DAT3D>(k);

         // Loop over middle data dimension
         for(int j=0; j < spRes->cpu(cpuId)->dim(fwdDim)->dim<Dimensions::Data:DAT2D>(k); ++j)
         {
            // Extract "physical" index of middle data dimension
            j_ = spRes->cpu(cpuId)->dim(fwdDim)->idx<Dimensions::Data:DAT2D>(j,k);

            // Loop over forward data dimension
            for(int i=0; i < spRes->cpu(cpuId)->dim(fwdDim)->dim<Dimensions::Data:DATF1D>(j,k); ++i)
            {
               // Extract "physical" index of forward data dimension
               i_ = spRes->cpu(cpuId)->dim(fwdDim)->idx<Dimensions::Data:DATF1D>(i,j,k);

               // Create key as (1D, 2D, 3D)
               key = std::tr1::make_tuple(i_, j_, k_);

               // Add key to remote set
               remoteKeys.insert(key);
            }
         }
      }

      // Extract map of shared indexes (stored as keys)
      std::map<Coordinate,Coordinate>  sharedMap;
      MpiConverterTools<Dimensions::THREED>::extractShared(sharedMap, localIdxMap, remoteKeys);

      // Create the datatype
      MpiConverterTools<Dimensions::THREED>::buildType(data, sharedMap, type);
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
         template <typename TFwd> void buildFwdDatatype(SharedResolution spRes, const Dimensions::Transform::Id fwdDim, TFwd &data, MPI_Datatype &type, const int cpuId);

         /**
          * @brief Build a backward MPI Datatype
          *
          * @param spRes   Shared Resolution
          * @param fwdDim  Dimension index for forward transform
          * @param data    Input data
          * @param type    Created MPI data type
          * @param cpuId   ID of the CPU
          */
         template <typename TBwd> void buildBwdDatatype(SharedResolution spRes, const Dimensions::Transform::Id fwdDim, TBwd &data, MPI_Datatype &type, const int cpuId);

         /**
          * @brief Extract shared indexes
          *
          * @param sharedMap     Storage for the shared index map
          * @param localIdxMap   Local node key to indexes map 
          * @param remoteKeys    Remote node index keys
          */
         void extractShared(std::map<Coordinate,Coordinate>& sharedMap, const std::map<Coordinate,Coordinate>& localIdxMap, const std::set<Coordinate>& remoteKeys);

         /**
          * @brief Create type
          *
          * @param data       The concerned data
          * @param sharedMap  Shared index map
          * @param type       MPI datatype storage
          */
         template <typename TData> void buildType(TData &data, const std::map<Coordinate,Coordinate>& sharedMap, MPI_Datatype &type);
         
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

   template <> void MpiConverterTools<Dimensions::TWOD>::extractShared(std::mapMpiConverterTools<Dimensions::TWOD>::<Coordinate,MpiConverterTools<Dimensions::TWOD>::Coordinate>& sharedMap, const std::map<MpiConverterTools<Dimensions::TWOD>::Coordinate,MpiConverterTools<Dimensions::TWOD>::Coordinate>& localIdxMap, const std::set<Coordinate>& remoteKeys)
   {
      // List of local index keys 
      std::set<Coordinate>  localKeys;

      // Create map iterator
      std::map<Coordinate,Coordinate>::const_iterator mapIt;
      
      // Extract the set of local keys
      for(mapIt = localIdxMap.begin(); mapIt != localIdxMap.end(); ++mapIt)
      {
         localKeys.insert(mapIt->first);
      }
      
      // Storage for the shared keys
      std::set<Coordinate> sharedKeys;

      // Create the list of common indexes
      std::set_intersection(localKeys.begin(), localKeys.end(), remoteKeys.begin(), remoteKeys.end(), std::inserter(sharedKeys, sharedKeys.begin()));

      // Clear the shared map
      sharedMap.clear();

      // Fill shared map
      std::set<Coordinate>::iterator sit;
      for(sit = sharedKeys.begin(); sit != sharedKeys.end(); sit++)
      {
         sharedMap.insert(std::make_pair(*sit, localIdxMap.find(*sit)->second));
      }
   }

   template <> template <typename TData> void MpiConverterTools<Dimensions::TWOD>::buildType(TData &data, const std::map<MpiConverterTools<Dimensions::TWOD>::Coordinate,MpiConverterTools<Dimensions::TWOD>::Coordinate>& sharedMap, MPI_Datatype &type)
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
      MPI_Type_create_hindexed(nElements, blocks, displ, MpiTypes::type<typename TFwdA::PointType>(), &type);

      // Commit MPI datatype
      MPI_Type_commit(&type);
   }

   template <> template <typename TFwd> void MpiConverterTools<Dimensions::TWOD>::buildFwdDatatype(SharedResolution spRes, const Dimensions::Transform::Id fwdDim, TFwd &data, MPI_Datatype &type, const int cpuId)
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
      for(int j=0; j < spRes->cpu(FrameworkMacro::id())->dim(fwdDim)->dim<Dimensions::Data:DAT2D>(); ++j)
      {
         // Extract "physical" index of slow data dimension
         j_ = spRes->cpu(FrameworkMacro::id())->dim(fwdDim)->idx<Dimensions::Data:DAT2D>(j);

         // Loop over forward data dimension
         for(int i=0; i < spRes->cpu(FrameworkMacro::id())->dim(fwdDim)->dim<Dimensions::Data:DATF1D>(j); ++i)
         {
            // Extract "physical" index of forward data dimension
            i_ = spRes->cpu(FrameworkMacro::id())->dim(fwdDim)->idx<Dimensions::Data:DATF1D>(i,j);

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
      for(int j=0; j < spRes->cpu(cpuId)->dim(fwdDim)->dim<Dimensions::Data:DAT2D>(); ++j)
      {
         // Extract "physical" index of slow data dimension
         j_ = spRes->cpu(cpuId)->dim(fwdDim)->idx<Dimensions::Data:DAT2D>(j);

         // Loop over backward data dimension
         for(int i=0; i < spRes->cpu(cpuId)->dim(Dimensions::Transform::jump<fwdDim,1>::id)->dim<Dimensions::Data:DATB1D>(j); ++i)
         {
            // Extract "physical" index of backward data dimension
            i_ = spRes->cpu(cpuId)->dim(Dimensions::Transform::jump<fwdDim,1>::id)->idx<Dimensions::Data:DATB1D>(i,j);

            // Create key as (2D, 1D)
            key = std::tr1::make_tuple(j_, i_);

            // Add key to remote set
            remoteKeys.insert(key);
         }
      }

      // Extract map of shared indexes (stored as keys)
      std::map<Coordinate,Coordinate>  sharedMap;
      MpiConverterTool<Dimensions::TWOD>::extractShared(sharedMap, localIdxMap, remoteKeys);

      // Create the datatype
      MpiConverterTool<Dimensions::TWOD>::buildType(data, sharedMap, type);
   }

   template <> template <typename TBwd void MpiConverterTools<Dimensions::TWOD>::buildBwdDatatype(SharedResolution spRes, const Dimensions::Transform::Id fwdDim, TBwd &data, MPI_Datatype &type, const int cpuId)
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
      for(int j=0; j < spRes->cpu(FrameworkMacro::id())->dim(Dimensions::Transform::jump<fwdDim,1>::id)->dim<Dimensions::Data:DAT2D>(); ++j)
      {
         // Extract "physical" index of slow data dimension
         j_ = spRes->cpu(FrameworkMacro::id())->dim(Dimensions::Transform::jump<fwdDim,1>::id)->idx<Dimensions::Data:DAT2D>(j);

         // Loop over backward data dimension
         for(int i=0; i < spRes->cpu(FrameworkMacro::id())->dim(Dimensions::Transform::jump<fwdDim,1>::id)->dim<Dimensions::Data:DATB1D>(j); ++i)
         {
            // Extract "physical" index of backward data dimension
            i_ = spRes->cpu(FrameworkMacro::id())->dim(Dimensions::Transform::jump<fwdDim,1>::id)->idx<Dimensions::Data:DATB1D>(i,j);

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
      for(int j=0; j < spRes->cpu(cpuId)->dim(fwdDim)->dim<Dimensions::Data:DAT2D>(); ++j)
      {
         // Extract "physical" index of slow data dimension
         j_ = spRes->cpu(cpuId)->dim(fwdDim)->idx<Dimensions::Data:DAT2D>(j);

         // Loop ver forward data dimension
         for(int i=0; i < spRes->cpu(cpuId)->dim(fwdDim)->dim<Dimensions::Data:DATF1D>(j); ++i)
         {
            // Extract "physical" index of forward data dimension
            i_ = spRes->cpu(cpuId)->dim(fwdDim)->idx<Dimensions::Data:DATF1D>(i,j);

            // Create key as (1D, 2D)
            key = std::tr1::make_pair(i_, j_);

            // Add key to remote set
            remoteKeys.insert(key);
         }
      }

      // Extract map of shared indexes (stored as keys)
      std::map<Coordinate,Coordinate>  sharedMap;
      MpiConverterTools<Dimensions::TWOD>::extractShared(sharedMap, localIdxMap, remoteKeys);

      // Create the datatype
      MpiConverterTools<Dimensions::TWOD>::buildType(data, sharedMap, type);
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
         template <typename TFwd> void buildFwdDatatype(SharedResolution spRes, const Dimensions::Transform::Id fwdDim, TFwd &data, MPI_Datatype &type, const int cpuId);

         /**
          * @brief Build a backward MPI Datatype
          *
          * @param spRes   Shared Resolution
          * @param fwdDim  Dimension index for forward transform
          * @param data    Input data
          * @param type    Created MPI data type
          * @param cpuId   ID of the CPU
          */
         template <typename TBwd> void buildBwdDatatype(SharedResolution spRes, const Dimensions::Transform::Id fwdDim, TBwd &data, MPI_Datatype &type, const int cpuId);

         /**
          * @brief Extract shared indexes
          *
          * @param sharedMap     Storage for the shared index map
          * @param localIdxMap   Local node key to indexes map 
          * @param remoteKeys    Remote node index keys
          */
         void extractShared(std::map<Coordinate,Coordinate>& sharedMap, const std::map<Coordinate,Coordinate>& localIdxMap, const std::set<Coordinate>& remoteKeys);

         /**
          * @brief Create type
          *
          * @param data       The concerned data
          * @param sharedMap  Shared index map
          * @param type       MPI datatype storage
          */
         template <typename TData> void buildType(TData &data, const std::map<Coordinate,Coordinate>& sharedMap, MPI_Datatype &type);
         
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

   template <> void MpiConverterTools<Dimensions::ONED>::extractShared(std::map<MpiConverterTools::<Dimensions::ONED>::Coordinate,MpiConverterTools::<Dimensions::ONED>::Coordinate>& sharedMap, const std::map<MpiConverterTools::<Dimensions::ONED>::Coordinate,MpiConverterTools::<Dimensions::ONED>::Coordinate>& localIdxMap, const std::set<MpiConverterTools::<Dimensions::ONED>::Coordinate>& remoteKeys)
   {
      // List of local index keys 
      std::set<Coordinate>  localKeys;

      // Create map iterator
      std::map<Coordinate,Coordinate>::const_iterator mapIt;
      
      // Extract the set of local keys
      for(mapIt = localIdxMap.begin(); mapIt != localIdxMap.end(); ++mapIt)
      {
         localKeys.insert(mapIt->first);
      }
      
      // Storage for the shared keys
      std::set<Coordinate> sharedKeys;

      // Create the list of common indexes
      std::set_intersection(localKeys.begin(), localKeys.end(), remoteKeys.begin(), remoteKeys.end(), std::inserter(sharedKeys, sharedKeys.begin()));

      // Clear the shared map
      sharedMap.clear();

      // Fill shared map
      std::set<Coordinate>::iterator sit;
      for(sit = sharedKeys.begin(); sit != sharedKeys.end(); sit++)
      {
         sharedMap.insert(std::make_pair(*sit, localIdxMap.find(*sit)->second));
      }
   }

   template <> template <typename TData> void MpiConverterTools<Dimensions::ONED>::buildType(TData &data, const std::map<MpiConverterTools::<Dimensions::ONED>::Coordinate,MpiConverterTools::<Dimensions::ONED>::Coordinate>& sharedMap, MPI_Datatype &type)
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
      MPI_Type_create_hindexed(nElements, blocks, displ, MpiTypes::type<typename TFwdA::PointType>(), &type);

      // Commit MPI datatype
      MPI_Type_commit(&type);
   }

   template <> template <typename TFwd> void MpiConverterTools<Dimensions::ONED>::buildFwdDatatype(SharedResolution spRes, const Dimensions::Transform::Id fwdDim, TFwd &data, MPI_Datatype &type, const int cpuId)
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
      for(int i=0; i < spRes->cpu(FrameworkMacro::id())->dim(fwdDim)->dim<Dimensions::Data:DATF1D>(); ++i)
      {
         // Extract "physical" index of forward data dimension
         i_ = spRes->cpu(FrameworkMacro::id())->dim(fwdDim)->idx<Dimensions::Data:DATF1D>(i);

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
      for(int i=0; i < spRes->cpu(cpuId)->dim(Dimensions::Transform::jump<fwdDim,1>::id)->dim<Dimensions::Data:DATB1D>(); ++i)
      {
         // Extract "physical" index of backward data dimension
         i_ = spRes->cpu(cpuId)->dim(Dimensions::Transform::jump<fwdDim,1>::id)->idx<Dimensions::Data:DATB1D>(i);

         // Set key
         key = i_;

         // Add key to remote set
         remoteKeys.insert(key);
      }

      // Extract map of shared indexes (stored as keys)
      std::map<Coordinate,Coordinate>  sharedMap;
      MpiConverterTools::<Dimensions::ONED>::extractShared(sharedMap, localIdxMap, remoteKeys);

      // Create the datatype
      MpiConverterTools::<Dimensions::ONED>::buildType(data, sharedMap, type);
   }

   template <> template <typename TBwd void MpiConverterTools<Dimensions::ONED>::buildBwdDatatype(SharedResolution spRes, const Dimensions::Transform::Id fwdDim, TBwd &data, MPI_Datatype &type, const int cpuId)
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
      for(int i=0; i < spRes->cpu(FrameworkMacro::id())->dim(Dimensions::Transform::jump<fwdDim,1>::id)->dim<Dimensions::Data:DATB1D>(); ++i)
      {
         // Extract "physical" index of backward data dimension
         i_ = spRes->cpu(FrameworkMacro::id())->dim(Dimensions::Transform::jump<fwdDim,1>::id)->idx<Dimensions::Data:DATB1D>(i);

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
      for(int i=0; i < spRes->cpu(cpuId)->dim(fwdDim)->dim<Dimensions::Data:DATF1D>(); ++i)
      {
         // Extract "physical" index of forward data dimension
         i_ = spRes->cpu(cpuId)->dim(Dimensions::Transform::jump<fwdDim,1>::id)->idx<Dimensions::Data:DATF1D>(i);

         // Set key
         key = i_;

         // Add key to remote set
         remoteKeys.insert(key);
      }

      // Extract map of shared indexes (stored as keys)
      std::map<Coordinate,Coordinate>  sharedMap;
      MpiConverterTools::<Dimensions::ONED>::extractShared(sharedMap, localIdxMap, remoteKeys);

      // Create the datatype
      MpiConverterTools::<Dimensions::ONED>::buildType(data, sharedMap, type);
   }

}
}

#endif // MPICONVERTERTOOLS_HPP
