/** 
 * @file MpiConverterTools.cpp
 * @brief Source of the tools for the MPI converter
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Communicators/Converters/MpiConverterTools.hpp"

// Project includes
//

namespace QuICC {

namespace Parallel {

   // 
   // Three dimensional
   //

   void MpiConverterTools<Dimensions::THREED>::buildFwdCpuMap(std::map<Coordinate,Coordinate>& sharedMap, const std::map<Coordinate,Coordinate>& localIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, const int cpuId)
   {
      // List of remote keys
      std::set<Coordinate>  remoteKeys;

      // Storage for the simulation wide indexes
      int i_, j_, k_;

      // Storage for the generated key
      Coordinate key;

      // MPI error code
      int ierr;

      // Number of coordinates
      int nCoords = -1;

      //
      // Create the list of remote indexes in next transform
      //

      // Remote is also local CPU
      MatrixI  matRemote;
      if(cpuId == FrameworkMacro::transformId(fwdDim))
      {
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

                  // Create key as (2D, 3D, 1D) indexes (i.e. data gets transposed during communication)
                  key = spRes->counter()->makeKey(Dimensions::jump(fwdDim,1), i_, j_, k_);

                  // Add key to remote set
                  remoteKeys.insert(key);
               }
            }
         }

         // Convert remote keys to matrix to send througg MPI
         matRemote.resize(3, remoteKeys.size());
         int col = 0; 
         for(std::set<Coordinate>::iterator it = remoteKeys.begin(); it != remoteKeys.end(); ++it)
         {
            matRemote(0,col) = std::tr1::get<0>(*it);
            matRemote(1,col) = std::tr1::get<1>(*it);
            matRemote(2,col) = std::tr1::get<2>(*it);
            col++;
         }
   
         // Broadcast size
         nCoords = remoteKeys.size();
         FrameworkMacro::syncTransform(fwdDim);
         ierr = MPI_Bcast(&nCoords, 1, MPI_INT, cpuId, FrameworkMacro::transformComm(fwdDim));
         FrameworkMacro::check(ierr, 951);

         // Broadcast data
         FrameworkMacro::syncTransform(fwdDim);
         ierr = MPI_Bcast(matRemote.data(), matRemote.size(), MPI_INT, cpuId, FrameworkMacro::transformComm(fwdDim)); 
         FrameworkMacro::check(ierr, 952);

      // Remote CPU needs to generate list 
      } else
      {
         // Get size
         FrameworkMacro::syncTransform(fwdDim);
         ierr = MPI_Bcast(&nCoords, 1, MPI_INT, cpuId, FrameworkMacro::transformComm(fwdDim));
         FrameworkMacro::check(ierr, 953);

         // Resize storage
         matRemote.resize(3, nCoords);

         // Get remote keys as matrix
         FrameworkMacro::syncTransform(fwdDim);
         ierr = MPI_Bcast(matRemote.data(), matRemote.size(), MPI_INT, cpuId, FrameworkMacro::transformComm(fwdDim)); 
         FrameworkMacro::check(ierr, 954);

         // Convert matrix to remoteKeys set
         for(int i = 0; i < matRemote.cols(); i++)
         {
            key = std::tr1::make_tuple(matRemote(0,i), matRemote(1,i), matRemote(2,i));
            remoteKeys.insert(key);
         }
      }

      // Extract map of shared indexes (stored as keys)
      sharedMap.clear();
      MpiConverterTools<Dimensions::THREED>::extractShared(sharedMap, localIdxMap, remoteKeys);
   }

   void MpiConverterTools<Dimensions::THREED>::buildBwdCpuMap(std::map<Coordinate,Coordinate>& sharedMap, const std::map<Coordinate,Coordinate>&  localIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, const int cpuId)
   {
      // List of remote keys
      std::set<Coordinate>  remoteKeys;
      
      // Storage for the simulation wide indexes
      int i_, j_, k_;

      // Storage for the generated key
      Coordinate key;

      // MPI error code
      int ierr;

      // Number of coordinates
      int nCoords = -1;

      //
      // Create the list of remote indexes
      //

      // Remote is also local CPU
      MatrixI  matRemote;
      if(cpuId == FrameworkMacro::transformId(fwdDim))
      {
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

                  // Create key as (1D, 2D, 3D)
                  key = spRes->counter()->makeKey(fwdDim, i_, j_, k_);

                  // Add key to remote set
                  remoteKeys.insert(key);
               }
            }
         }

         // Convert remote keys to matrix to send through MPI
         matRemote.resize(3, remoteKeys.size());
         int col = 0; 
         for(std::set<Coordinate>::iterator it = remoteKeys.begin(); it != remoteKeys.end(); ++it)
         {
            matRemote(0,col) = std::tr1::get<0>(*it);
            matRemote(1,col) = std::tr1::get<1>(*it);
            matRemote(2,col) = std::tr1::get<2>(*it);
            col++;
         }

         // Broadcast size
         nCoords = remoteKeys.size();
         FrameworkMacro::syncTransform(fwdDim);
         ierr = MPI_Bcast(&nCoords, 1, MPI_INT, cpuId, FrameworkMacro::transformComm(fwdDim));
         FrameworkMacro::check(ierr, 961);

         // Broadcast data
         FrameworkMacro::syncTransform(fwdDim);
         ierr = MPI_Bcast(matRemote.data(), matRemote.size(), MPI_INT, cpuId, FrameworkMacro::transformComm(fwdDim)); 
         FrameworkMacro::check(ierr, 962);

      // Remote CPU needs to generate list 
      } else
      {
         // Get size
         FrameworkMacro::syncTransform(fwdDim);
         ierr = MPI_Bcast(&nCoords, 1, MPI_INT, cpuId, FrameworkMacro::transformComm(fwdDim));
         FrameworkMacro::check(ierr, 963);

         matRemote.resize(3, nCoords);

         // Get remot ekeys as matrix
         FrameworkMacro::syncTransform(fwdDim);
         ierr = MPI_Bcast(matRemote.data(), matRemote.size(), MPI_INT, cpuId, FrameworkMacro::transformComm(fwdDim)); 
         FrameworkMacro::check(ierr, 964);

         // Convert matrix to remoteKeys set
         for(int i = 0; i < matRemote.cols(); i++)
         {
            key = std::tr1::make_tuple(matRemote(0,i), matRemote(1,i), matRemote(2,i));
            remoteKeys.insert(key);
         }
      }

      // Extract map of shared indexes (stored as keys)
      sharedMap.clear();
      MpiConverterTools<Dimensions::THREED>::extractShared(sharedMap, localIdxMap, remoteKeys);
   }

   void MpiConverterTools<Dimensions::THREED>::extractShared(std::map<MpiConverterTools<Dimensions::THREED>::Coordinate,MpiConverterTools<Dimensions::THREED>::Coordinate>& sharedMap, const std::map<MpiConverterTools<Dimensions::THREED>::Coordinate,MpiConverterTools<Dimensions::THREED>::Coordinate>& localIdxMap, const std::set<MpiConverterTools<Dimensions::THREED>::Coordinate>& remoteKeys)
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

   // 
   // Two dimensional
   //

   void MpiConverterTools<Dimensions::TWOD>::buildFwdCpuMap(std::map<Coordinate,Coordinate>& sharedMap, const std::map<Coordinate,Coordinate>& localIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, const int cpuId)
   {
      // List of remote keys
      std::set<Coordinate>  remoteKeys;
      
      // Storage for the simulation wide indexes
      int i_, j_;

      // Storage for the generated key
      Coordinate key;

      // MPI error code
      int ierr;

      // Number of coordinates
      int nCoords = -1;

      //
      // Create the list of remote indexes in next transform
      //

      // Remote is also local CPU
      MatrixI  matRemote;
      if(cpuId == FrameworkMacro::transformId(fwdDim))
      {
         // Loop over middle data dimension
         for(int j=0; j < spRes->cpu()->dim(Dimensions::jump(fwdDim,1))->dim<Dimensions::Data::DAT2D>(); ++j)
         {
            // Extract "physical" index of middle data dimension
            j_ = spRes->cpu()->dim(Dimensions::jump(fwdDim,1))->idx<Dimensions::Data::DAT2D>(j);

            // Loop over backward data dimension
            for(int i=0; i < spRes->cpu()->dim(Dimensions::jump(fwdDim,1))->dim<Dimensions::Data::DATB1D>(); ++i)
            {
               // Extract "physical" index of backward data dimension
               i_ = spRes->cpu()->dim(Dimensions::jump(fwdDim,1))->idx<Dimensions::Data::DATB1D>(i);

               // Create key as (2D, 3D, 1D) indexes (i.e. data gets transposed during communication)
               key = spRes->counter()->makeKey(Dimensions::jump(fwdDim,1), i_, j_);

               // Add key to remote set
               remoteKeys.insert(key);
            }
         }

         // Convert remote keys to matrix to send througg MPI
         matRemote.resize(2, remoteKeys.size());
         int col = 0; 
         for(std::set<Coordinate>::iterator it = remoteKeys.begin(); it != remoteKeys.end(); ++it)
         {
            matRemote(0,col) = it->first;
            matRemote(1,col) = it->second;
            col++;
         }
   
         // Broadcast size
         nCoords = remoteKeys.size();
         FrameworkMacro::syncTransform(fwdDim);
         ierr = MPI_Bcast(&nCoords, 1, MPI_INT, cpuId, FrameworkMacro::transformComm(fwdDim));
         FrameworkMacro::check(ierr, 951);

         // Broadcast data
         FrameworkMacro::syncTransform(fwdDim);
         ierr = MPI_Bcast(matRemote.data(), matRemote.size(), MPI_INT, cpuId, FrameworkMacro::transformComm(fwdDim)); 
         FrameworkMacro::check(ierr, 952);

      // Remote CPU needs to generate list 
      } else
      {
         // Get size
         FrameworkMacro::syncTransform(fwdDim);
         ierr = MPI_Bcast(&nCoords, 1, MPI_INT, cpuId, FrameworkMacro::transformComm(fwdDim));
         FrameworkMacro::check(ierr, 953);

         matRemote.resize(2, nCoords);

         // Get remote keys as matrix
         FrameworkMacro::syncTransform(fwdDim);
         ierr = MPI_Bcast(matRemote.data(), matRemote.size(), MPI_INT, cpuId, FrameworkMacro::transformComm(fwdDim)); 
         FrameworkMacro::check(ierr, 954);

         // Convert matrix to remoteKeys set
         for(int i = 0; i < matRemote.cols(); i++)
         {
            key = std::make_pair(matRemote(0,i), matRemote(1,i));
            remoteKeys.insert(key);
         }
      }

      // Extract map of shared indexes (stored as keys)
      sharedMap.clear();
      MpiConverterTools<Dimensions::TWOD>::extractShared(sharedMap, localIdxMap, remoteKeys);
   }

   void MpiConverterTools<Dimensions::TWOD>::buildBwdCpuMap(std::map<Coordinate,Coordinate>& sharedMap, const std::map<Coordinate,Coordinate>& localIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, const int cpuId)
   {
      // List of remote keys
      std::set<Coordinate>  remoteKeys;
      
      // Storage for the simulation wide indexes
      int i_, j_;

      // Storage for the generated key
      Coordinate key;

      // MPI error code
      int ierr;

      // Number of coordinates
      int nCoords = -1;

      //
      // Create the list of remote indexes
      //

      // Remote is also local CPU
      MatrixI  matRemote;
      if(cpuId == FrameworkMacro::transformId(fwdDim))
      {
         // Loop over middle data dimension
         for(int j = 0; j < spRes->cpu()->dim(fwdDim)->dim<Dimensions::Data::DAT2D>(); ++j)
         {
            // Extract "physical" index of middle data dimension
            j_ = spRes->cpu()->dim(fwdDim)->idx<Dimensions::Data::DAT2D>(j);

            // Loop over forward data dimension
            for(int i = 0; i < spRes->cpu()->dim(fwdDim)->dim<Dimensions::Data::DATF1D>(); ++i)
            {
               // Extract "physical" index of forward data dimension
               i_ = spRes->cpu()->dim(fwdDim)->idx<Dimensions::Data::DATF1D>(i);

               // Create key as (1D, 2D, 3D)
               key = spRes->counter()->makeKey(fwdDim, i_, j_);

               // Add key to remote set
               remoteKeys.insert(key);
            }
         }

         // Convert remote keys to matrix to send through MPI
         matRemote.resize(2, remoteKeys.size());
         int col = 0; 
         for(std::set<Coordinate>::iterator it = remoteKeys.begin(); it != remoteKeys.end(); ++it)
         {
            matRemote(0,col) = it->first;
            matRemote(1,col) = it->second;
            col++;
         }

         // Broadcast size
         nCoords = remoteKeys.size();
         FrameworkMacro::syncTransform(fwdDim);
         ierr = MPI_Bcast(&nCoords, 1, MPI_INT, cpuId, FrameworkMacro::transformComm(fwdDim));
         FrameworkMacro::check(ierr, 961);

         // Broadcast data
         FrameworkMacro::syncTransform(fwdDim);
         ierr = MPI_Bcast(matRemote.data(), matRemote.size(), MPI_INT, cpuId, FrameworkMacro::transformComm(fwdDim)); 
         FrameworkMacro::check(ierr, 962);

      // Remote CPU needs to generate list 
      } else
      {
         // Get size
         FrameworkMacro::syncTransform(fwdDim);
         ierr = MPI_Bcast(&nCoords, 1, MPI_INT, cpuId, FrameworkMacro::transformComm(fwdDim));
         FrameworkMacro::check(ierr, 963);

         matRemote.resize(2, nCoords);

         // Get remot ekeys as matrix
         FrameworkMacro::syncTransform(fwdDim);
         ierr = MPI_Bcast(matRemote.data(), matRemote.size(), MPI_INT, cpuId, FrameworkMacro::transformComm(fwdDim)); 
         FrameworkMacro::check(ierr, 964);

         // Convert matrix to remoteKeys set
         for(int i = 0; i < matRemote.cols(); i++)
         {
            key = std::make_pair(matRemote(0,i), matRemote(1,i));
            remoteKeys.insert(key);
         }
      }

      // Extract map of shared indexes (stored as keys)
      sharedMap.clear();
      MpiConverterTools<Dimensions::TWOD>::extractShared(sharedMap, localIdxMap, remoteKeys);
   }

   void MpiConverterTools<Dimensions::TWOD>::extractShared(std::map<MpiConverterTools<Dimensions::TWOD>::Coordinate,MpiConverterTools<Dimensions::TWOD>::Coordinate>& sharedMap, const std::map<MpiConverterTools<Dimensions::TWOD>::Coordinate,MpiConverterTools<Dimensions::TWOD>::Coordinate>& localIdxMap, const std::set<Coordinate>& remoteKeys)
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

   // 
   // One dimensional
   //

   void MpiConverterTools<Dimensions::ONED>::extractShared(std::map<MpiConverterTools<Dimensions::ONED>::Coordinate,MpiConverterTools<Dimensions::ONED>::Coordinate>& sharedMap, const std::map<MpiConverterTools<Dimensions::ONED>::Coordinate,MpiConverterTools<Dimensions::ONED>::Coordinate>& localIdxMap, const std::set<MpiConverterTools<Dimensions::ONED>::Coordinate>& remoteKeys)
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

}
}
