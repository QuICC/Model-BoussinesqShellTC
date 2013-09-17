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

namespace GeoMHDiSCC {

namespace Parallel {

   // 
   // Three dimensional
   //

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
