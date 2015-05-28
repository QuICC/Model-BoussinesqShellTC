/** 
 * @file IForwardGrouper.cpp
 * @brief Source of the implementation of the equation wise forward transform grouper
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes
//
#include <set>

// Class include
//
#include "TransformGroupers/IForwardGrouper.hpp"

// Project includes
//
#include "TransformConfigurators/TransformStepsMacro.h"

namespace GeoMHDiSCC {

namespace Transform {

   IForwardGrouper::IForwardGrouper()
      : split(Splitting::Locations::NONE)
   {
   }

   IForwardGrouper::~IForwardGrouper()
   {
   }

   ArrayI IForwardGrouper::namePacks1D(const std::vector<IntegratorTree>& integratorTree)
   {
      // Create list of packet sizes
      std::set<int>  list;

      std::vector<IntegratorTree>::const_iterator treeIt;
      for(treeIt = integratorTree.begin(); treeIt != integratorTree.end(); ++treeIt)
      {
         int counter = treeIt->nEdges2D();
         list.insert(counter);

         this->mNamedPacks1D.insert(std::make_pair(std::make_pair(treeIt->name(),treeIt->comp()), counter));
      }

      // Initialise the number of packs
      ArrayI packs(list.size());

      // Set packet sizes
      std::set<int>::iterator it;
      int i = 0;
      for(it = list.begin(); it != list.end(); ++it, ++i)
      {
         packs(i) = *it;
      }

      return packs;
   }

   ArrayI IForwardGrouper::namePacks2D(const  std::vector<IntegratorTree>& integratorTree)
   {  
      // Create list of packet sizes
      std::set<int>  list;

      // Loop over all edges
      std::vector<IntegratorTree>::const_iterator treeIt;
      for(treeIt = integratorTree.begin(); treeIt != integratorTree.end(); ++treeIt)
      {
         int counter = treeIt->nEdges3D();
         list.insert(counter);

         this->mNamedPacks2D.insert(std::make_pair(std::make_pair(treeIt->name(),treeIt->comp()), counter));
      }

      // Initialise the number of packs
      ArrayI packs(list.size());

      // Set packet sizes
      std::set<int>::iterator it;
      int i = 0;
      for(it = list.begin(); it != list.end(); ++it, ++i)
      {
         packs(i) = *it;
      }

      return packs;
   }

   ArrayI IForwardGrouper::groupPacks1D(const std::vector<IntegratorTree>& integratorTree)
   {
      // Initialise the number of packs
      ArrayI packs = this->namePacks1D(integratorTree);

      // Resize packs to single value
      packs.resize(1);
      packs.setConstant(0);

      for(std::map<FieldIdType, int>::const_iterator it = this->mNamedPacks1D.begin(); it != this->mNamedPacks1D.end(); ++it)
      {
         packs(0) = packs(0) + it->second;
      }

      return packs;
   }

   ArrayI IForwardGrouper::groupPacks2D(const std::vector<IntegratorTree>& integratorTree)
   {  
      // Initialise the number of packs
      ArrayI packs = this->namePacks2D(integratorTree);

      // Resize packs to single value
      packs.resize(1);
      packs.setConstant(0);

      for(std::map<FieldIdType, int>::const_iterator it = this->mNamedPacks2D.begin(); it != this->mNamedPacks2D.end(); ++it)
      {
         packs(0) = packs(0) + it->second;
      }

      return packs;
   }

}
}
