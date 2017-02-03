/** 
 * @file IBackwardGrouper2D.cpp
 * @brief Source of the implementation of the equation wise forward transform grouper in 2D space
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
#include "TransformGroupers/IBackwardGrouper2D.hpp"

// Project includes
//
#include "TransformConfigurators/TransformStepsMacro.h"

namespace QuICC {

namespace Transform {

   IBackwardGrouper2D::IBackwardGrouper2D()
      : split(Splitting::Locations::NONE)
   {
   }

   IBackwardGrouper2D::~IBackwardGrouper2D()
   {
   }

   ArrayI IBackwardGrouper2D::namePacks1D(const std::vector<TransformTree>& projectorTree)
   {
      // Create list of packet sizes
      std::set<int>  list;

      // Loop over all edges
      std::vector<TransformTree>::const_iterator treeIt;
      for(treeIt = projectorTree.begin(); treeIt != projectorTree.end(); ++treeIt)
      {
         int counter = treeIt->nEdges(0);
         list.insert(counter);

         this->mNamedPacks1D.insert(std::make_pair(std::make_pair(treeIt->name(),treeIt->comp<FieldComponents::Spectral::Id>()), counter));
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

   ArrayI IBackwardGrouper2D::groupPacks1D(const std::vector<TransformTree>& projectorTree)
   {
      // Initialise the number of packs
      ArrayI packs = this->namePacks1D(projectorTree);

      // Resize packs to single value
      packs.resize(1);
      packs.setConstant(0);

      for(std::map<FieldIdType, int>::const_iterator it = this->mNamedPacks1D.begin(); it != this->mNamedPacks1D.end(); ++it)
      {
         packs(0) = packs(0) + it->second;
      }

      return packs;
   }

}
}
