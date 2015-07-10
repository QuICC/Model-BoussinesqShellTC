/** 
 * @file IBackwardGrouper3D.cpp
 * @brief Source of the implementation of the equation wise forward transform grouper in 3D space
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
#include "TransformGroupers/IBackwardGrouper3D.hpp"

// Project includes
//
#include "TransformConfigurators/TransformStepsMacro.h"

namespace GeoMHDiSCC {

namespace Transform {

   IBackwardGrouper3D::IBackwardGrouper3D()
      : IBackwardGrouper2D()
   {
   }

   IBackwardGrouper3D::~IBackwardGrouper3D()
   {
   }

   ArrayI IBackwardGrouper3D::namePacks2D(const std::vector<ProjectorTree>& projectorTree)
   {  
      // Create list of packet sizes
      std::set<int>  list;

      // Loop over all edges
      std::vector<ProjectorTree>::const_iterator treeIt;
      for(treeIt = projectorTree.begin(); treeIt != projectorTree.end(); ++treeIt)
      {
         int counter = treeIt->nPartEdges();
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

   ArrayI IBackwardGrouper3D::groupPacks2D(const std::vector<ProjectorTree>& projectorTree)
   {  
      // Initialise the number of packs
      ArrayI packs = this->namePacks2D(projectorTree);

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
