/** 
 * @file IForwardGrouper3D.cpp
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
#include "TransformGroupers/IForwardGrouper3D.hpp"

// Project includes
//
#include "TransformConfigurators/TransformStepsMacro.h"

namespace GeoMHDiSCC {

namespace Transform {

   IForwardGrouper3D::IForwardGrouper3D()
      : IForwardGrouper2D()
   {
   }

   IForwardGrouper3D::~IForwardGrouper3D()
   {
   }

   ArrayI IForwardGrouper3D::namePacks2D(const  std::vector<TransformTree>& integratorTree)
   {  
      // Create list of packet sizes
      std::set<int>  list;

      // Loop over all edges
      std::vector<TransformTree>::const_iterator treeIt;
      for(treeIt = integratorTree.begin(); treeIt != integratorTree.end(); ++treeIt)
      {
         int counter = treeIt->nEdges(0);
         list.insert(counter);

         this->mNamedPacks2D.insert(std::make_pair(std::make_pair(treeIt->name(),treeIt->comp<FieldComponents::Physical::Id>()), counter));
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

   ArrayI IForwardGrouper3D::groupPacks2D(const std::vector<TransformTree>& integratorTree)
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
