/** 
 * @file TransformTreeTools.cpp
 * @brief Source of the implementation of tools to work with projector transform trees
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Debug includes
//
#include "Debug/DebuggerMacro.h"

// Configuration includes
//

// System includes
//
#include <set>
#include <vector>
#include <map>
#include <tr1/tuple>

// External includes
//

// Class include
//
#include "TransformConfigurators/TransformTreeTools.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

   void TransformTreeTools::growTree(std::map<PhysicalNames::Id, std::vector<TransformPath> >::const_iterator nameIt, std::set<int>::const_iterator compIt, const int dim, std::vector<int>& path, TransformTreeEdge &rPrev)
   {
      // Typedefs to simplify notation
      typedef std::map<int, int> OpMap;
      typedef std::multimap<int, std::tr1::tuple<std::vector<int>, FieldType::Id, Arithmetics::Id> > OpMulti;

      bool endRecursion = false;

      // Extract unique operators
      OpMap op;
      std::pair<OpMap::iterator,bool> opPairIt;
      OpMulti multi;
      for(std::vector<TransformPath>::const_iterator branchIt = nameIt->second.begin(); branchIt != nameIt->second.end(); ++branchIt)
      {
         // Check path
         bool isSame = true;
         for(int i = 0; i < dim; i++)
         {
            isSame = isSame && (branchIt->edge(i).opId() == path.at(i));
         }

         // Count path multiplicity
         if(branchIt->startId() == *compIt && isSame)
         {
            opPairIt = op.insert(std::make_pair(branchIt->edge(dim).opId(), 0));
            opPairIt.first->second += 1;

            // Check for end of recursion
            if(dim == branchIt->nEdges()-1)
            {
               multi.insert(std::make_pair(branchIt->edge(dim).opId(), std::tr1::make_tuple(branchIt->edge(dim).outId(), branchIt->fieldId(), branchIt->edge(dim).arithId())));
               endRecursion = true;
            }
         }
      }

      // Check for end of recursion
      if(endRecursion)
      {
         // Create edges
         for(OpMap::const_iterator opIt = op.begin(); opIt != op.end(); ++opIt)
         {
            std::pair<OpMulti::const_iterator, OpMulti::const_iterator> multiRange = multi.equal_range(opIt->first);
            for(OpMulti::const_iterator multiIt = multiRange.first; multiIt != multiRange.second; ++multiIt)
            {
               TransformTreeEdge &rNext = rPrev.addEdge(opIt->first, opIt->second);

               rNext.setEnd(std::tr1::get<0>(multiIt->second), std::tr1::get<1>(multiIt->second), std::tr1::get<2>(multiIt->second));
            }
         }

      } else
      {
         // Create edges
         for(OpMap::const_iterator opIt = op.begin(); opIt != op.end(); ++opIt)
         {
            TransformTreeEdge &rNext = rPrev.addEdge(opIt->first, opIt->second);

            // Add last element to path
            path.push_back(opIt->first);

            // Grow the tree
            growTree(nameIt, compIt, dim+1, path, rNext);

            // Remove last element from path
            path.pop_back();
         }
      }
   }

   void TransformTreeTools::generateTrees(std::vector<TransformTree>& rTrees, const std::map<PhysicalNames::Id, std::vector<TransformPath> >& branches)
   {
      // Loop over all physical fields
      std::map<PhysicalNames::Id, std::vector<TransformPath> >::const_iterator nameIt;
      for(nameIt = branches.begin(); nameIt != branches.end(); ++nameIt)
      {
         // Generate list of unique components
         std::vector<TransformPath>::const_iterator branchIt;
         std::set<int> components;
         for(branchIt = nameIt->second.begin(); branchIt != nameIt->second.end(); ++branchIt)
         {
            components.insert(branchIt->startId());
         }

         // Initialise component trees
         std::vector<int> path;
         std::set<int>::const_iterator compIt;
         for(compIt = components.begin(); compIt != components.end(); ++compIt)
         {
            // Initialize tree
            rTrees.push_back(TransformTree(nameIt->first, *compIt));
            
            // Grow the tree recursively
            path.clear();
            growTree(nameIt, compIt, 0, path, rTrees.back().rRoot()); 
         }
      }
   }
}
}
