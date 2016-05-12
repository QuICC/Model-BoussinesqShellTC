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
#include "Exceptions/Exception.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   void TransformTreeTools::generateTrees(std::vector<TransformTree>& rTrees, const std::map<PhysicalNames::Id, std::vector<TransformPath> >& branches)
   {
      //
      // First stage: Collapse calculations with same input into same tree branch
      //
      buildTrees(rTrees, branches);

      //
      // Second stage: Set proper recovery and hold flags on the tree nodes
      //
      finalizeTrees(rTrees);

      // Third stage: Combine arithmetic operations into single tree branch
      optimizeTrees(rTrees);
   }

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

   void TransformTreeTools::buildTrees(std::vector<TransformTree>& rTrees, const std::map<PhysicalNames::Id, std::vector<TransformPath> >& branches)
   {
      //
      // Construct the trees from the transfrom paths
      //
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

   void TransformTreeTools::setInputInfoEdge(TransformTreeEdge& edge)
   {
      // Get range of edges
      TransformTreeEdge::EdgeType_range rangeIt = edge.rEdgeRange();

      // Initialize recover and hold
      int recover = 0;
      int hold = std::distance(rangeIt.first, rangeIt.second) - 1;

      // Loop over edges
      for(TransformTreeEdge::EdgeType_iterator edgeIt = rangeIt.first; edgeIt != rangeIt.second; ++edgeIt, ++recover, --hold)
      {
         // Set recover and hold info
         edgeIt->setInputInfo(recover, hold);

         // Recursively go to next level
         setInputInfoEdge(*edgeIt);
      }
   }

   void TransformTreeTools::setOutputOrder(TransformTreeEdge& edge, std::set<std::pair<FieldType::Id,std::vector<int> > >& outIds)
   {
      // Get range of edges
      TransformTreeEdge::EdgeType_range rangeIt = edge.rEdgeRange();

      // Check for last recurrence
      if(std::distance(rangeIt.first->edgeRange().first, rangeIt.first->edgeRange().second) > 0)
      {
         // Loop over edges
         for(TransformTreeEdge::EdgeType_iterator edgeIt = rangeIt.first; edgeIt != rangeIt.second; ++edgeIt)
         {
            // Recursively go to next level
            setOutputOrder(*edgeIt, outIds);
         }
      
      // Reached last level
      } else
      {
         // Loop over edges
         for(TransformTreeEdge::EdgeType_iterator edgeIt = rangeIt.first; edgeIt != rangeIt.second; ++edgeIt)
         {
            // Try to add to unique output fields
            std::pair<std::set<std::pair<FieldType::Id,std::vector<int> > >::iterator,bool> outIt = outIds.insert(std::make_pair(edgeIt->fieldId(), edgeIt->outIds()));

            // If first entry, replace arithmetics with SET or SETNEG
            if(outIt.second)
            {
               if(edgeIt->arithId() == Arithmetics::ADD)
               {
                  edgeIt->setArithId(Arithmetics::SET);

               } else if(edgeIt->arithId() == Arithmetics::SUB)
               {
                  edgeIt->setArithId(Arithmetics::SETNEG);
               } else
               {
                  throw Exception("Unknown arithmetic operation during transform tree construction!");
               }
            } else
            {
                  throw Exception("tHIS HAPPENS!");
            }
         }
      }
   }

   void TransformTreeTools::finalizeTrees(std::vector<TransformTree>& rTrees)
   {
      // Loop over all trees
      std::set<std::pair<FieldType::Id,std::vector<int> > > outIds;
      PhysicalNames::Id curName = static_cast<PhysicalNames::Id>(0);
      for(std::vector<TransformTree>::iterator treeIt = rTrees.begin(); treeIt != rTrees.end(); ++treeIt)
      {
         // Recursively set hold and recovery flags
         setInputInfoEdge(treeIt->rRoot());

         // Recursively order operations by replacing first by SET or SETNEG
         if(curName != treeIt->name())
         {
            outIds.clear();
         }
         setOutputOrder(treeIt->rRoot(), outIds);
      }
   }

   void TransformTreeTools::optimizeTrees(std::vector<TransformTree>& rTrees)
   {
   }
}
}
