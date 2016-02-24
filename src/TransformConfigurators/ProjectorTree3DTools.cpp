/** 
 * @file ProjectorTree3DTools.cpp
 * @brief Source of the implementation of tools to work with projector transform trees for 3D space
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Debug includes
//
#include "Debug/DebuggerMacro.h"

// Configuration includes
//

// System includes
//
#include <tr1/tuple>

// External includes
//

// Class include
//
#include "TransformConfigurators/ProjectorTree3DTools.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

   void ProjectorTree3DTools::generateTrees(std::vector<ProjectorTree3D>& rTrees, const std::map<PhysicalNames::Id, std::vector<ProjectorBranch3D> >& branches)
   {
      // Debugging info 
      DebuggerMacro_enter("Generating projector trees", 1);

      // Loop over all physical fields
      std::map<PhysicalNames::Id, std::vector<ProjectorBranch3D> >::const_iterator nameIt;
      for(nameIt = branches.begin(); nameIt != branches.end(); ++nameIt)
      {
         // Loop over branches
         std::vector<ProjectorBranch3D>::const_iterator branchIt;
         std::set<FieldComponents::Spectral::Id> components;
         for(branchIt = nameIt->second.begin(); branchIt != nameIt->second.end(); ++branchIt)
         {
            // Generate list of unique components
            components.insert(branchIt->specId());
         }

         // Initialise component trees
         std::set<FieldComponents::Spectral::Id>::const_iterator compIt;
         for(compIt = components.begin(); compIt != components.end(); ++compIt)
         {
            // Initilize tree
            rTrees.push_back(ProjectorTree3D(nameIt->first, *compIt));

            DebuggerMacro_showValue("Field: ", 2, nameIt->first);
            DebuggerMacro_showValue("Component: ", 2, *compIt);

            // Extract unique 1D operators
            std::map<ProjSpecId, int> op1D;
            std::pair<std::map<ProjSpecId,int>::iterator,bool> op1DPairIt;
            for(branchIt = nameIt->second.begin(); branchIt != nameIt->second.end(); ++branchIt)
            {
               if(branchIt->specId() == *compIt)
               {
                  op1DPairIt = op1D.insert(std::make_pair(branchIt->projSpecId(), 0));
                  op1DPairIt.first->second += 1;
               }
            }

            // Create 1D edges
            std::map<ProjSpecId,int>::const_iterator op1DIt;
            for(op1DIt = op1D.begin(); op1DIt != op1D.end(); ++op1DIt)
            {
               ProjectorSpecEdge &rEdge1D = rTrees.back().addEdge(op1DIt->first, op1DIt->second);
               DebuggerMacro_showValue("Edge operator: ", 3, op1DIt->first);
               DebuggerMacro_showValue("Edge weight: ", 3, op1DIt->second);

               // Extract unique 2D Operators
               std::map<ProjPartId, int> op2D;
               std::pair<std::map<ProjPartId,int>::iterator,bool> op2DPairIt;
               for(branchIt = nameIt->second.begin(); branchIt != nameIt->second.end(); ++branchIt)
               {
                  if(branchIt->specId() == *compIt && branchIt->projSpecId() == op1DIt->first)
                  {
                     op2DPairIt = op2D.insert(std::make_pair(branchIt->projPartId(), 0));
                     op2DPairIt.first->second += 1;
                  }
               }

               // Create 2D edges
               std::map<ProjPartId,int>::const_iterator op2DIt;
               for(op2DIt = op2D.begin(); op2DIt != op2D.end(); ++op2DIt)
               {
                  ProjectorPartEdge &rEdge2D = rEdge1D.addEdge(op2DIt->first, op2DIt->second);
                  DebuggerMacro_showValue("Edge operator: ", 4, op2DIt->first);
                  DebuggerMacro_showValue("Edge weight: ", 4, op2DIt->second);

                  // Extract unique 3D operators
                  std::map<ProjPhysId, int> op3D;
                  std::multimap<ProjPhysId, std::tr1::tuple<std::vector<FieldComponents::Physical::Id>,FieldType::Id,Arithmetics::Id> > op3DPhys;
                  std::pair<std::map<ProjPhysId,int>::iterator,bool> op3DPairIt;
                  for(branchIt = nameIt->second.begin(); branchIt != nameIt->second.end(); ++branchIt)
                  {
                     if(branchIt->specId() == *compIt && branchIt->projSpecId() == op1DIt->first && branchIt->projPartId() == op2DIt->first)
                     {
                        op3DPairIt = op3D.insert(std::make_pair(branchIt->projPhysId(), 0));
                        op3DPairIt.first->second += 1;
                        op3DPhys.insert(std::make_pair(branchIt->projPhysId(), std::tr1::make_tuple(branchIt->physId(), branchIt->fieldId(), branchIt->arithId())));
                     }
                  }

                  // Create 3D edges
                  std::map<ProjPhysId,int>::const_iterator op3DIt;
                  std::multimap<ProjPhysId, std::tr1::tuple<std::vector<FieldComponents::Physical::Id>,FieldType::Id,Arithmetics::Id> >::iterator op3DPhysIt;
                  std::pair<std::multimap<ProjPhysId, std::tr1::tuple<std::vector<FieldComponents::Physical::Id>,FieldType::Id,Arithmetics::Id> >::iterator,std::multimap<ProjPhysId, std::tr1::tuple<std::vector<FieldComponents::Physical::Id>,FieldType::Id,Arithmetics::Id> >::iterator> op3DPhysRange;
                  for(op3DIt = op3D.begin(); op3DIt != op3D.end(); ++op3DIt)
                  {
                     op3DPhysRange = op3DPhys.equal_range(op3DIt->first);
                     for(op3DPhysIt = op3DPhysRange.first; op3DPhysIt != op3DPhysRange.second; ++op3DPhysIt)
                     {
                        ProjectorPhysEdge &rEdge3D = rEdge2D.addEdge(op3DIt->first, op3DIt->second);

                        rEdge3D.setPhysical(std::tr1::get<0>(op3DPhysIt->second), std::tr1::get<1>(op3DPhysIt->second), std::tr1::get<2>(op3DPhysIt->second));
                        DebuggerMacro_showValue("Edge operator: ", 5, op3DIt->first);
                        DebuggerMacro_showValue("Edge weight: ", 5, op3DIt->second);
                        DebuggerMacro_showValue("Physical component: ", 6, rEdge3D.physId());
                        DebuggerMacro_showValue("Field type: ", 6, rEdge3D.fieldId());
                        DebuggerMacro_showValue("Arithmetic: ", 6, rEdge3D.arithId());
                     }
                  }
               }
            }
         }
      }

      // Debugging info 
      DebuggerMacro_leave("Projector trees done", 1);
   }

}
}
