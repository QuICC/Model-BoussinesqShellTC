/** 
 * @file ProjectorTree2DTools.cpp
 * @brief Source of the implementation of tools to work with projector transform trees for 2D space
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
#include "TransformConfigurators/ProjectorTree2DTools.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

   void ProjectorTree2DTools::generateTrees(std::vector<ProjectorTree2D>& rTrees, const std::map<PhysicalNames::Id, std::vector<ProjectorBranch2D> >& branches)
   {
      // Debugging info 
      DebuggerMacro_enter("Generating projector trees", 1);

      // Loop over all physical fields
      std::map<PhysicalNames::Id, std::vector<ProjectorBranch2D> >::const_iterator nameIt;
      for(nameIt = branches.begin(); nameIt != branches.end(); ++nameIt)
      {
         // Loop over branches
         std::vector<ProjectorBranch2D>::const_iterator branchIt;
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
            rTrees.push_back(ProjectorTree2D(nameIt->first, *compIt));

            DebuggerMacro_showValue("Field: ", 2, nameIt->first);
            DebuggerMacro_showValue("Component: ", 2, *compIt);

            // Extract unique spectral operators
            std::map<ProjSpecId, int> opSpec;
            std::pair<std::map<ProjSpecId,int>::iterator,bool> opSpecPairIt;
            for(branchIt = nameIt->second.begin(); branchIt != nameIt->second.end(); ++branchIt)
            {
               if(branchIt->specId() == *compIt)
               {
                  opSpecPairIt = opSpec.insert(std::make_pair(branchIt->projSpecId(), 0));
                  opSpecPairIt.first->second += 1;
               }
            }

            // Create spectral edges
            std::map<ProjSpecId,int>::const_iterator opSpecIt;
            for(opSpecIt = opSpec.begin(); opSpecIt != opSpec.end(); ++opSpecIt)
            {
               ProjectorSpecEdge &rEdgeSpec = rTrees.back().addEdge(opSpecIt->first, opSpecIt->second);
               DebuggerMacro_showValue("Edge operator: ", 3, opSpecIt->first);
               DebuggerMacro_showValue("Edge weight: ", 3, opSpecIt->second);

               // Extract unique physical Operators
               std::map<ProjPhysId, int> opPhys;
               std::multimap<ProjPhysId, std::tr1::tuple<FieldComponents::Physical::Id,FieldType::Id,Arithmetics::Id> > op3DPhys;
               std::pair<std::map<ProjPhysId,int>::iterator,bool> opPhysPairIt;
               for(branchIt = nameIt->second.begin(); branchIt != nameIt->second.end(); ++branchIt)
               {
                  if(branchIt->specId() == *compIt && branchIt->projSpecId() == opSpecIt->first)
                  {
                     opPhysPairIt = opPhys.insert(std::make_pair(branchIt->projPhysId(), 0));
                     opPhysPairIt.first->second += 1;
                     op3DPhys.insert(std::make_pair(branchIt->projPhysId(), std::tr1::make_tuple(branchIt->physId(), branchIt->fieldId(), branchIt->arithId())));
                  }
               }

               // Create physical edges
               std::map<ProjPhysId,int>::const_iterator opPhysIt;
               std::multimap<ProjPhysId, std::tr1::tuple<FieldComponents::Physical::Id,FieldType::Id,Arithmetics::Id> >::iterator op3DPhysIt;
               std::pair<std::multimap<ProjPhysId, std::tr1::tuple<FieldComponents::Physical::Id,FieldType::Id,Arithmetics::Id> >::iterator,std::multimap<ProjPhysId, std::tr1::tuple<FieldComponents::Physical::Id,FieldType::Id,Arithmetics::Id> >::iterator> op3DPhysRange;
               for(opPhysIt = opPhys.begin(); opPhysIt != opPhys.end(); ++opPhysIt)
               {
                  op3DPhysRange = op3DPhys.equal_range(opPhysIt->first);
                  for(op3DPhysIt = op3DPhysRange.first; op3DPhysIt != op3DPhysRange.second; ++op3DPhysIt)
                  {
                     ProjectorPhysEdge &rEdgePhys = rEdgeSpec.addEdge(opPhysIt->first, opPhysIt->second);

                     rEdgePhys.setPhysical(std::tr1::get<0>(op3DPhysIt->second), std::tr1::get<1>(op3DPhysIt->second), std::tr1::get<2>(op3DPhysIt->second));
                     DebuggerMacro_showValue("Edge operator: ", 4, opPhysIt->first);
                     DebuggerMacro_showValue("Edge weight: ", 4, opPhysIt->second);
                     DebuggerMacro_showValue("Physical component: ", 4, rEdgePhys.physId());
                     DebuggerMacro_showValue("Field type: ", 4, rEdgePhys.fieldId());
                     DebuggerMacro_showValue("Arithmetic: ", 4, rEdgePhys.arithId());
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
