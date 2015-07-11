/** 
 * @file IntegratorTree2DTools.cpp
 * @brief Source of the implementation of tools to work with integrator transform trees for 2D space
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
#include "TransformConfigurators/IntegratorTree2DTools.hpp"

// Intgect includes
//

namespace GeoMHDiSCC {

namespace Transform {

   void IntegratorTree2DTools::generateTrees(std::vector<IntegratorTree2D>& rTrees, const std::map<PhysicalNames::Id, std::vector<IntegratorBranch2D> >& branches)
   {
      // Debugging info 
      DebuggerMacro_enter("Generating integrator trees", 1);

      // Loop over all physical fields
      std::map<PhysicalNames::Id, std::vector<IntegratorBranch2D> >::const_iterator nameIt;
      for(nameIt = branches.begin(); nameIt != branches.end(); ++nameIt)
      {
         // Loop over branches
         std::vector<IntegratorBranch2D>::const_iterator branchIt;
         std::set<FieldComponents::Physical::Id> components;
         for(branchIt = nameIt->second.begin(); branchIt != nameIt->second.end(); ++branchIt)
         {
            // Generate list of unique components
            components.insert(branchIt->physId());
         }

         // Initialise component trees
         std::set<FieldComponents::Physical::Id>::const_iterator compIt;
         for(compIt = components.begin(); compIt != components.end(); ++compIt)
         {
            // Initilize tree
            rTrees.push_back(IntegratorTree2D(nameIt->first, *compIt));

            DebuggerMacro_showValue("Field: ", 2, nameIt->first);
            DebuggerMacro_showValue("Component: ", 2, *compIt);

            // Extract unique physical operators
            std::map<IntgPhysId, int> opPhys;
            std::pair<std::map<IntgPhysId,int>::iterator,bool> opPhysPairIt;
            for(branchIt = nameIt->second.begin(); branchIt != nameIt->second.end(); ++branchIt)
            {
               if(branchIt->physId() == *compIt)
               {
                  opPhysPairIt = opPhys.insert(std::make_pair(branchIt->intgPhysId(), 0));
                  opPhysPairIt.first->second += 1;
               }
            }

            // Create physical edges
            std::map<IntgPhysId,int>::const_iterator opPhysIt;
            for(opPhysIt = opPhys.begin(); opPhysIt != opPhys.end(); ++opPhysIt)
            {
               IntegratorPhysEdge &rEdgePhys = rTrees.back().addEdge(opPhysIt->first, opPhysIt->second);
               DebuggerMacro_showValue("Edge operator: ", 3, opPhysIt->first);
               DebuggerMacro_showValue("Edge weight: ", 3, opPhysIt->second);

               // Extract unique spectral Operators
               std::map<IntgSpecId, int> opSpec;
               std::map<IntgSpecId, std::tr1::tuple<FieldComponents::Spectral::Id,FieldType::Id,Arithmetics::Id> > op1DSpec;
               std::pair<std::map<IntgSpecId,int>::iterator,bool> opSpecPairIt;
               for(branchIt = nameIt->second.begin(); branchIt != nameIt->second.end(); ++branchIt)
               {
                  if(branchIt->physId() == *compIt && branchIt->intgPhysId() == opPhysIt->first)
                  {
                     opSpecPairIt = opSpec.insert(std::make_pair(branchIt->intgSpecId(), 0));
                     opSpecPairIt.first->second += 1;
                     op1DSpec.insert(std::make_pair(branchIt->intgSpecId(), std::tr1::make_tuple(branchIt->specId(), branchIt->fieldId(), branchIt->arithId())));
                  }
               }

               // Create spectral edges
               std::map<IntgSpecId,int>::const_iterator opSpecIt;
               for(opSpecIt = opSpec.begin(); opSpecIt != opSpec.end(); ++opSpecIt)
               {
                  IntegratorSpecEdge &rEdgeSpec = rEdgePhys.addEdge(opSpecIt->first, opSpecIt->second);

                  rEdgeSpec.setSpectral(std::tr1::get<0>(op1DSpec.find(opSpecIt->first)->second), std::tr1::get<1>(op1DSpec.find(opSpecIt->first)->second), std::tr1::get<2>(op1DSpec.find(opSpecIt->first)->second));
                  DebuggerMacro_showValue("Edge operator: ", 4, opSpecIt->first);
                  DebuggerMacro_showValue("Edge weight: ", 4, opSpecIt->second);
                  DebuggerMacro_showValue("Spectral component: ", 4, rEdgeSpec.specId());
                  DebuggerMacro_showValue("Field type: ", 4, rEdgeSpec.fieldId());
                  DebuggerMacro_showValue("Arithmetic: ", 4, rEdgeSpec.arithId());
               }
            }
         }
      }

      // Debugging info 
      DebuggerMacro_leave("Integrator trees done", 1);
   }

}
}
