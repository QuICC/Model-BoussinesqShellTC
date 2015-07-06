/** 
 * @file IntegratorTree3DTools.cpp
 * @brief Source of the implementation of tools to work with integrator transform trees for 3D space
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
#include "TransformConfigurators/IntegratorTree3DTools.hpp"

// Intgect includes
//

namespace GeoMHDiSCC {

namespace Transform {

   void IntegratorTree3DTools::generateTrees(std::vector<IntegratorTree3D>& rTrees, const std::map<PhysicalNames::Id, std::vector<IntegratorBranch3D> >& branches)
   {
      // Debugging info 
      DebuggerMacro_enter("Generating integrator trees", 1);

      // Loop over all physical fields
      std::map<PhysicalNames::Id, std::vector<IntegratorBranch3D> >::const_iterator nameIt;
      for(nameIt = branches.begin(); nameIt != branches.end(); ++nameIt)
      {
         // Loop over branches
         std::vector<IntegratorBranch3D>::const_iterator branchIt;
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
            rTrees.push_back(IntegratorTree3D(nameIt->first, *compIt));

            DebuggerMacro_showValue("Field: ", 2, nameIt->first);
            DebuggerMacro_showValue("Component: ", 2, *compIt);

            // Extract unique 3D operators
            std::map<IntgPhysId, int> op3D;
            std::pair<std::map<IntgPhysId,int>::iterator,bool> op3DPairIt;
            for(branchIt = nameIt->second.begin(); branchIt != nameIt->second.end(); ++branchIt)
            {
               if(branchIt->physId() == *compIt)
               {
                  op3DPairIt = op3D.insert(std::make_pair(branchIt->intgPhysId(), 0));
                  op3DPairIt.first->second += 1;
               }
            }

            // Create 3D edges
            std::map<IntgPhysId,int>::const_iterator op3DIt;
            for(op3DIt = op3D.begin(); op3DIt != op3D.end(); ++op3DIt)
            {
               IntegratorPhysEdge &rEdge3D = rTrees.back().addEdge(op3DIt->first, op3DIt->second);
               DebuggerMacro_showValue("Edge operator: ", 3, op3DIt->first);
               DebuggerMacro_showValue("Edge weight: ", 3, op3DIt->second);

               // Extract unique 2D Operators
               std::map<IntgPartId, int> op2D;
               std::pair<std::map<IntgPartId,int>::iterator,bool> op2DPairIt;
               for(branchIt = nameIt->second.begin(); branchIt != nameIt->second.end(); ++branchIt)
               {
                  if(branchIt->physId() == *compIt && branchIt->intgPhysId() == op3DIt->first)
                  {
                     op2DPairIt = op2D.insert(std::make_pair(branchIt->intgPartId(), 0));
                     op2DPairIt.first->second += 1;
                  }
               }

               // Create 2D edges
               std::map<IntgPartId,int>::const_iterator op2DIt;
               for(op2DIt = op2D.begin(); op2DIt != op2D.end(); ++op2DIt)
               {
                  IntegratorPartEdge &rEdge2D = rEdge3D.addEdge(op2DIt->first, op2DIt->second);
                  DebuggerMacro_showValue("Edge operator: ", 4, op2DIt->first);
                  DebuggerMacro_showValue("Edge weight: ", 4, op2DIt->second);

                  // Extract unique 1D operators
                  std::map<IntgSpecId, int> op1D;
                  std::map<IntgSpecId, std::tr1::tuple<FieldComponents::Spectral::Id,FieldType::Id,Arithmetics::Id> > op1DSpec;
                  std::pair<std::map<IntgSpecId,int>::iterator,bool> op1DPairIt;
                  for(branchIt = nameIt->second.begin(); branchIt != nameIt->second.end(); ++branchIt)
                  {
                     if(branchIt->physId() == *compIt && branchIt->intgPhysId() == op3DIt->first && branchIt->intgPartId() == op2DIt->first)
                     {
                        op1DPairIt = op1D.insert(std::make_pair(branchIt->intgSpecId(), 0));
                        op1DPairIt.first->second += 1;
                        op1DSpec.insert(std::make_pair(branchIt->intgSpecId(), std::tr1::make_tuple(branchIt->specId(), branchIt->fieldId(), branchIt->arithId())));
                     }
                  }

                  // Create 1D edges
                  std::map<IntgSpecId,int>::const_iterator op1DIt;
                  for(op1DIt = op1D.begin(); op1DIt != op1D.end(); ++op1DIt)
                  {
                     IntegratorSpecEdge &rEdge1D = rEdge2D.addEdge(op1DIt->first, op1DIt->second);

                     rEdge1D.setSpectral(std::tr1::get<0>(op1DSpec.find(op1DIt->first)->second), std::tr1::get<1>(op1DSpec.find(op1DIt->first)->second), std::tr1::get<2>(op1DSpec.find(op1DIt->first)->second));
                     DebuggerMacro_showValue("Edge operator: ", 5, op1DIt->first);
                     DebuggerMacro_showValue("Edge weight: ", 5, op3DIt->second);
                     DebuggerMacro_showValue("Spectral component: ", 6, rEdge1D.specId());
                     DebuggerMacro_showValue("Field type: ", 6, rEdge1D.fieldId());
                     DebuggerMacro_showValue("Arithmetic: ", 6, rEdge1D.arithId());
                  }
               }
            }
         }
      }

      // Debugging info 
      DebuggerMacro_leave("Integrator trees done", 1);
   }

}
}
