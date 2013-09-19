/** 
 * @file RequirementTools.cpp
 * @brief Source of the requirement tools to work with variables and equations
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "Variables/RequirementTools.hpp"

// Project includes
//

namespace GeoMHDiSCC {

   RequirementTools::RequirementTools()
   {
   }

   RequirementTools::~RequirementTools()
   {
   }

   void RequirementTools::initVariables(VariableRequirement& varInfo, std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>& rScalarVars, std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>& rVectorVars, const std::vector<Equations::SharedIScalarEquation>& scalarEqs, const std::vector<Equations::SharedIVectorEquation>& vectorEqs, SharedResolution spRes)
   {
      // Iterator over info
      VariableRequirement::const_iterator infoIt;

      //
      // Identify the required variables
      //

      // Loop over all scalar equations
      std::vector<Equations::SharedIScalarEquation>::const_iterator scalEqIt;
      for(scalEqIt = scalarEqs.begin(); scalEqIt < scalarEqs.end(); scalEqIt++)
      {
         varInfo.merge((*scalEqIt)->requirements());
      }

      // Loop over all vector equations
      std::vector<Equations::SharedIVectorEquation>::const_iterator vectEqIt;
      for(vectEqIt = vectorEqs.begin(); vectEqIt < vectorEqs.end(); vectEqIt++)
      {
         varInfo.merge((*vectEqIt)->requirements());
      }

      // 
      // Create the required variables
      //

      // Initialise variables
      for(infoIt = varInfo.begin(); infoIt != varInfo.end(); infoIt++)
      {
         // Check if spectral variable is required
         if(infoIt->second.needSpectral())
         {
            // Separate scalar and vector fields
            if(infoIt->second.isScalar())
            {
               // Create the shared scalar variable
               rScalarVars.insert(std::make_pair(infoIt->first, Datatypes::SharedScalarVariableType(new Datatypes::ScalarVariableType(spRes))));
            } else
            {
               // Create the shared vector variable
               rVectorVars.insert(std::make_pair(infoIt->first, Datatypes::SharedVectorVariableType(new Datatypes::VectorVariableType(spRes))));
            }

            // Initialise the physical values if required
            if(infoIt->second.needPhysical())
            {
               // Separate scalar and vector fields
               if(infoIt->second.isScalar())
               {
                  rScalarVars.at(infoIt->first)->initPhysical();
               } else
               {
                  rVectorVars.at(infoIt->first)->initPhysical();
               }
            }

            // Initialise the physical differential values if required (gradient or curl)
            if(infoIt->second.needPhysicalDiff())
            {
               // Separate scalar and vector fields
               if(infoIt->second.isScalar())
               {
                  rScalarVars.at(infoIt->first)->initPhysicalDiff();
               } else
               {
                  rVectorVars.at(infoIt->first)->initPhysicalDiff();
               }
            }

            // Separate scalar and vector fields
            if(infoIt->second.isScalar())
            {
               // Initialise to zero
               rScalarVars.at(infoIt->first)->setZeros();

               #ifdef GEOMHDISCC_STORAGEPROFILE
                  StorageProfilerMacro_update(Debug::StorageProfiler::VARIABLES, rScalarVars.at(infoIt->first)->requiredStorage());
               #endif // GEOMHDISCC_STORAGEPROFILE
            } else
            {
               // Initialise to zero
               rVectorVars.at(infoIt->first)->setZeros();

               #ifdef GEOMHDISCC_STORAGEPROFILE
                  StorageProfilerMacro_update(Debug::StorageProfiler::VARIABLES, rVectorVars.at(infoIt->first)->requiredStorage());
               #endif // GEOMHDISCC_STORAGEPROFILE
            }
         }
      }
   }

   void RequirementTools::mapEquationVariables(std::set<PhysicalNames::Id>& nonInfo, std::vector<Equations::SharedIScalarEquation>& rScalarEqs, std::vector<Equations::SharedIVectorEquation>& rVectorEqs, const std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>& scalarVars, const std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>& vectorVars)
   {
      // Loop over all scalar variables
      std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>::const_iterator scalIt;
      for(scalIt = scalarVars.begin(); scalIt != scalarVars.end(); scalIt++)
      {
         // Loop over scalar equations
         std::vector<Equations::SharedIScalarEquation>::iterator scalEqIt;
         for(scalEqIt = rScalarEqs.begin(); scalEqIt < rScalarEqs.end(); scalEqIt++)
         {
            // Set scalar variable as unknown scalar field
            if((*scalEqIt)->name() == scalIt->first)
            {
               (*scalEqIt)->setUnknown(scalarVars.at(scalIt->first));

               // Finish initialisation of equation
               (*scalEqIt)->init();

               // Check for nonlinear requirements
               if((*scalEqIt)->couplingInfo(FieldComponents::Spectral::SCALAR).hasNonlinear())
               {
                  nonInfo.insert((*scalEqIt)->name());
               }
            }

            // Set scalar variable as additional scalar field
            if((*scalEqIt)->requirements(scalIt->first).needPhysical() || (*scalEqIt)->requirements(scalIt->first).needPhysicalDiff())
            {
               (*scalEqIt)->setField(scalIt->first, scalarVars.at(scalIt->first));
            }
         }

         // Loop over vector equations
         std::vector<Equations::SharedIVectorEquation>::iterator vectEqIt;
         for(vectEqIt = rVectorEqs.begin(); vectEqIt < rVectorEqs.end(); vectEqIt++)
         {
            // Set scalar variable as additional scalar field
            if((*vectEqIt)->requirements(scalIt->first).needPhysical() || (*vectEqIt)->requirements(scalIt->first).needPhysicalDiff())
            {
               (*vectEqIt)->setField(scalIt->first, scalarVars.at(scalIt->first));
            }
         }
      }

      // Loop over all vector variables
      std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>::const_iterator vectIt;
      for(vectIt = vectorVars.begin(); vectIt != vectorVars.end(); vectIt++)
      {
         // Loop over scalar equations
         std::vector<Equations::SharedIScalarEquation>::iterator scalEqIt;
         for(scalEqIt = rScalarEqs.begin(); scalEqIt < rScalarEqs.end(); scalEqIt++)
         {
            // Set vector variable as additional vector field
            if((*scalEqIt)->requirements(vectIt->first).needPhysical() || (*scalEqIt)->requirements(vectIt->first).needPhysicalDiff())
            {
               (*scalEqIt)->setField(vectIt->first, vectorVars.at(vectIt->first));
            }
         }

         // Loop over vector equations
         std::vector<Equations::SharedIVectorEquation>::iterator vectEqIt;
         for(vectEqIt = rVectorEqs.begin(); vectEqIt < rVectorEqs.end(); vectEqIt++)
         {
            // Set vector variable as unknown vector field
            if((*vectEqIt)->name() == vectIt->first)
            {
               (*vectEqIt)->setUnknown(vectorVars.at(vectIt->first));

               // Finish initialisation of equation
               (*vectEqIt)->init();

               // Check for nonlinear requirements
               if((*vectEqIt)->couplingInfo(FieldComponents::Spectral::ONE).hasNonlinear())
               {
                  nonInfo.insert((*vectEqIt)->name());
               }
            }

            // Set vector variable as additional vector field
            if((*vectEqIt)->requirements(vectIt->first).needPhysical() || (*vectEqIt)->requirements(vectIt->first).needPhysicalDiff())
            {
               (*vectEqIt)->setField(vectIt->first, vectorVars.at(vectIt->first));
            }
         }
      }
   }
}
