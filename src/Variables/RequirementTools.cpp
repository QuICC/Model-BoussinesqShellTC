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
#include "TransformConfigurators/TransformStepsMacro.h"
#include "TransformConfigurators/ProjectorBranch.hpp"
#include "TransformConfigurators/ProjectorTree.hpp"
#include "TransformConfigurators/ProjectorTreeTools.hpp"

namespace GeoMHDiSCC {

   RequirementTools::RequirementTools()
   {
   }

   RequirementTools::~RequirementTools()
   {
   }

   void RequirementTools::initVariables(std::vector<Transform::ProjectorTree>& projectorTree, std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>& rScalarVars, std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>& rVectorVars, const std::vector<Equations::SharedIScalarEquation>& scalarEqs, const std::vector<Equations::SharedIVectorEquation>& vectorEqs, SharedResolution spRes)
   {
      VariableRequirement varInfo;

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
      // Create the required variables and corresponding transform branches
      //
      std::vector<Transform::ProjectorBranch> tmpBranches;
      std::map<PhysicalNames::Id, std::vector<Transform::ProjectorBranch> > branches;

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
               rScalarVars.find(infoIt->first)->second->initSpectral(infoIt->second.spectralIds());
            } else
            {
               // Create the shared vector variable
               rVectorVars.insert(std::make_pair(infoIt->first, Datatypes::SharedVectorVariableType(new Datatypes::VectorVariableType(spRes))));
               rVectorVars.find(infoIt->first)->second->initSpectral(infoIt->second.spectralIds());
            }

            // Initialise transform branch
            branches.insert(std::make_pair(infoIt->first, std::vector<Transform::ProjectorBranch>()));

            // Initialise the physical values if required
            if(infoIt->second.needPhysical())
            {
               // Separate scalar and vector fields
               if(infoIt->second.isScalar())
               {
                  rScalarVars.at(infoIt->first)->initPhysical(infoIt->second.mapPhysicalComps());
                  tmpBranches = Transform::TransformSteps::backwardScalar(infoIt->second.mapPhysicalComps());
                  branches.find(infoIt->first)->second.insert(branches.find(infoIt->first)->second.end(),tmpBranches.begin(), tmpBranches.end());
               } else
               {
                  rVectorVars.at(infoIt->first)->initPhysical(infoIt->second.mapPhysicalComps());
                  tmpBranches = Transform::TransformSteps::backwardVector(infoIt->second.mapPhysicalComps());
                  branches.find(infoIt->first)->second.insert(branches.find(infoIt->first)->second.end(),tmpBranches.begin(), tmpBranches.end());
               }
            }

            // Initialise the physical gradient values if required
            if(infoIt->second.needPhysicalGradient())
            {
               // Separate scalar and vector fields
               if(infoIt->second.isScalar())
               {
                  rScalarVars.at(infoIt->first)->initPhysicalGradient(FieldComponents::Spectral::SCALAR, infoIt->second.mapGradientComps(FieldComponents::Spectral::SCALAR));
                  tmpBranches = Transform::TransformSteps::backwardGradient(infoIt->second.mapGradientComps(FieldComponents::Spectral::SCALAR));
                  branches.find(infoIt->first)->second.insert(branches.find(infoIt->first)->second.end(),tmpBranches.begin(), tmpBranches.end());
               } else
               {
                  std::vector<FieldComponents::Spectral::Id>::const_iterator it;
                  for(it = infoIt->second.spectralIds().begin(); it != infoIt->second.spectralIds().end(); ++it)
                  {
                     rVectorVars.at(infoIt->first)->initPhysicalGradient(*it, infoIt->second.mapGradientComps(*it));
                     tmpBranches = Transform::TransformSteps::backwardVGradient(*it, infoIt->second.mapGradientComps(*it));
                     branches.find(infoIt->first)->second.insert(branches.find(infoIt->first)->second.end(),tmpBranches.begin(), tmpBranches.end());
                  }
               }
            }

            // Initialise the physical curl values if required
            if(infoIt->second.needPhysicalCurl())
            {
               // Separate scalar and vector fields
               if(infoIt->second.isScalar())
               {
                  throw Exception("Can't initialise curl on scalar field!");
               } else
               {
                  rVectorVars.at(infoIt->first)->initPhysicalCurl(infoIt->second.mapCurlComps());
                  tmpBranches = Transform::TransformSteps::backwardCurl(infoIt->second.mapCurlComps());
                  branches.find(infoIt->first)->second.insert(branches.find(infoIt->first)->second.end(),tmpBranches.begin(), tmpBranches.end());
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

      // Create the projector tree(s)
      Transform::ProjectorTreeTools::generateTrees(projectorTree, branches);
   }

   void RequirementTools::mapEquationVariables(std::vector<Transform::IntegratorTree>& integratorTree, std::vector<Equations::SharedIScalarEquation>& rScalarEqs, std::vector<Equations::SharedIVectorEquation>& rVectorEqs, const std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>& scalarVars, const std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>& vectorVars)
   {
      std::vector<Transform::IntegratorBranch> tmpBranches;
      std::map<PhysicalNames::Id, std::vector<Transform::IntegratorBranch> > branches;

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
                  // Initialise transform branch
                  branches.insert(std::make_pair(scalIt->first, std::vector<Transform::IntegratorBranch>()));

                  // Create scalar forward transform
                  tmpBranches = Transform::TransformSteps::forwardScalar();
                  branches.find(scalIt->first)->second.insert(branches.find(scalIt->first)->second.end(),tmpBranches.begin(), tmpBranches.end());
               }
            }

            // Set scalar variable as additional scalar field
            if((*scalEqIt)->requirements(scalIt->first).needPhysical() || (*scalEqIt)->requirements(scalIt->first).needPhysicalGradient())
            {
               (*scalEqIt)->setField(scalIt->first, scalarVars.at(scalIt->first));
            }
         }

         // Loop over vector equations
         std::vector<Equations::SharedIVectorEquation>::iterator vectEqIt;
         for(vectEqIt = rVectorEqs.begin(); vectEqIt < rVectorEqs.end(); vectEqIt++)
         {
            // Set scalar variable as additional scalar field
            if((*vectEqIt)->requirements(scalIt->first).needPhysical() || (*vectEqIt)->requirements(scalIt->first).needPhysicalGradient())
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
            if((*scalEqIt)->requirements(vectIt->first).needPhysical() || (*scalEqIt)->requirements(vectIt->first).needPhysicalGradient() || (*scalEqIt)->requirements(vectIt->first).needPhysicalCurl())
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
                  // Initialise transform branch
                  branches.insert(std::make_pair(vectIt->first, std::vector<Transform::IntegratorBranch>()));

                  // Create scalar forward transform
                  tmpBranches = Transform::TransformSteps::forwardScalar();
                  branches.find(vectIt->first)->second.insert(branches.find(vectIt->first)->second.end(),tmpBranches.begin(), tmpBranches.end());
               }
            }

            // Set vector variable as additional vector field
            if((*vectEqIt)->requirements(vectIt->first).needPhysical() || (*vectEqIt)->requirements(vectIt->first).needPhysicalGradient() || (*vectEqIt)->requirements(vectIt->first).needPhysicalCurl())
            {
               (*vectEqIt)->setField(vectIt->first, vectorVars.at(vectIt->first));
            }
         }
      }

      // Create the integrator tree(s)
      Transform::IntegratorTreeTools::generateTrees(integratorTree, branches);
   }
}
