/** \file SpatialSimulationBase.cpp
 *  \brief Source of the building block for the implementation of a simulation
 */

// Debug includes
//
#include "Debug/PrepMacros/ProfilerMacro.h"
#include "Debug/PrepMacros/StorageProfilerMacro.h"

// Configuration includes
//
#include "Base/PrepMacros/SmartPointerMacro.h"
#include "Base/PrepMacros/FrameworkMacro.h"
#include "Simulation/PrepMacros/SpatialSchemeMacro.h"
#include "Simulation/PrepMacros/EquationParametersMacro.h"

// System includes
//

// External includes
//

// Class include
//
#include "Simulation/System/SpatialSimulationBase.hpp"

// Project includes
//
#include "Base/Resolutions/Resolution.hpp"
#include "Base/IO/Ascii/FormatToolbox.hpp"
#include "Base/LoadSplitter/LoadSplitter.hpp"
#include "Base/IO/Control/ConfigParts/PhysicalPart.hpp"
#include "Simulation/IO/Hdf5/StateFileReader.hpp"
#include "Base/LoadSplitter/Algorithms/SplittingDescription.hpp"


namespace EPMPhoenix {

   SpatialSimulationBase::SpatialSimulationBase()
      : SystemBase()
   {
   }

   void SpatialSimulationBase::initSpatial(const int nCpu, const ArrayI& dim)
   {
      // Initialise the base system
      this->initSystem();

      // Initialise the workflow
      FrameworkMacro::setup(nCpu);

      // Create the load splitter
      LoadSplitter splitter(FrameworkMacro::id(), nCpu);

      // Initialise the load splitter
      splitter.init<Code::SpatialScheme>(dim);

      // Get best splitting resolution object
      std::pair<SharedResolution, SplittingDescription>  best = splitter.bestSplitting();

      // Store the shared resolution object
      this->mspRes = best.first;

      // Initialise the transform grouper
      TransformGrouperMacro::setGrouper(best.second, this->mspFwdGrouper, this->mspBwdGrouper);
   }

   void SpatialSimulationBase::initVariables()
   {
      // Iterator over info
      std::map<PhysicalNames::Id, std::pair<bool,TriBool> >::const_iterator infoIt;

      // Storage for the variable info
      std::map<PhysicalNames::Id, std::pair<bool,TriBool> > varInfo;
      // storage for equation variable info
      std::map<PhysicalNames::Id, std::pair<bool,TriBool> > eqInfo;

      //
      // Identify the required variables
      //

      // Loop over all scalar equations
      std::vector<SharedScalarEquation>::iterator scalEqIt;
      for(scalEqIt = this->mScalarEquations.begin(); scalEqIt < this->mScalarEquations.end(); scalEqIt++)
      {
         eqInfo = (*scalEqIt)->requirements();
         for(infoIt = eqInfo.begin(); infoIt != eqInfo.end(); infoIt++)
         {
            if(varInfo.count(infoIt->first) == 0)
            {
               varInfo.insert(*infoIt);
            } else
            {
               varInfo.find(infoIt->first)->second.second += infoIt->second.second;
            }
         }
      }

      // Loop over all vector equations
      std::vector<SharedVectorEquation>::iterator vectEqIt;
      for(vectEqIt = this->mVectorEquations.begin(); vectEqIt < this->mVectorEquations.end(); vectEqIt++)
      {
         eqInfo = (*vectEqIt)->requirements();
         for(infoIt = eqInfo.begin(); infoIt != eqInfo.end(); infoIt++)
         {
            if(varInfo.count(infoIt->first) == 0)
            {
               varInfo.insert(*infoIt);
            } else
            {
               varInfo.find(infoIt->first)->second.second += infoIt->second.second;
            }
         }
      }

      // 
      // Create the required variables
      //

      // Initialise variables
      for(infoIt = varInfo.begin(); infoIt != varInfo.end(); infoIt++)
      {
         // Check if variable is required
         if(infoIt->second.second(0))
         {
            // Separate scalar and vector fields
            if(infoIt->second.first)
            {
               // Create the shared scalar variable
               this->mScalarVariables.insert(std::make_pair(infoIt->first, Code::SharedScalarVariable(new Code::ScalarVariable(this->spRes()))));
            } else
            {
               // Create the shared vector variable
               this->mVectorVariables.insert(std::make_pair(infoIt->first, Code::SharedVectorVariable(new Code::VectorVariable(this->spRes()))));
            }

            // Initialise the rtp values if required
            if(infoIt->second.second(1))
            {
               // Separate scalar and vector fields
               if(infoIt->second.first)
               {
                  this->mScalarVariables.at(infoIt->first)->initialisePhysical();
               } else
               {
                  this->mVectorVariables.at(infoIt->first)->initialisePhysical();
               }
            }

            // Initialise the differential values if required (gradient or curl)
            if(infoIt->second.second(2))
            {
               // Separate scalar and vector fields
               if(infoIt->second.first)
               {
                  this->mScalarVariables.at(infoIt->first)->initialisePhysicalDiff();
               } else
               {
                  this->mVectorVariables.at(infoIt->first)->initialisePhysicalDiff();
               }
            }

            // Separate scalar and vector fields
            if(infoIt->second.first)
            {
               // Initialise to zero
               this->mScalarVariables.at(infoIt->first)->initialiseZeros();

               #ifdef EPMPHOENIX_STORAGEPROFILE
                  StorageProfilerMacro_update(StorageProfiler::VARIABLES, this->mScalarVariables.at(infoIt->first)->requiredStorage());
               #endif // EPMPHOENIX_STORAGEPROFILE
            } else
            {
               // Initialise to zero
               this->mVectorVariables.at(infoIt->first)->initialiseZeros();

               #ifdef EPMPHOENIX_STORAGEPROFILE
                  StorageProfilerMacro_update(StorageProfiler::VARIABLES, this->mVectorVariables.at(infoIt->first)->requiredStorage());
               #endif // EPMPHOENIX_STORAGEPROFILE
            }
         }
      }
      
      // Initialise the transform coordinator
      this->transformCoordinator().initTransforms(*Code::SpatialScheme::spSetup1D(this->mspRes), *Code::SpatialScheme::spSetup2D(this->mspRes), *Code::SpatialScheme::spSetup3D(this->mspRes), varInfo);

      // Initialise the communicator
      this->transformCoordinator().initCommunicator(this->mspRes);

      // Get the buffer pack sizes
      ArrayI packs1DFwd = this->mspFwdGrouper->packs1D(varInfo);
      ArrayI packs2DFwd = this->mspFwdGrouper->packs2D(varInfo);
      ArrayI packs1DBwd = this->mspBwdGrouper->packs1D(varInfo);
      ArrayI packs2DBwd = this->mspBwdGrouper->packs2D(varInfo);

      // Initialise the converters
      this->transformCoordinator().communicator().initConverter(this->mspRes, packs1DFwd, packs1DBwd, packs2DFwd, packs2DBwd, this->mspFwdGrouper->split);
   }

   void SpatialSimulationBase::setupEquations()
   {
      // Loop over all scalar variables
      std::map<PhysicalNames::Id, Code::SharedScalarVariable>::iterator scalIt;
      for(scalIt = this->mScalarVariables.begin(); scalIt != this->mScalarVariables.end(); scalIt++)
      {
         // Loop over scalar equations
         std::vector<SharedScalarEquation>::iterator scalEqIt;
         for(scalEqIt = this->mScalarEquations.begin(); scalEqIt < this->mScalarEquations.end(); scalEqIt++)
         {
            // Set scalar variable as unknown scalar field
            if((*scalEqIt)->name() == scalIt->first)
            {
               (*scalEqIt)->setUnknown(this->mScalarVariables.at(scalIt->first));

               // Finish initialisation of equation
               (*scalEqIt)->init();
            }

            // Set scalar variable as additional scalar field
            if((*scalEqIt)->requirements(scalIt->first)(1) || (*scalEqIt)->requirements(scalIt->first)(2))
            {
               (*scalEqIt)->setField(scalIt->first, this->mScalarVariables.at(scalIt->first));
            }
         }

         // Loop over vector equations
         std::vector<SharedVectorEquation>::iterator vectEqIt;
         for(vectEqIt = this->mVectorEquations.begin(); vectEqIt < this->mVectorEquations.end(); vectEqIt++)
         {
            // Set scalar variable as additional scalar field
            if((*vectEqIt)->requirements(scalIt->first)(1) || (*vectEqIt)->requirements(scalIt->first)(2))
            {
               (*vectEqIt)->setField(scalIt->first, this->mScalarVariables.at(scalIt->first));
            }
         }
      }

      // Loop over all vector variables
      std::map<PhysicalNames::Id, Code::SharedVectorVariable>::iterator vectIt;
      for(vectIt = this->mVectorVariables.begin(); vectIt != this->mVectorVariables.end(); vectIt++)
      {
         // Loop over scalar equations
         std::vector<SharedScalarEquation>::iterator scalEqIt;
         for(scalEqIt = this->mScalarEquations.begin(); scalEqIt < this->mScalarEquations.end(); scalEqIt++)
         {
            // Set vector variable as additional vector field
            if((*scalEqIt)->requirements(vectIt->first)(1) || (*scalEqIt)->requirements(vectIt->first)(2))
            {
               (*scalEqIt)->setField(vectIt->first, this->mVectorVariables.at(vectIt->first));
            }
         }

         // Loop over vector equations
         std::vector<SharedVectorEquation>::iterator vectEqIt;
         for(vectEqIt = this->mVectorEquations.begin(); vectEqIt < this->mVectorEquations.end(); vectEqIt++)
         {
            // Set vector variable as unknown vector field
            if((*vectEqIt)->name() == vectIt->first)
            {
               (*vectEqIt)->setUnknown(this->mVectorVariables.at(vectIt->first));

               // Finish initialisation of equation
               (*vectEqIt)->init();
            }

            // Set vector variable as additional vector field
            if((*vectEqIt)->requirements(vectIt->first)(1) || (*vectEqIt)->requirements(vectIt->first)(2))
            {
               (*vectEqIt)->setField(vectIt->first, this->mVectorVariables.at(vectIt->first));
            }
         }
      }
   }

   void SpatialSimulationBase::computeNonlinear()
   {
      // Compute backward transform
      ProfilerMacro_start(ProfilerMacro::BWDTRANSFORM);
      this->mspBwdGrouper->transform(this->mScalarVariables, this->mVectorVariables, this->transformCoordinator());
      ProfilerMacro_stop(ProfilerMacro::BWDTRANSFORM);

      // compute nonlinear interaction and forward transform
      this->mspFwdGrouper->transform(this->mScalarEquations, this->mVectorEquations, this->transformCoordinator());
   }

}
