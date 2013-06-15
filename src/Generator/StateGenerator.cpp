/** \file StateGenerator.cpp
 *  \brief Source of the high level state generator
 */

// Debug includes
//

// Configuration includes
//
#include "Debug/DebuggerMacro.h"
#include "Profiler/ProfilerMacro.h"
#include "StorageProfiler/StorageProfilerMacro.h"
#include "Framework/FrameworkMacro.h"

// System includes
//

// External includes
//

// Class include
//
#include "Generator/StateGenerator.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"
#include "IoTools/Formatter.hpp"
#include "SpectralOperators/BoundaryConditions.hpp"

namespace GeoMHDiSCC {

   StateGenerator::StateGenerator()
   {
   }

   StateGenerator::~StateGenerator()
   {
   }

   void StateGenerator::initBase()
   {
      // Debug statement
      DebuggerMacro_enter("initBase",0);

      // Make sure to catch raised exception in initialisation steps
      try{
         // Initialise the IO system
         this->mSimIoCtrl.init();

         // Initialise the equation parameters
         this->mspEqParams->init(this->mSimIoCtrl.configPhysical());

         // Get number CPU from configuration file
         int nCpu = this->mSimIoCtrl.configNCpu();

         // Initialise the workflow
         FrameworkMacro::setup(nCpu);
      }
      catch(Exception &e)
      {
         std::cout << e.what() << std::endl;

         throw -1;
      }

      // Make sure nodes are synchronised after initialisation
      FrameworkMacro::synchronize();

      // Debug statement
      DebuggerMacro_leave("initBase",0);
   }

   void StateGenerator::init(const SharedSimulationBoundary spBcs)
   {
      // Debug statement
      DebuggerMacro_enter("init",0);

      // Storage for the variable info
      VariableRequirement varInfo;

      // Initialise the variables and set general variable requirements
      this->initVariables(varInfo);

      // Storage for the nonlinear requirement info
      std::set<PhysicalNames::Id>   nonInfo;

      // Map variables to the equations and set nonlinear requirements
      this->mapEquationVariables(nonInfo);

      // Initialise the transform coordinator
      this->initTransformCoordinator(varInfo, nonInfo);

      // Initialise the equations (generate operators, etc)
      this->setupEquations(spBcs);

      // Setup output files (ASCII diagnostics, state files, etc)
      this->setupOutput();

      // Cleanup IO control
      this->mSimIoCtrl.cleanup();

      // Print message to signal successful completion of initialisation step
      if(FrameworkMacro::allowsIO())
      {
         IoTools::Formatter::printLine(std::cout, '-');
         IoTools::Formatter::printCentered(std::cout, "State generator initialisation successfull", '*');
         IoTools::Formatter::printLine(std::cout, '-');
         IoTools::Formatter::printNewline(std::cout);
      }

      // Debug statement
      DebuggerMacro_leave("init",0);
   }

   void StateGenerator::run()
   {
      // Debug statement
      DebuggerMacro_enter("run",0);

      // Compute backward transform
      ProfilerMacro_start(ProfilerMacro::BWDTRANSFORM);
      this->mspBwdGrouper->transform(this->mScalarVariables, this->mVectorVariables, this->mTransformCoordinator);
      ProfilerMacro_stop(ProfilerMacro::BWDTRANSFORM);

      // compute nonlinear interaction and forward transform
      this->mspFwdGrouper->transform(this->mScalarEquations, this->mVectorEquations, this->mTransformCoordinator);

      // Synchronise computation nodes
      FrameworkMacro::synchronize();

      // Debug statement
      DebuggerMacro_leave("run",0);
   }

   void StateGenerator::finalize()
   {
      // Debug statement
      DebuggerMacro_enter("finalize",0);

      // Print profiling infos (if required)
      ProfilerMacro_printInfo();

      // Print storage profiling infos (if required)
      StorageProfilerMacro_printInfo();

      // Debug statement
      DebuggerMacro_leave("finalize",0);
   }

   void StateGenerator::setInitialState(IoVariable::SharedIVariableHdf5Reader spInitFile)
   {
      // Loop over all scalars
      std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>::iterator scalIt;
      for(scalIt = this->mScalarVariables.begin(); scalIt != this->mScalarVariables.end(); scalIt++)
      {
         spInitFile->addScalar((*scalIt));
      }

      // Loop over all vector variables
      std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>::iterator vectIt;
      for(vectIt = this->mVectorVariables.begin(); vectIt != this->mVectorVariables.end(); vectIt++)
      {
         spInitFile->addVector((*vectIt));
      }

      // Initialise file
      spInitFile->init();

      // Read in data
      spInitFile->read();

      // Finalise file
      spInitFile->finalize();
   }

   void StateGenerator::writeOutput()
   {
      // Debug statement
      DebuggerMacro_enter("postRun",1);

      // Get time and timestep from configuration
      Array tstep = this->mSimIoCtrl.configTimestepping();

      // Set to zero if value is negative
      if(tstep(0) < 0)
      {
         tstep(0) = 0;
      }
      if(tstep(1) < 0)
      {
         tstep(1) = 0;
      }

      // Write final ASCII output
      this->mSimIoCtrl.writeAscii(tstep(0), tstep(1));

      // Write final state file
      this->mSimIoCtrl.writeHdf5(tstep(0), tstep(1));

      // Synchronise all nodes of simulation
      FrameworkMacro::synchronize();

      // Debug statement
      DebuggerMacro_leave("postRun",1);
   }

   void StateGenerator::initVariables(VariableRequirement& varInfo)
   {
      // Iterator over info
      VariableRequirement::const_iterator infoIt;

      //
      // Identify the required variables
      //

      // Loop over all scalar equations
      std::vector<Equations::SharedIScalarEquation>::iterator scalEqIt;
      for(scalEqIt = this->mScalarEquations.begin(); scalEqIt < this->mScalarEquations.end(); scalEqIt++)
      {
         varInfo.merge((*scalEqIt)->requirements());
      }

      // Loop over all vector equations
      std::vector<Equations::SharedIVectorEquation>::iterator vectEqIt;
      for(vectEqIt = this->mVectorEquations.begin(); vectEqIt < this->mVectorEquations.end(); vectEqIt++)
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
               this->mScalarVariables.insert(std::make_pair(infoIt->first, Datatypes::SharedScalarVariableType(new Datatypes::ScalarVariableType(this->mspRes))));
            } else
            {
               // Create the shared vector variable
               this->mVectorVariables.insert(std::make_pair(infoIt->first, Datatypes::SharedVectorVariableType(new Datatypes::VectorVariableType(this->mspRes))));
            }

            // Initialise the physical values if required
            if(infoIt->second.needPhysical())
            {
               // Separate scalar and vector fields
               if(infoIt->second.isScalar())
               {
                  this->mScalarVariables.at(infoIt->first)->initPhysical();
               } else
               {
                  this->mVectorVariables.at(infoIt->first)->initPhysical();
               }
            }

            // Initialise the physical differential values if required (gradient or curl)
            if(infoIt->second.needPhysicalDiff())
            {
               // Separate scalar and vector fields
               if(infoIt->second.isScalar())
               {
                  this->mScalarVariables.at(infoIt->first)->initPhysicalDiff();
               } else
               {
                  this->mVectorVariables.at(infoIt->first)->initPhysicalDiff();
               }
            }

            // Separate scalar and vector fields
            if(infoIt->second.isScalar())
            {
               // Initialise to zero
               this->mScalarVariables.at(infoIt->first)->setZeros();

               #ifdef EPMPHOENIX_STORAGEPROFILE
                  StorageProfilerMacro_update(StorageProfiler::VARIABLES, this->mScalarVariables.at(infoIt->first)->requiredStorage());
               #endif // EPMPHOENIX_STORAGEPROFILE
            } else
            {
               // Initialise to zero
               this->mVectorVariables.at(infoIt->first)->setZeros();

               #ifdef EPMPHOENIX_STORAGEPROFILE
                  StorageProfilerMacro_update(StorageProfiler::VARIABLES, this->mVectorVariables.at(infoIt->first)->requiredStorage());
               #endif // EPMPHOENIX_STORAGEPROFILE
            }
         }
      }
   }

   void StateGenerator::mapEquationVariables(std::set<PhysicalNames::Id>& nonInfo)
   {
      // Loop over all scalar variables
      std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>::iterator scalIt;
      for(scalIt = this->mScalarVariables.begin(); scalIt != this->mScalarVariables.end(); scalIt++)
      {
         // Loop over scalar equations
         std::vector<Equations::SharedIScalarEquation>::iterator scalEqIt;
         for(scalEqIt = this->mScalarEquations.begin(); scalEqIt < this->mScalarEquations.end(); scalEqIt++)
         {
            // Set scalar variable as unknown scalar field
            if((*scalEqIt)->name() == scalIt->first)
            {
               (*scalEqIt)->setUnknown(this->mScalarVariables.at(scalIt->first));

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
               (*scalEqIt)->setField(scalIt->first, this->mScalarVariables.at(scalIt->first));
            }
         }

         // Loop over vector equations
         std::vector<Equations::SharedIVectorEquation>::iterator vectEqIt;
         for(vectEqIt = this->mVectorEquations.begin(); vectEqIt < this->mVectorEquations.end(); vectEqIt++)
         {
            // Set scalar variable as additional scalar field
            if((*vectEqIt)->requirements(scalIt->first).needPhysical() || (*vectEqIt)->requirements(scalIt->first).needPhysicalDiff())
            {
               (*vectEqIt)->setField(scalIt->first, this->mScalarVariables.at(scalIt->first));
            }
         }
      }

      // Loop over all vector variables
      std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>::iterator vectIt;
      for(vectIt = this->mVectorVariables.begin(); vectIt != this->mVectorVariables.end(); vectIt++)
      {
         // Loop over scalar equations
         std::vector<Equations::SharedIScalarEquation>::iterator scalEqIt;
         for(scalEqIt = this->mScalarEquations.begin(); scalEqIt < this->mScalarEquations.end(); scalEqIt++)
         {
            // Set vector variable as additional vector field
            if((*scalEqIt)->requirements(vectIt->first).needPhysical() || (*scalEqIt)->requirements(vectIt->first).needPhysicalDiff())
            {
               (*scalEqIt)->setField(vectIt->first, this->mVectorVariables.at(vectIt->first));
            }
         }

         // Loop over vector equations
         std::vector<Equations::SharedIVectorEquation>::iterator vectEqIt;
         for(vectEqIt = this->mVectorEquations.begin(); vectEqIt < this->mVectorEquations.end(); vectEqIt++)
         {
            // Set vector variable as unknown vector field
            if((*vectEqIt)->name() == vectIt->first)
            {
               (*vectEqIt)->setUnknown(this->mVectorVariables.at(vectIt->first));

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
               (*vectEqIt)->setField(vectIt->first, this->mVectorVariables.at(vectIt->first));
            }
         }
      }
   }
      
   void StateGenerator::initTransformCoordinator(const VariableRequirement& varInfo, const std::set<PhysicalNames::Id>& nonInfo)
   {
      // Initialise the transform coordinator
      this->mTransformCoordinator.initTransforms(this->mspRes, varInfo);

      // Initialise the communicator
      this->mTransformCoordinator.initCommunicator(this->mspRes);

      // Get the buffer pack sizes
      ArrayI packs1DFwd = this->mspFwdGrouper->packs1D(varInfo, nonInfo);
      ArrayI packs2DFwd = this->mspFwdGrouper->packs2D(varInfo, nonInfo);
      ArrayI packs1DBwd = this->mspBwdGrouper->packs1D(varInfo);
      ArrayI packs2DBwd = this->mspBwdGrouper->packs2D(varInfo);

      // Initialise the converters
      this->mTransformCoordinator.communicator().initConverter(this->mspRes, packs1DFwd, packs1DBwd, packs2DFwd, packs2DBwd, this->mspFwdGrouper->split);
   }

   void StateGenerator::setupEquations(const SharedSimulationBoundary spBcs)
   {
      // Create iterators over scalar equations
      std::vector<Equations::SharedIScalarEquation>::iterator scalEqIt;
      // Create iterators over vector equations
      std::vector<Equations::SharedIVectorEquation>::iterator vectEqIt;

      // Loop over all scalar equations
      for(scalEqIt = this->mScalarEquations.begin(); scalEqIt < this->mScalarEquations.end(); ++scalEqIt)
      {
         (*scalEqIt)->initSpectralMatrices(spBcs);
      }

      // Loop over all vector equations
      for(vectEqIt = this->mVectorEquations.begin(); vectEqIt < this->mVectorEquations.end(); ++vectEqIt)
      {
         (*vectEqIt)->initSpectralMatrices(spBcs);
      }
   }

   void StateGenerator::setupOutput()
   {
      // Loop over all files added to the simulation control
      SimulationIoControl::hdf5_iterator  fIt;
      for(fIt = this->mSimIoCtrl.beginHdf5(); fIt != this->mSimIoCtrl.endHdf5(); ++fIt)
      {
         // Loop over all scalars
         std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>::iterator scalIt;
         for(scalIt = this->mScalarVariables.begin(); scalIt != this->mScalarVariables.end(); scalIt++)
         {
            (*fIt)->addScalar((*scalIt));
         }

         // Loop over all vector variables
         std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>::iterator vectIt;
         for(vectIt = this->mVectorVariables.begin(); vectIt != this->mVectorVariables.end(); vectIt++)
         {
            (*fIt)->addVector((*vectIt));
         }

         if((*fIt)->space() == Dimensions::Space::PHYSICAL)
         {
            (*fIt)->setMesh(this->mTransformCoordinator.mesh());
         }
      }

      // init the output writers
      this->mSimIoCtrl.initWriters();
   }

}
