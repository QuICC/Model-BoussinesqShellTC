/** \file Simulation.cpp
 *  \brief Source of the high level simulation base 
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
#include <algorithm>

// External includes
//

// Class include
//
#include "Simulation/SimulationBase.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"
#include "IoTools/Formatter.hpp"
#include "SpectralOperators/BoundaryConditions.hpp"

namespace GeoMHDiSCC {

   SimulationBase::SimulationBase()
      : mExecutionTimer(true), mSimRunCtrl(), mDiagnostics()
   {
   }

   SimulationBase::~SimulationBase()
   {
   }

   void SimulationBase::initBase()
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

         // Initialise additional things depending on implementation
         this->initAdditionalBase();
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

   void SimulationBase::init(const SharedSimulationBoundary spBcs)
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

      // Sort the equations by type: time/solver/direct
      this->sortEquations();

      // Setup output files (ASCII diagnostics, state files, etc)
      this->setupOutput();

      // Initialise the diagnostics
      this->mDiagnostics.init(this->mTransformCoordinator.mesh(), this->mScalarVariables, this->mVectorVariables);

      // Cleanup IO control
      this->mSimIoCtrl.cleanup();

      // Print message to signal successful completion of initialisation step
      if(FrameworkMacro::allowsIO())
      {
         IoTools::Formatter::printLine(std::cout, '-');
         IoTools::Formatter::printCentered(std::cout, "Simulation initialisation successfull", '*');
         IoTools::Formatter::printLine(std::cout, '-');
         IoTools::Formatter::printNewline(std::cout);
      }

      // Debug statement
      DebuggerMacro_leave("init",0);
   }

   void SimulationBase::run()
   {
      // Final initialisation of the solvers
      this->initSolvers();

      // Debug statement
      DebuggerMacro_enter("run",0);

      // Stop timer and update initialisation time
      this->mExecutionTimer.stop();
      this->mExecutionTimer.update(ExecutionTimer::INIT);

      // Start timer
      this->mExecutionTimer.start();

      // Execute pre-run steps
      this->preRun();

      // Stop pre-run timing
      this->mExecutionTimer.stop();
      this->mExecutionTimer.update(ExecutionTimer::PRERUN);

      // Start timer
      this->mExecutionTimer.start();

      // Do main loop
      this->mainRun();

      // Stop main loop timing
      this->mExecutionTimer.stop();
      this->mExecutionTimer.update(ExecutionTimer::RUN);

      // Start timer
      this->mExecutionTimer.start();

      // Execute post-run operations
      this->postRun();

      // Synchronise computation nodes
      FrameworkMacro::synchronize();

      // Stop post-run timing
      this->mExecutionTimer.stop();
      this->mExecutionTimer.update(ExecutionTimer::POSTRUN);

      // Synchronise computation nodes
      FrameworkMacro::synchronize();

      // Debug statement
      DebuggerMacro_leave("run",0);
   }

   void SimulationBase::finalize()
   {
      // Debug statement
      DebuggerMacro_enter("finalize",0);

      // Print execution timer infos
      this->mExecutionTimer.printInfo(std::cout);

      // Print profiling infos (if required)
      ProfilerMacro_printInfo();

      // Print storage profiling infos (if required)
      StorageProfilerMacro_printInfo();

      // Debug statement
      DebuggerMacro_leave("finalize",0);
   }

   void SimulationBase::setInitialState(IoVariable::SharedIVariableHdf5Reader spInitFile)
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

   void SimulationBase::addOutputFile(int spOutFile)//SharedAscii spOutFile)
   {
      /// \mhdBug Fake implementation
   }

   void SimulationBase::addOutputFile(IoVariable::SharedIVariableHdf5NWriter spOutFile)
   {
      this->mSimIoCtrl.addOutputFile(spOutFile);
   }

   void SimulationBase::initSolvers()
   {
      // Debug statement
      DebuggerMacro_enter("initSolvers",1);

      // Init linear solver for diagnostic equations
      this->mLinearCoordinator.init(this->mScalarDiagnosticRange, this->mVectorDiagnosticRange);

      // Debug statement
      DebuggerMacro_leave("initSolvers",1);
   }

   void SimulationBase::computeNonlinear()
   {
      // Debug statement
      DebuggerMacro_enter("computeNonlinear",2);

      // Compute backward transform
      ProfilerMacro_start(ProfilerMacro::BWDTRANSFORM);
      this->mspBwdGrouper->transform(this->mScalarVariables, this->mVectorVariables, this->mTransformCoordinator);
      ProfilerMacro_stop(ProfilerMacro::BWDTRANSFORM);

      // compute nonlinear interaction and forward transform
      this->mspFwdGrouper->transform(this->mScalarEquations, this->mVectorEquations, this->mTransformCoordinator);

      // Debug statement
      DebuggerMacro_leave("computeNonlinear",2);
   }

   void SimulationBase::solveDiagnosticEquations()
   {
      // Debug statement
      DebuggerMacro_enter("solveDiagnosticEquations",2);

      DebuggerMacro_start("Solve diagnotic",3);
      ProfilerMacro_start(ProfilerMacro::DIAGNOSTICEQUATION);
      this->mLinearCoordinator.solve(this->mScalarDiagnosticRange, this->mVectorDiagnosticRange);
      ProfilerMacro_stop(ProfilerMacro::DIAGNOSTICEQUATION);
      DebuggerMacro_stop("Solve diagnostic t = ",3);

      // Debug statement
      DebuggerMacro_leave("solveDiagnosticEquations",2);
   }

   void SimulationBase::initVariables(VariableRequirement& varInfo)
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

   void SimulationBase::mapEquationVariables(std::set<PhysicalNames::Id>& nonInfo)
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
      
   void SimulationBase::initTransformCoordinator(const VariableRequirement& varInfo, const std::set<PhysicalNames::Id>& nonInfo)
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

   void SimulationBase::setupEquations(const SharedSimulationBoundary spBcs)
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

   void SimulationBase::setupOutput()
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

      // Allow for implementation specific tuning
      this->tuneOutput();

      // init the output writers
      this->mSimIoCtrl.initWriters();
   }

   void SimulationBase::tuneOutput()
   {
   }

   void SimulationBase::addConfigurationPart(IoConfig::SharedConfigurationReader spCfgFile)
   {
      // Empty to simplify implementations but can't be called from derived class
   }

   void SimulationBase::initAdditionalBase()
   {
      // Empty to simplify implementations but can't be called from derived class
   }

   void SimulationBase::sortEquations()
   {
      // Sort scalar equations
      std::stable_sort(this->mScalarEquations.begin(), this->mScalarEquations.end(), sortScalarEquationType);

      // Set the ranges for the different types
      int group = 1;
      std::vector<Equations::SharedIScalarEquation>::iterator  scalEqIt;
      for(scalEqIt = this->mScalarEquations.begin(); scalEqIt != this->mScalarEquations.end(); ++scalEqIt)
      {
         // Set time equation range
         if(group == 1 && computeScalarEquationType(*scalEqIt) == 2)
         {
            this->mScalarPrognosticRange = std::make_pair(this->mScalarEquations.begin(), scalEqIt);
            group++;
         }
         // Set solver and direct equation ranges
         if(group == 2 && computeScalarEquationType(*scalEqIt) == 3)
         {
            // Set solver equation range
            this->mScalarDiagnosticRange = std::make_pair(this->mScalarPrognosticRange.second, scalEqIt);

            // Set direct equation range
            this->mScalarDirectRange = std::make_pair(scalEqIt, this->mScalarEquations.end());

            // Exit loop
            break;
         }
      }

      // Sort vector equations
      std::stable_sort(this->mVectorEquations.begin(), this->mVectorEquations.end(), sortVectorEquationType);

      // Set the ranges for the different types
      group = 1;
      std::vector<Equations::SharedIVectorEquation>::iterator  vectEqIt;
      for(vectEqIt = this->mVectorEquations.begin(); vectEqIt != this->mVectorEquations.end(); ++vectEqIt)
      {
         // Set time equation range
         if(group == 1 && computeVectorEquationType(*vectEqIt) == 2)
         {
            this->mVectorPrognosticRange = std::make_pair(this->mVectorEquations.begin(), vectEqIt);
            group++;
         }
         // Set solver and direct equation ranges
         if(group == 2 && computeVectorEquationType(*vectEqIt) == 3)
         {
            // Set solver equation range
            this->mVectorDiagnosticRange = std::make_pair(this->mVectorPrognosticRange.second, vectEqIt);

            // Set direct equation range
            this->mVectorDirectRange = std::make_pair(vectEqIt, this->mVectorEquations.end());

            // Exit loop
            break;
         }
      }
   }

   int computeScalarEquationType(Equations::SharedIScalarEquation eqA)
   {
      return static_cast<int>(eqA->couplingInfo(FieldComponents::Spectral::SCALAR).equationType());
   }

   bool sortScalarEquationType(Equations::SharedIScalarEquation eqA, Equations::SharedIScalarEquation eqB)
   {
      return computeScalarEquationType(eqA) < computeScalarEquationType(eqB);
   }

   int computeVectorEquationType(Equations::SharedIVectorEquation eqA)
   { 
      return static_cast<int>(eqA->couplingInfo(FieldComponents::Spectral::ONE).equationType());
   }

   bool sortVectorEquationType(Equations::SharedIVectorEquation eqA, Equations::SharedIVectorEquation eqB)
   {
      return computeVectorEquationType(eqA) < computeVectorEquationType(eqB);
   }

}
