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

#include <iostream>
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

      // Sort the equations by type: time/solver/trivial
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

      // Init trivial solver for trivial equations
      this->mTrivialCoordinator.init(this->mScalarTrivialRange, this->mVectorTrivialRange);

      // Init linear solver for trivial equations
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

   void SimulationBase::solveTrivialEquations()
   {
      // Debug statement
      DebuggerMacro_enter("solveTrivialEquations",2);

      DebuggerMacro_start("Solve trivial",3);
      ProfilerMacro_start(ProfilerMacro::TRIVIALEQUATION);
      this->mTrivialCoordinator.solve(this->mScalarTrivialRange, this->mVectorTrivialRange);
      ProfilerMacro_stop(ProfilerMacro::TRIVIALEQUATION);
      DebuggerMacro_stop("Solve trivial t = ",3);

      // Debug statement
      DebuggerMacro_leave("solveDiagnosticEquations",2);
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

      // Initialise scalar ranges: prognostic is full, the two other are empty
      this->mScalarPrognosticRange = std::make_pair(this->mScalarEquations.begin(), this->mScalarEquations.end());
      this->mScalarDiagnosticRange = std::make_pair(this->mScalarEquations.end(), this->mScalarEquations.end());
      this->mScalarTrivialRange = std::make_pair(this->mScalarEquations.end(), this->mScalarEquations.end());

      // Determine the ranges for the different types
      int group = 1;
      std::vector<Equations::SharedIScalarEquation>::iterator  scalEqIt;
      for(scalEqIt = this->mScalarEquations.begin(); scalEqIt != this->mScalarEquations.end(); ++scalEqIt)
      {
         // Transition from prognostic to the two others
         if(group == 1 && computeScalarEquationType(*scalEqIt) > 1)
         {
            // No prognostic equation present
            if(scalEqIt == this->mScalarEquations.begin())
            {
               this->mScalarPrognosticRange = std::make_pair(this->mScalarEquations.end(), this->mScalarEquations.end());
            
            // With prognostic equations present
            } else
            {
               this->mScalarPrognosticRange = std::make_pair(this->mScalarEquations.begin(), scalEqIt);
            }

            // Transitions to diagnostic equations
            if(computeScalarEquationType(*scalEqIt) == 2)
            {
               // Set diagnostic equation range
               this->mScalarDiagnosticRange = std::make_pair(scalEqIt, this->mScalarEquations.end());

            // Transitions to trivial equations
            } else
            {
               // Set trivial equation range
               this->mScalarTrivialRange = std::make_pair(scalEqIt, this->mScalarEquations.end());
               break;
            }
            group++;
         }

         // Transition from diagnostic to trivial
         if(group == 2 && computeScalarEquationType(*scalEqIt) == 3)
         {
            // Set diagnostic equation range
            this->mScalarDiagnosticRange.second = scalEqIt;

            // Set trivial equation range
            this->mScalarTrivialRange = std::make_pair(scalEqIt, this->mScalarEquations.end());

            // Exit loop
            break;
         }
      }

      // Sort vector equations
      std::stable_sort(this->mVectorEquations.begin(), this->mVectorEquations.end(), sortVectorEquationType);

      // Initialise vector ranges: prognostic is full, the two other are empty
      this->mVectorPrognosticRange = std::make_pair(this->mVectorEquations.begin(), this->mVectorEquations.end());
      this->mVectorDiagnosticRange = std::make_pair(this->mVectorEquations.end(), this->mVectorEquations.end());
      this->mVectorTrivialRange = std::make_pair(this->mVectorEquations.end(), this->mVectorEquations.end());

      // Set the ranges for the different types
      group = 1;
      std::vector<Equations::SharedIVectorEquation>::iterator  vectEqIt;
      for(vectEqIt = this->mVectorEquations.begin(); vectEqIt != this->mVectorEquations.end(); ++vectEqIt)
      {
         // Transition from prognostic to the two others
         if(group == 1 && computeVectorEquationType(*vectEqIt) > 1)
         {
            // No prognostic equation present
            if(vectEqIt == this->mVectorEquations.begin())
            {
               this->mVectorPrognosticRange = std::make_pair(this->mVectorEquations.end(), this->mVectorEquations.end());
            
            // With prognostic equations present
            } else
            {
               this->mVectorPrognosticRange = std::make_pair(this->mVectorEquations.begin(), vectEqIt);
            }

            // Transitions to diagnostic equations
            if(computeVectorEquationType(*vectEqIt) == 2)
            {
               // Set diagnostic equation range
               this->mVectorDiagnosticRange = std::make_pair(vectEqIt, this->mVectorEquations.end());

            // Transitions to trivial equations
            } else
            {
               // Set trivial equation range
               this->mVectorTrivialRange = std::make_pair(vectEqIt, this->mVectorEquations.end());
               break;
            }
            group++;
         }

         // Transition from diagnostic to trivial
         if(group == 2 && computeVectorEquationType(*vectEqIt) == 3)
         {
            // Set diagnostic equation range
            this->mVectorDiagnosticRange.second = vectEqIt;

            // Set trivial equation range
            this->mVectorTrivialRange = std::make_pair(vectEqIt, this->mVectorEquations.end());

            // Exit loop
            break;
         }
      }

      // Identifiy the solver indexes by analysing the coupling between the equations
      DebuggerMacro_enter("identifyCoupling_Prognostic",1);
      this->identifyCoupling(this->mScalarPrognosticRange, this->mVectorPrognosticRange);
      DebuggerMacro_leave("identifyCoupling_Prognostic",1);

      DebuggerMacro_enter("identifyCoupling_Diagnostic",1);
      this->identifyCoupling(this->mScalarDiagnosticRange, this->mVectorPrognosticRange);
      DebuggerMacro_leave("identifyCoupling_Diagnostic",1);

      DebuggerMacro_enter("identifyCoupling_Trivial",1);
      this->identifyCoupling(this->mScalarTrivialRange, this->mVectorPrognosticRange);
      DebuggerMacro_leave("identifyCoupling_Trivial",1);
   }

   void SimulationBase::identifyCoupling(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
   {
      // Iterators for scalar equations
      std::vector<Equations::SharedIScalarEquation>::iterator  scalEqIt;
      std::vector<Equations::SharedIScalarEquation>::iterator  doneSEqIt;

      // Current solver indexes for real and complex solvers
      int dIdx = 0;
      int zIdx = 0;

      // Coupling flag
      bool coupled = false;

      // Field identification
      FieldComponents::Spectral::Id compId = FieldComponents::Spectral::SCALAR;
      SpectralFieldId   fieldId;
      
      // Loop over the scalar equations
      for(scalEqIt = scalEq.first; scalEqIt != scalEq.second; ++scalEqIt)
      {
         // Build field identity
         fieldId = std::make_pair((*scalEqIt)->name(), compId);

         // Loop over already identified equations
         for(doneSEqIt = scalEq.first; doneSEqIt != scalEqIt; ++doneSEqIt)
         {
            // loop over the implicit range
            Equations::CouplingInformation::FieldId_iterator fIt;
            Equations::CouplingInformation::FieldId_range fRange = (*doneSEqIt)->couplingInfo(compId).implicitRange();
            for(fIt = fRange.first; fIt != fRange.second; ++fIt)
            {
               // Check if field is in implicit range
               if(*fIt == fieldId)
               {
                  // Set the solver index
                  (*scalEqIt)->setSolverIndex(compId, (*doneSEqIt)->couplingInfo(compId).solverIndex());
                  
                  // Debug statements
                  DebuggerMacro_showValue("Identified coupled scalar solver: ", 2, (*scalEqIt)->name());
                  DebuggerMacro_showValue("---> solver index: ", 2, (*scalEqIt)->couplingInfo(compId).solverIndex());
                  DebuggerMacro_showValue("---> is complex? ", 2, (*scalEqIt)->couplingInfo(compId).isComplex());
                  
                  // Set coupling flag and break out
                  coupled = true;
                  break;
               }
            }

            // Break out of loop if already identified
            if(coupled)
            {
               break;
            }
         }

         // All checked and equation is not coupled
         if(!coupled)
         {
            // Set solver index for complex equation
            if((*scalEqIt)->couplingInfo(compId).isComplex())
            {
               // Set complex solver index
               (*scalEqIt)->setSolverIndex(compId, zIdx);

               // Increment complex solver index
               zIdx++;

            // Set solver index for real equation
            } else
            {
               // Set real solver index
               (*scalEqIt)->setSolverIndex(compId, dIdx);
               
               // Increment real solver index
               dIdx++;
            }

            // Debug statements
            DebuggerMacro_showValue("Identified first scalar solver: ", 2, (*scalEqIt)->name());
            DebuggerMacro_showValue("---> solver index: ", 2, (*scalEqIt)->couplingInfo(compId).solverIndex());
            DebuggerMacro_showValue("---> is complex? ", 2, (*scalEqIt)->couplingInfo(compId).isComplex());
         }

         // Reset coupling flag
         coupled = false;
      }

      // Iterators for vector equations
      std::vector<Equations::SharedIVectorEquation>::iterator  vectEqIt;
      std::vector<Equations::SharedIVectorEquation>::iterator  doneVEqIt;
      
      // Loop over the vector equations
      bool coupledOne = false;
      bool coupledTwo = false;
      for(vectEqIt = vectEq.first; vectEqIt != vectEq.second; ++vectEqIt)
      {
         // Loop over the (identified) scalar equations
         for(doneSEqIt = scalEq.first; doneSEqIt != scalEq.second; ++doneSEqIt)
         {
            // loop over the implicit range
            Equations::CouplingInformation::FieldId_iterator fIt;
            Equations::CouplingInformation::FieldId_range fRange = (*doneSEqIt)->couplingInfo(FieldComponents::Spectral::SCALAR).implicitRange();
            for(fIt = fRange.first; fIt != fRange.second; ++fIt)
            {
               // Check if field's first component is in implicit range
               compId = FieldComponents::Spectral::ONE;
               if(!coupledOne && *fIt == std::make_pair((*vectEqIt)->name(), compId))
               {
                  // Set the solver index
                  (*vectEqIt)->setSolverIndex(compId, (*doneSEqIt)->couplingInfo(FieldComponents::Spectral::SCALAR).solverIndex());
                  
                  // Debug statements
                  DebuggerMacro_showValue("Identified coupled vector solver: ", 2, (*vectEqIt)->name());
                  DebuggerMacro_showValue("---> component: ", 2, compId);
                  DebuggerMacro_showValue("---> solver index: ", 2, (*vectEqIt)->couplingInfo(compId).solverIndex());
                  DebuggerMacro_showValue("---> is complex? ", 2, (*vectEqIt)->couplingInfo(compId).isComplex());
                  
                  // Set coupling flag and break out
                  coupledOne = true;
               }

               // Check if field's second component is in implicit range
               compId = FieldComponents::Spectral::TWO;
               if(!coupledTwo && *fIt == std::make_pair((*vectEqIt)->name(), compId))
               {
                  // Set the solver index
                  (*vectEqIt)->setSolverIndex(compId, (*doneSEqIt)->couplingInfo(FieldComponents::Spectral::SCALAR).solverIndex());
                  
                  // Debug statements
                  DebuggerMacro_showValue("Identified coupled vector solver: ", 2, (*vectEqIt)->name());
                  DebuggerMacro_showValue("---> component: ", 2, compId);
                  DebuggerMacro_showValue("---> solver index: ", 2, (*vectEqIt)->couplingInfo(compId).solverIndex());
                  DebuggerMacro_showValue("---> is complex? ", 2, (*vectEqIt)->couplingInfo(compId).isComplex());
                  
                  // Set coupling flag and break out
                  coupledTwo = true;
               }

               // Break out of loop if already identified
               if(coupledOne && coupledTwo)
               {
                  break;
               }
            }

            // Break out of loop if already identified
            if(coupledOne && coupledTwo)
            {
               break;
            }
         }

         if(!coupledOne || !coupledTwo)
         {
            // Loop over the (identified) vector equations
            for(doneVEqIt = vectEq.first; doneVEqIt != vectEqIt; ++doneVEqIt)
            {
               // loop over the implicit range of first component
               FieldComponents::Spectral::Id doneId = FieldComponents::Spectral::ONE;
               Equations::CouplingInformation::FieldId_iterator fIt;
               Equations::CouplingInformation::FieldId_range fRange = (*doneVEqIt)->couplingInfo(doneId).implicitRange();
               for(fIt = fRange.first; fIt != fRange.second; ++fIt)
               {
                  // Check if field's first component is in implicit range
                  compId = FieldComponents::Spectral::ONE;
                  if(!coupledOne && *fIt == std::make_pair((*vectEqIt)->name(), compId))
                  {
                     // Set the solver index
                     (*vectEqIt)->setSolverIndex(compId, (*doneVEqIt)->couplingInfo(doneId).solverIndex());

                     // Debug statements
                     DebuggerMacro_showValue("Identified coupled vector solver: ", 2, (*vectEqIt)->name());
                     DebuggerMacro_showValue("---> component: ", 2, compId);
                     DebuggerMacro_showValue("---> solver index: ", 2, (*vectEqIt)->couplingInfo(compId).solverIndex());
                     DebuggerMacro_showValue("---> is complex? ", 2, (*vectEqIt)->couplingInfo(compId).isComplex());

                     // Set coupling flag and break out
                     coupledOne = true;
                  }

                  // Check if field's first component is in implicit range
                  compId = FieldComponents::Spectral::TWO;
                  if(!coupledTwo && *fIt == std::make_pair((*vectEqIt)->name(), compId))
                  {
                     // Set the solver index
                     (*vectEqIt)->setSolverIndex(compId, (*doneVEqIt)->couplingInfo(doneId).solverIndex());

                     // Debug statements
                     DebuggerMacro_showValue("Identified coupled vector solver: ", 2, (*vectEqIt)->name());
                     DebuggerMacro_showValue("---> component: ", 2, compId);
                     DebuggerMacro_showValue("---> solver index: ", 2, (*vectEqIt)->couplingInfo(compId).solverIndex());
                     DebuggerMacro_showValue("---> is complex? ", 2, (*vectEqIt)->couplingInfo(compId).isComplex());

                     // Set coupling flag and break out
                     coupledTwo = true;
                  }

                  // Break out of loop if already identified
                  if(coupledOne && coupledTwo)
                  {
                     break;
                  }
               }

               // loop over the implicit range of second component
               doneId = FieldComponents::Spectral::TWO;
               fRange = (*doneVEqIt)->couplingInfo(doneId).implicitRange();
               for(fIt = fRange.first; fIt != fRange.second; ++fIt)
               {
                  // Check if field's first component is in implicit range
                  compId = FieldComponents::Spectral::ONE;
                  if(!coupledOne && *fIt == std::make_pair((*vectEqIt)->name(), compId))
                  {
                     // Set the solver index
                     (*vectEqIt)->setSolverIndex(compId, (*doneVEqIt)->couplingInfo(doneId).solverIndex());

                     // Debug statements
                     DebuggerMacro_showValue("Identified coupled vector solver: ", 2, (*vectEqIt)->name());
                     DebuggerMacro_showValue("---> component: ", 2, compId);
                     DebuggerMacro_showValue("---> solver index: ", 2, (*vectEqIt)->couplingInfo(compId).solverIndex());
                     DebuggerMacro_showValue("---> is complex? ", 2, (*vectEqIt)->couplingInfo(compId).isComplex());

                     // Set coupling flag and break out
                     coupledOne = true;
                  }

                  // Check if field's first component is in implicit range
                  compId = FieldComponents::Spectral::TWO;
                  if(!coupledTwo && *fIt == std::make_pair((*vectEqIt)->name(), compId))
                  {
                     // Set the solver index
                     (*vectEqIt)->setSolverIndex(compId, (*doneVEqIt)->couplingInfo(doneId).solverIndex());

                     // Debug statements
                     DebuggerMacro_showValue("Identified coupled vector solver: ", 2, (*vectEqIt)->name());
                     DebuggerMacro_showValue("---> component: ", 2, compId);
                     DebuggerMacro_showValue("---> solver index: ", 2, (*vectEqIt)->couplingInfo(compId).solverIndex());
                     DebuggerMacro_showValue("---> is complex? ", 2, (*vectEqIt)->couplingInfo(compId).isComplex());

                     // Set coupling flag and break out
                     coupledTwo = true;
                  }

                  // Break out of loop if already identified
                  if(coupledOne && coupledTwo)
                  {
                     break;
                  }
               }

               // Break out of loop if already identified
               if(coupledOne && coupledTwo)
               {
                  break;
               }
            }
         }

         // All checked and equation is not coupled
         if(!coupledOne || !coupledTwo)
         {
            if(!coupledOne)
            {
               compId = FieldComponents::Spectral::ONE;

               // Set solver index for complex equation
               if((*vectEqIt)->couplingInfo(compId).isComplex())
               {
                  // Set complex solver index
                  (*vectEqIt)->setSolverIndex(compId, zIdx);

                  // Increment complex solver index
                  zIdx++;

               // Set solver index for real equation
               } else
               {
                  // Set real solver index
                  (*vectEqIt)->setSolverIndex(compId, dIdx);

                  // Increment real solver index
                  dIdx++;
               }

               // Debug statements
               DebuggerMacro_showValue("Identified first vector solver: ", 2, (*vectEqIt)->name());
               DebuggerMacro_showValue("---> component: ", 2, compId);
               DebuggerMacro_showValue("---> solver index: ", 2, (*vectEqIt)->couplingInfo(compId).solverIndex());
               DebuggerMacro_showValue("---> is complex? ", 2, (*vectEqIt)->couplingInfo(compId).isComplex());
            }

            if(!coupledTwo)
            {
               compId = FieldComponents::Spectral::TWO;

               // Set solver index for complex equation
               if((*vectEqIt)->couplingInfo(compId).isComplex())
               {
                  // Set complex solver index
                  (*vectEqIt)->setSolverIndex(compId, zIdx);

                  // Increment complex solver index
                  zIdx++;

               // Set solver index for real equation
               } else
               {
                  // Set real solver index
                  (*vectEqIt)->setSolverIndex(compId, dIdx);

                  // Increment real solver index
                  dIdx++;
               }

               // Debug statements
               DebuggerMacro_showValue("Identified first vector solver: ", 2, (*vectEqIt)->name());
               DebuggerMacro_showValue("---> component: ", 2, compId);
               DebuggerMacro_showValue("---> solver index: ", 2, (*vectEqIt)->couplingInfo(compId).solverIndex());
               DebuggerMacro_showValue("---> is complex? ", 2, (*vectEqIt)->couplingInfo(compId).isComplex());
            }
         }

         // Reset coupling flag
         coupledOne = false;
         coupledTwo = false;
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
