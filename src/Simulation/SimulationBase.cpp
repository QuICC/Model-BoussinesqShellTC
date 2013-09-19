/** 
 * @file Simulation.cpp
 * @brief Source of the high level simulation base 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Debug includes
//
#include "Debug/DebuggerMacro.h"

// Configuration includes
//
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
#include "Variables/RequirementTools.hpp"
#include "TransformCoordinators/TransformCoordinatorTools.hpp"

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
      RequirementTools::initVariables(varInfo, this->mScalarVariables, this->mVectorVariables, this->mScalarEquations, this->mVectorEquations, this->mspRes);

      // Storage for the nonlinear requirement info
      std::set<PhysicalNames::Id>   nonInfo;

      // Map variables to the equations and set nonlinear requirements
      RequirementTools::mapEquationVariables(nonInfo, this->mScalarEquations, this->mVectorEquations, this->mScalarVariables, this->mVectorVariables);

      // Initialise the transform coordinator
      this->initTransformCoordinator(varInfo, nonInfo);

      // Initialise the equations (generate operators, etc)
      this->setupEquations(spBcs);

      // Sort the equations by type: time/solver/trivial
      this->sortEquations();

      // Setup output files (ASCII diagnostics, state files, etc)
      this->setupOutput();

      // Get timestep information from configuration file
      Array tstep = this->mSimIoCtrl.configTimestepping();

      // Initialise the diagnostics
      this->mDiagnostics.init(this->mTransformCoordinator.mesh(), this->mScalarVariables, this->mVectorVariables, tstep);

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

   void SimulationBase::setInitialState(IoVariable::SharedStateFileReader spInitFile)
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

      // Addition operations on initial state file
      this->tuneInitialState(spInitFile);

      // Finalise file
      spInitFile->finalize();
   }

   void SimulationBase::tuneInitialState(IoVariable::SharedStateFileReader spInitFile)
   {
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
      DebuggerMacro_enter("solveTrivialEquations",3);

      DebuggerMacro_start("Solve trivial",4);
      ProfilerMacro_start(ProfilerMacro::TRIVIALEQUATION);
      this->mTrivialCoordinator.solve(this->mScalarTrivialRange, this->mVectorTrivialRange, this->mScalarVariables, this->mVectorVariables);
      ProfilerMacro_stop(ProfilerMacro::TRIVIALEQUATION);
      DebuggerMacro_stop("Solve trivial t = ",4);

      // Debug statement
      DebuggerMacro_leave("solveDiagnosticEquations",3);
   }

   void SimulationBase::solveDiagnosticEquations()
   {
      // Debug statement
      DebuggerMacro_enter("solveDiagnosticEquations",3);

      DebuggerMacro_start("Solve diagnotic",4);
      ProfilerMacro_start(ProfilerMacro::DIAGNOSTICEQUATION);
      this->mLinearCoordinator.solve(this->mScalarDiagnosticRange, this->mVectorDiagnosticRange, this->mScalarVariables, this->mVectorVariables);
      ProfilerMacro_stop(ProfilerMacro::DIAGNOSTICEQUATION);
      DebuggerMacro_stop("Solve diagnostic t = ",4);

      // Debug statement
      DebuggerMacro_leave("solveDiagnosticEquations",3);
   }
      
   void SimulationBase::initTransformCoordinator(const VariableRequirement& varInfo, const std::set<PhysicalNames::Id>& nonInfo)
   {
      // Extract the run options for the equation parameters
      std::map<NonDimensional::Id,MHDFloat> runOptions;
      std::vector<NonDimensional::Id>::iterator it;
      for(it = this->mspEqParams->ids().begin(); it != this->mspEqParams->ids().end(); ++it)
      {
         runOptions.insert(std::make_pair(*it, this->mspEqParams->nd(*it)));
      }

      // Initialise the transform coordinator
      Transform::TransformCoordinatorTools::init(this->mTransformCoordinator, this->mspFwdGrouper, this->mspBwdGrouper, varInfo, nonInfo, this->mspRes, runOptions);
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

      // Initialise empty scalar ranges
      this->mScalarPrognosticRange = std::make_pair(this->mScalarEquations.end(), this->mScalarEquations.end());
      this->mScalarDiagnosticRange = std::make_pair(this->mScalarEquations.end(), this->mScalarEquations.end());
      this->mScalarTrivialRange = std::make_pair(this->mScalarEquations.end(), this->mScalarEquations.end());

      // Determine the ranges for the different types
      std::vector<Equations::SharedIScalarEquation>::iterator  scalEqIt;

      scalEqIt = std::find_if(this->mScalarEquations.begin(), this->mScalarEquations.end(), scalarIsPrognostic);
      if(scalEqIt != this->mScalarEquations.end())
      {
         this->mScalarPrognosticRange = std::equal_range(this->mScalarEquations.begin(), this->mScalarEquations.end(), *scalEqIt, sortScalarEquationType);
      }

      scalEqIt = std::find_if(this->mScalarEquations.begin(), this->mScalarEquations.end(), scalarIsDiagnostic);
      if(scalEqIt != this->mScalarEquations.end())
      {
         this->mScalarDiagnosticRange = std::equal_range(this->mScalarEquations.begin(), this->mScalarEquations.end(), *scalEqIt, sortScalarEquationType);
      }

      scalEqIt = std::find_if(this->mScalarEquations.begin(), this->mScalarEquations.end(), scalarIsTrivial);
      if(scalEqIt != this->mScalarEquations.end())
      {
         this->mScalarTrivialRange = std::equal_range(this->mScalarEquations.begin(), this->mScalarEquations.end(), *scalEqIt, sortScalarEquationType);
      }

      // Sort vector equations
      std::stable_sort(this->mVectorEquations.begin(), this->mVectorEquations.end(), sortVectorEquationType);

      // Initialise empty vector ranges
      this->mVectorPrognosticRange = std::make_pair(this->mVectorEquations.end(), this->mVectorEquations.end());
      this->mVectorDiagnosticRange = std::make_pair(this->mVectorEquations.end(), this->mVectorEquations.end());
      this->mVectorTrivialRange = std::make_pair(this->mVectorEquations.end(), this->mVectorEquations.end());

      // Determine the ranges for the different types
      std::vector<Equations::SharedIVectorEquation>::iterator  vectEqIt;

      vectEqIt = std::find_if(this->mVectorEquations.begin(), this->mVectorEquations.end(), vectorIsPrognostic);
      if(vectEqIt != this->mVectorEquations.end())
      {
         this->mVectorPrognosticRange = std::equal_range(this->mVectorEquations.begin(), this->mVectorEquations.end(), *vectEqIt, sortVectorEquationType);
      }

      vectEqIt = std::find_if(this->mVectorEquations.begin(), this->mVectorEquations.end(), vectorIsDiagnostic);
      if(vectEqIt != this->mVectorEquations.end())
      {
         this->mVectorDiagnosticRange = std::equal_range(this->mVectorEquations.begin(), this->mVectorEquations.end(), *vectEqIt, sortVectorEquationType);
      }

      vectEqIt = std::find_if(this->mVectorEquations.begin(), this->mVectorEquations.end(), vectorIsTrivial);
      if(vectEqIt != this->mVectorEquations.end())
      {
         this->mVectorTrivialRange = std::equal_range(this->mVectorEquations.begin(), this->mVectorEquations.end(), *vectEqIt, sortVectorEquationType);
      }

      // Set the ranges for the different types

      // Identifiy the solver indexes by analysing the coupling between the equations
      DebuggerMacro_enter("identifyCoupling_Prognostic",1);
      DebuggerMacro_showValue("---> Prognostic scalar equations: ",1, this->mScalarPrognosticRange.second - this->mScalarPrognosticRange.first);
      DebuggerMacro_showValue("---> Prognostic vector equations: ",1, this->mVectorPrognosticRange.second - this->mVectorPrognosticRange.first);
      this->identifyCoupling(this->mScalarPrognosticRange, this->mVectorPrognosticRange);
      DebuggerMacro_leave("identifyCoupling_Prognostic",1);

      DebuggerMacro_enter("identifyCoupling_Diagnostic",1);
      DebuggerMacro_showValue("---> Diagnostic scalar equations: ",1, this->mScalarDiagnosticRange.second - this->mScalarDiagnosticRange.first);
      DebuggerMacro_showValue("---> Diagnostic vector equations: ",1, this->mVectorDiagnosticRange.second - this->mVectorDiagnosticRange.first);
      this->identifyCoupling(this->mScalarDiagnosticRange, this->mVectorDiagnosticRange);
      DebuggerMacro_leave("identifyCoupling_Diagnostic",1);

      DebuggerMacro_enter("identifyCoupling_Trivial",1);
      DebuggerMacro_showValue("---> Trivial scalar equations: ",1, this->mScalarTrivialRange.second - this->mScalarTrivialRange.first);
      DebuggerMacro_showValue("---> Trivial vector equations: ",1, this->mVectorTrivialRange.second - this->mVectorTrivialRange.first);
      this->identifyCoupling(this->mScalarTrivialRange, this->mVectorTrivialRange);
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
      for(vectEqIt = vectEq.first; vectEqIt != vectEq.second; ++vectEqIt)
      {
         // Get coupled counter for each component
         ArrayI counter((*vectEqIt)->nSpectral());
         counter.setConstant(0);

         // Loop over the (identified) scalar equations
         for(doneSEqIt = scalEq.first; doneSEqIt != scalEq.second; ++doneSEqIt)
         {
            // loop over the implicit range
            Equations::CouplingInformation::FieldId_iterator fIt;
            Equations::CouplingInformation::FieldId_range fRange = (*doneSEqIt)->couplingInfo(FieldComponents::Spectral::SCALAR).implicitRange();
            for(fIt = fRange.first; fIt != fRange.second; ++fIt)
            {
               Equations::IVectorEquation::SpectralComponent_iterator compIt;
               Equations::IVectorEquation::SpectralComponent_range  compRange = (*vectEqIt)->spectralRange();
               int i = 0;
               for(compIt = compRange.first; compIt != compRange.second; ++compIt, ++i)
               {
                  // Check if field's first component is in implicit range
                  if(counter(i) == 0 && *fIt == std::make_pair((*vectEqIt)->name(), *compIt))
                  {
                     // Set the solver index
                     (*vectEqIt)->setSolverIndex(*compIt, (*doneSEqIt)->couplingInfo(FieldComponents::Spectral::SCALAR).solverIndex());

                     // Debug statements
                     DebuggerMacro_showValue("Identified coupled vector solver: ", 2, (*vectEqIt)->name());
                     DebuggerMacro_showValue("---> component: ", 2, *compIt);
                     DebuggerMacro_showValue("---> solver index: ", 2, (*vectEqIt)->couplingInfo(*compIt).solverIndex());
                     DebuggerMacro_showValue("---> is complex? ", 2, (*vectEqIt)->couplingInfo(*compIt).isComplex());

                     // Set coupling flag and break out
                     counter(i) = 1;
                  }
               }

               // Break out of loop if already identified
               if(counter.sum() == counter.size())
               {
                  break;
               }
            }

            // Break out of loop if already identified
            if(counter.sum() == counter.size())
            {
               break;
            }
         }

         if(counter.sum() < counter.size())
         {
            // Loop over the (identified) vector equations
            for(doneVEqIt = vectEq.first; doneVEqIt != vectEqIt; ++doneVEqIt)
            {
               Equations::IVectorEquation::SpectralComponent_iterator doneIt;
               Equations::IVectorEquation::SpectralComponent_range  doneRange = (*doneVEqIt)->spectralRange();
               for(doneIt = doneRange.first; doneIt != doneRange.second; ++doneIt)
               {
                  // loop over the implicit range of identified components
                  Equations::CouplingInformation::FieldId_iterator fIt;
                  Equations::CouplingInformation::FieldId_range fRange = (*doneVEqIt)->couplingInfo(*doneIt).implicitRange();
                  for(fIt = fRange.first; fIt != fRange.second; ++fIt)
                  {
                     Equations::IVectorEquation::SpectralComponent_iterator compIt;
                     Equations::IVectorEquation::SpectralComponent_range  compRange = (*vectEqIt)->spectralRange();
                     int i = 0;
                     for(compIt = compRange.first; compIt != compRange.second; ++compIt, ++i)
                     {
                        // Check if field's first component is in implicit range
                        if(counter(i) == 0 && *fIt == std::make_pair((*vectEqIt)->name(), *compIt))
                        {
                           // Set the solver index
                           (*vectEqIt)->setSolverIndex(*compIt, (*doneVEqIt)->couplingInfo(*doneIt).solverIndex());

                           // Debug statements
                           DebuggerMacro_showValue("Identified coupled vector solver: ", 2, (*vectEqIt)->name());
                           DebuggerMacro_showValue("---> component: ", 2, *compIt);
                           DebuggerMacro_showValue("---> solver index: ", 2, (*vectEqIt)->couplingInfo(*compIt).solverIndex());
                           DebuggerMacro_showValue("---> is complex? ", 2, (*vectEqIt)->couplingInfo(*compIt).isComplex());

                           // Set coupling flag and break out
                           counter(i) = 1;
                        }
                     }


                     // Break out of loop if already identified
                     if(counter.sum() == counter.size())
                     {
                        break;
                     }
                  }

                  // Break out of loop if already identified
                  if(counter.sum() == counter.size())
                  {
                     break;
                  }
               }
            }
         }

         // All checked and equation is not coupled
         if(counter.sum() < counter.size())
         {
            Equations::IVectorEquation::SpectralComponent_iterator compIt;
            Equations::IVectorEquation::SpectralComponent_range  compRange = (*vectEqIt)->spectralRange();
            int i = 0;
            for(compIt = compRange.first; compIt != compRange.second; ++compIt, ++i)
            {
               if(counter(i) == 0)
               {
                  // Set solver index for complex equation
                  if((*vectEqIt)->couplingInfo(*compIt).isComplex())
                  {
                     // Set complex solver index
                     (*vectEqIt)->setSolverIndex(*compIt, zIdx);

                     // Increment complex solver index
                     zIdx++;

                     // Set solver index for real equation
                  } else
                  {
                     // Set real solver index
                     (*vectEqIt)->setSolverIndex(*compIt, dIdx);

                     // Increment real solver index
                     dIdx++;
                  }

                  // Debug statements
                  DebuggerMacro_showValue("Identified first vector solver: ", 2, (*vectEqIt)->name());
                  DebuggerMacro_showValue("---> component: ", 2, *compIt);
                  DebuggerMacro_showValue("---> solver index: ", 2, (*vectEqIt)->couplingInfo(*compIt).solverIndex());
                  DebuggerMacro_showValue("---> is complex? ", 2, (*vectEqIt)->couplingInfo(*compIt).isComplex());
               }
            }
         }

         // Reset coupling flag
         counter.setConstant(0);
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

   bool scalarIsPrognostic(Equations::SharedIScalarEquation eqA)
   {
      return eqA->couplingInfo(FieldComponents::Spectral::SCALAR).equationType() == Equations::CouplingInformation::PROGNOSTIC;
   }

   bool scalarIsDiagnostic(Equations::SharedIScalarEquation eqA)
   {
      return eqA->couplingInfo(FieldComponents::Spectral::SCALAR).equationType() == Equations::CouplingInformation::DIAGNOSTIC;
   }

   bool scalarIsTrivial(Equations::SharedIScalarEquation eqA)
   {
      return eqA->couplingInfo(FieldComponents::Spectral::SCALAR).equationType() == Equations::CouplingInformation::TRIVIAL;
   }

   bool scalarIsWrapper(Equations::SharedIScalarEquation eqA)
   {
      return eqA->couplingInfo(FieldComponents::Spectral::SCALAR).equationType() == Equations::CouplingInformation::WRAPPER;
   }

   bool vectorIsPrognostic(Equations::SharedIVectorEquation eqA)
   {
      return eqA->couplingInfo(FieldComponents::Spectral::ONE).equationType() == Equations::CouplingInformation::PROGNOSTIC;
   }

   bool vectorIsDiagnostic(Equations::SharedIVectorEquation eqA)
   {
      return eqA->couplingInfo(FieldComponents::Spectral::ONE).equationType() == Equations::CouplingInformation::DIAGNOSTIC;
   }

   bool vectorIsTrivial(Equations::SharedIVectorEquation eqA)
   {
      return eqA->couplingInfo(FieldComponents::Spectral::ONE).equationType() == Equations::CouplingInformation::TRIVIAL;
   }

   bool vectorIsWrapper(Equations::SharedIVectorEquation eqA)
   {
      return eqA->couplingInfo(FieldComponents::Spectral::ONE).equationType() == Equations::CouplingInformation::WRAPPER;
   }
}
