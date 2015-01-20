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
#include "Variables/RequirementTools.hpp"
#include "TransformCoordinators/TransformCoordinatorTools.hpp"
#include "Equations/Tools/EquationTools.hpp"

namespace GeoMHDiSCC {

   SimulationBase::SimulationBase()
      : mExecutionTimer(true), mSimRunCtrl(), mDiagnostics(), mForwardIsNonlinear(true)
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

      // Transform projector tree
      std::vector<Transform::ProjectorTree> projectorTree;

      // Initialise the variables and set general variable requirements
      RequirementTools::initVariables(projectorTree, this->mScalarVariables, this->mVectorVariables, this->mScalarEquations, this->mVectorEquations, this->mspRes);

      // Transform integrator tree
      std::vector<Transform::IntegratorTree> integratorTree;

      // Map variables to the equations and set nonlinear requirements
      RequirementTools::mapEquationVariables(integratorTree, this->mScalarEquations, this->mVectorEquations, this->mScalarVariables, this->mVectorVariables, this->mForwardIsNonlinear);

      // Initialise the transform coordinator
      this->initTransformCoordinator(integratorTree, projectorTree);

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

      // Print simulation run infos
      this->mSimRunCtrl.printInfo(std::cout);

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

      // Forward state file time and timestep to diagnostic coordinator
      this->mDiagnostics.useStateTime(spInitFile->time(), spInitFile->timestep());

      // Finalise file
      spInitFile->finalize();
   }

   void SimulationBase::tuneInitialState(IoVariable::SharedStateFileReader spInitFile)
   {
   }

   void SimulationBase::addAsciiOutputFile(IoVariable::SharedIVariableAsciiEWriter spOutFile)
   {
      this->mSimIoCtrl.addAsciiOutputFile(spOutFile);
   }

   void SimulationBase::addHdf5OutputFile(IoVariable::SharedIVariableHdf5NWriter spOutFile)
   {
      this->mSimIoCtrl.addHdf5OutputFile(spOutFile);
   }

   void SimulationBase::initSolvers()
   {
      // Debug statement
      DebuggerMacro_enter("initSolvers",1);

      // Print message to signal start of timestepper building
      if(FrameworkMacro::allowsIO())
      {
         IoTools::Formatter::printCentered(std::cout, "(... Building Solvers ...)", ' ');
         IoTools::Formatter::printNewline(std::cout);
      }

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
      DebuggerMacro_start("computeNonlinear",3);

      // Compute backward transform
      ProfilerMacro_start(ProfilerMacro::BWDTRANSFORM);
      this->mspBwdGrouper->transform(this->mScalarVariables, this->mVectorVariables, this->mTransformCoordinator);
      ProfilerMacro_stop(ProfilerMacro::BWDTRANSFORM);

      // compute nonlinear interaction and forward transform
      this->mspFwdGrouper->transform(this->mScalarEquations, this->mVectorEquations, this->mTransformCoordinator);

      // Debug statement
      DebuggerMacro_stop("computeNonlinear t = ",3);
      DebuggerMacro_leave("computeNonlinear",2);
   }

   void SimulationBase::explicitTrivialEquations(const SolveTiming::Id time, const ExplicitTiming::Id expTime)
   {
      // Debug statement
      DebuggerMacro_enter("explicitTrivialEquations",3);

      DebuggerMacro_start("Explicit trivial",4);
      ProfilerMacro_start(ProfilerMacro::TRIVIALEQUATION);
      this->mTrivialCoordinator.setSolveTime(time);
      this->mTrivialCoordinator.setExplicitTime(expTime);
      this->mTrivialCoordinator.getExplicitInput(this->mScalarTrivialRange, this->mVectorTrivialRange, this->mScalarVariables, this->mVectorVariables);
      ProfilerMacro_stop(ProfilerMacro::TRIVIALEQUATION);
      DebuggerMacro_stop("Explicit trivial t = ",4);

      // Debug statement
      DebuggerMacro_leave("explicitTrivialEquations",3);
   }

   void SimulationBase::explicitDiagnosticEquations(const SolveTiming::Id time, const ExplicitTiming::Id expTime)
   {
      // Debug statement
      DebuggerMacro_enter("explicitDiagnosticEquations",3);

      DebuggerMacro_start("Explicit diagnostic",4);
      ProfilerMacro_start(ProfilerMacro::DIAGNOSTICEQUATION);
      this->mTrivialCoordinator.setSolveTime(time);
      this->mLinearCoordinator.setExplicitTime(expTime);
      this->mLinearCoordinator.getExplicitInput(this->mScalarDiagnosticRange, this->mVectorDiagnosticRange, this->mScalarVariables, this->mVectorVariables);
      ProfilerMacro_stop(ProfilerMacro::DIAGNOSTICEQUATION);
      DebuggerMacro_stop("Explicit diagnostic t = ",4);

      // Debug statement
      DebuggerMacro_leave("explicitDiagnosticEquations",3);
   }

   void SimulationBase::solveTrivialEquations(const SolveTiming::Id time)
   {
      // Debug statement
      DebuggerMacro_enter("solveTrivialEquations",3);

      DebuggerMacro_start("Solve trivial",4);
      ProfilerMacro_start(ProfilerMacro::TRIVIALEQUATION);
      this->mTrivialCoordinator.setSolveTime(time);
      this->mTrivialCoordinator.solve(this->mScalarTrivialRange, this->mVectorTrivialRange, this->mScalarVariables, this->mVectorVariables);
      ProfilerMacro_stop(ProfilerMacro::TRIVIALEQUATION);
      DebuggerMacro_stop("Solve trivial t = ",4);

      // Debug statement
      DebuggerMacro_leave("solveTrivialEquations",3);
   }

   void SimulationBase::solveDiagnosticEquations(const SolveTiming::Id time)
   {
      // Debug statement
      DebuggerMacro_enter("solveDiagnosticEquations",3);

      DebuggerMacro_start("Solve diagnostic",4);
      ProfilerMacro_start(ProfilerMacro::DIAGNOSTICEQUATION);
      this->mLinearCoordinator.setSolveTime(time);
      this->mLinearCoordinator.solve(this->mScalarDiagnosticRange, this->mVectorDiagnosticRange, this->mScalarVariables, this->mVectorVariables);
      ProfilerMacro_stop(ProfilerMacro::DIAGNOSTICEQUATION);
      DebuggerMacro_stop("Solve diagnostic t = ",4);

      // Debug statement
      DebuggerMacro_leave("solveDiagnosticEquations",3);
   }
      
   void SimulationBase::initTransformCoordinator(const std::vector<Transform::IntegratorTree>& integratorTree, const std::vector<Transform::ProjectorTree>& projectorTree)
   {
      // Extract the run options for the equation parameters
      std::map<NonDimensional::Id,MHDFloat> runOptions;
      std::vector<NonDimensional::Id>::iterator it;
      std::vector<NonDimensional::Id>  names = this->mspEqParams->ids();
      for(it = names.begin(); it != names.end(); ++it)
      {
         runOptions.insert(std::make_pair(*it, this->mspEqParams->nd(*it)));
      }

      // Initialise the transform coordinator
      Transform::TransformCoordinatorTools::init(this->mTransformCoordinator, this->mspFwdGrouper, this->mspBwdGrouper, integratorTree, projectorTree, this->mspRes, runOptions);
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
      // Loop over all ASCII files added to the simulation control
      SimulationIoControl::ascii_iterator  asciiIt;
      for(asciiIt = this->mSimIoCtrl.beginAscii(); asciiIt != this->mSimIoCtrl.endAscii(); ++asciiIt)
      {
         // Loop over all scalars
         std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>::iterator scalIt;
         for(scalIt = this->mScalarVariables.begin(); scalIt != this->mScalarVariables.end(); scalIt++)
         {
            (*asciiIt)->addScalar((*scalIt));
         }

         // Loop over all vector variables
         std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>::iterator vectIt;
         for(vectIt = this->mVectorVariables.begin(); vectIt != this->mVectorVariables.end(); vectIt++)
         {
            (*asciiIt)->addVector((*vectIt));
         }
      }

      // Loop over all HDF5 files added to the simulation control
      SimulationIoControl::hdf5_iterator  hdf5It;
      for(hdf5It = this->mSimIoCtrl.beginHdf5(); hdf5It != this->mSimIoCtrl.endHdf5(); ++hdf5It)
      {
         // Loop over all scalars
         std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>::iterator scalIt;
         for(scalIt = this->mScalarVariables.begin(); scalIt != this->mScalarVariables.end(); scalIt++)
         {
            (*hdf5It)->addScalar((*scalIt));
         }

         // Loop over all vector variables
         std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>::iterator vectIt;
         for(vectIt = this->mVectorVariables.begin(); vectIt != this->mVectorVariables.end(); vectIt++)
         {
            (*hdf5It)->addVector((*vectIt));
         }

         if((*hdf5It)->space() == Dimensions::Space::PHYSICAL)
         {
            (*hdf5It)->setMesh(this->mTransformCoordinator.mesh());
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
      Equations::Tools::sortByType(this->mScalarEquations, this->mScalarPrognosticRange, this->mScalarDiagnosticRange, this->mScalarTrivialRange);

      // Sort vector equations
      Equations::Tools::sortByType(this->mVectorEquations, this->mVectorPrognosticRange, this->mVectorDiagnosticRange, this->mVectorTrivialRange);

      // Identifiy the solver indexes by analysing the coupling between the equations
      Equations::Tools::identifySolver(this->mScalarPrognosticRange, this->mVectorPrognosticRange);
      Equations::Tools::identifySolver(this->mScalarDiagnosticRange, this->mVectorDiagnosticRange);
      Equations::Tools::identifySolver(this->mScalarTrivialRange, this->mVectorTrivialRange);
   }
}
