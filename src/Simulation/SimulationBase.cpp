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
#include "Python/PythonModelWrapper.hpp"
#include "Python/PythonTools.hpp"

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
      // Make sure to catch raised exception in initialisation steps
      try{
         // Initialise the IO system
         this->mSimIoCtrl.init();

         //
         // Extend/modify parameters with automatically computed values
         //
         PyObject *pValue;
         PyObject *pTmp = PythonTools::makeDict(this->mSimIoCtrl.configPhysical());

         // Call model operator Python routine
         PyObject *pArgs = PyTuple_New(1);
         PyTuple_SetItem(pArgs, 0, pTmp);
         PythonModelWrapper::setMethod((char *)"automatic_parameters");
         pValue = PythonModelWrapper::callMethod(pArgs);

         // Create storage
         PythonTools::getDict(this->mSimIoCtrl.rConfigPhysical(), pValue, true);
         Py_DECREF(pValue);
         Py_DECREF(pTmp);
         Py_DECREF(pArgs);

         // Cleanup Python interpreter
         PythonModelWrapper::cleanup();

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
   }

   void SimulationBase::init(const SharedSimulationBoundary spBcs)
   {
      StageTimer stage;

      // Transform projector tree
      std::vector<Transform::TransformTree> projectorTree;

      stage.start("initializing variables");

      // Initialise the variables and set general variable requirements
      RequirementTools::initVariables(projectorTree, this->mScalarVariables, this->mVectorVariables, this->mScalarEquations, this->mVectorEquations, this->mspRes);

      // Transform integrator tree
      std::vector<Transform::TransformTree> integratorTree;

      // Map variables to the equations and set nonlinear requirements
      RequirementTools::mapEquationVariables(integratorTree, this->mScalarEquations, this->mVectorEquations, this->mScalarVariables, this->mVectorVariables, this->mForwardIsNonlinear);

      stage.done();

      // Initialise the transform coordinator
      this->initTransformCoordinator(integratorTree, projectorTree);

      // Initialize Imposed fields
      this->initImposed();

      stage.start("setup equations");

      // Initialise the equations (generate operators, etc)
      this->setupEquations(spBcs);

      // Sort the equations by type: time/solver/trivial
      this->sortEquations();

      stage.done();
      stage.start("setup output files");

      // Setup output files (ASCII diagnostics, state files, etc)
      this->setupOutput();

      // Get timestep information from configuration file
      Array tstep = this->mSimIoCtrl.configTimestepping();

      stage.done();
      stage.start("initializing diagnostics");

      // Initialise the diagnostics
      this->mDiagnostics.init(this->mTransformCoordinator.mesh(), this->mScalarVariables, this->mVectorVariables, tstep);

      // Cleanup IO control
      this->mSimIoCtrl.cleanup();

      stage.done();
   }

   void SimulationBase::initImposed()
   {
      StageTimer stage;

      // Transform projector tree
      std::vector<Transform::TransformTree> projectorTree;

      stage.start("initializing imposed variables");

      // Initialise the variables and set general variable requirements
      RequirementTools::initImposedVariables(projectorTree, this->mImposedScalarVariables, this->mImposedVectorVariables, this->mScalarEquations, this->mVectorEquations, this->mspRes);

      // Transform integrator tree
      std::vector<Transform::TransformTree> integratorTree;

      // Map variables to the equations and set nonlinear requirements
      RequirementTools::mapImposedVariables(this->mScalarEquations, this->mVectorEquations, this->mImposedScalarVariables, this->mImposedVectorVariables);

      stage.done();

      // Extract the run options for the equation parameters
      std::map<NonDimensional::Id,MHDFloat> runOptions;
      std::vector<NonDimensional::Id>::iterator it;
      std::vector<NonDimensional::Id>  names = this->mspEqParams->ids();
      for(it = names.begin(); it != names.end(); ++it)
      {
         runOptions.insert(std::make_pair(*it, this->mspEqParams->nd(*it)));
      }

      Transform::TransformCoordinatorType imposedTransformCoordinator;

      // Initialise the transform coordinator
      Transform::TransformCoordinatorTools::init(imposedTransformCoordinator, this->mspImposedFwdGrouper, this->mspImposedBwdGrouper, integratorTree, projectorTree, this->mspRes, runOptions);

      // Compute physical space values if required
      this->mspImposedBwdGrouper->transform(this->mImposedScalarVariables, this->mImposedVectorVariables, imposedTransformCoordinator);

      this->mspImposedFwdGrouper.reset();
      this->mspImposedBwdGrouper.reset();
   }

   void SimulationBase::run()
   {
      // Final initialisation of the solvers
      this->initSolvers();

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

      // Synchronize computation nodes
      FrameworkMacro::synchronize();
      StageTimer::completed("Simulation initialisation successfull");

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

      // Stop post-run timing
      this->mExecutionTimer.stop();
      this->mExecutionTimer.update(ExecutionTimer::POSTRUN);

      // Synchronise computation nodes
      FrameworkMacro::synchronize();
   }

   void SimulationBase::finalize()
   {
      // Print simulation run infos
      this->mSimRunCtrl.printInfo(std::cout);

      // Print execution timer infos
      this->mExecutionTimer.printInfo(std::cout);

      // Print profiling infos (if required)
      ProfilerMacro_printInfo();

      // Print storage profiling infos (if required)
      StorageProfilerMacro_printInfo();
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

      // Loop over all imposed scalars
      for(scalIt = this->mImposedScalarVariables.begin(); scalIt != this->mImposedScalarVariables.end(); scalIt++)
      {
         spInitFile->addScalar((*scalIt));
      }

      // Loop over all imposed vector variables
      for(vectIt = this->mImposedVectorVariables.begin(); vectIt != this->mImposedVectorVariables.end(); vectIt++)
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
      StageTimer stage;
      stage.start("building trivial solvers");

      // Init trivial solver for trivial equations
      this->mTrivialCoordinator.init(this->mScalarTrivialRange, this->mVectorTrivialRange);

      stage.done();
      stage.start("building diagnostic solvers");

      // Init linear solver for trivial equations
      this->mLinearCoordinator.init(this->mScalarDiagnosticRange, this->mVectorDiagnosticRange);

      stage.done();
   }

   void SimulationBase::computeNonlinear()
   {
      // Compute backward transform
      ProfilerMacro_start(ProfilerMacro::BWDTRANSFORM);
      this->mspBwdGrouper->transform(this->mScalarVariables, this->mVectorVariables, this->mTransformCoordinator);
      ProfilerMacro_stop(ProfilerMacro::BWDTRANSFORM);

      // compute nonlinear interaction and forward transform
      this->mspFwdGrouper->transform(this->mScalarEquations, this->mVectorEquations, this->mTransformCoordinator);
   }

   void SimulationBase::explicitTrivialEquations(const ModelOperator::Id opId)
   {
      ProfilerMacro_start(ProfilerMacro::TRIVIALEQUATION);
      this->mTrivialCoordinator.getExplicitInput(opId, this->mScalarTrivialRange, this->mVectorTrivialRange, this->mScalarVariables, this->mVectorVariables);
      ProfilerMacro_stop(ProfilerMacro::TRIVIALEQUATION);
   }

   void SimulationBase::explicitDiagnosticEquations(const ModelOperator::Id opId)
   {
      ProfilerMacro_start(ProfilerMacro::DIAGNOSTICEQUATION);
      this->mLinearCoordinator.getExplicitInput(opId, this->mScalarDiagnosticRange, this->mVectorDiagnosticRange, this->mScalarVariables, this->mVectorVariables);
      ProfilerMacro_stop(ProfilerMacro::DIAGNOSTICEQUATION);
   }

   void SimulationBase::solveTrivialEquations(const SolveTiming::Id time)
   {
      ProfilerMacro_start(ProfilerMacro::TRIVIALEQUATION);
      this->mTrivialCoordinator.setSolveTime(time);
      this->mTrivialCoordinator.solve(this->mScalarTrivialRange, this->mVectorTrivialRange, this->mScalarVariables, this->mVectorVariables);
      ProfilerMacro_stop(ProfilerMacro::TRIVIALEQUATION);
   }

   void SimulationBase::solveDiagnosticEquations(const SolveTiming::Id time)
   {
      ProfilerMacro_start(ProfilerMacro::DIAGNOSTICEQUATION);
      this->mLinearCoordinator.setSolveTime(time);
      this->mLinearCoordinator.solve(this->mScalarDiagnosticRange, this->mVectorDiagnosticRange, this->mScalarVariables, this->mVectorVariables);
      ProfilerMacro_stop(ProfilerMacro::DIAGNOSTICEQUATION);
   }
      
   void SimulationBase::initTransformCoordinator(const std::vector<Transform::TransformTree>& integratorTree, const std::vector<Transform::TransformTree>& projectorTree)
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

         // Loop over all imposed scalars
         for(scalIt = this->mImposedScalarVariables.begin(); scalIt != this->mImposedScalarVariables.end(); scalIt++)
         {
            (*hdf5It)->addScalar((*scalIt));
         }

         // Loop over all imposed vector variables
         for(vectIt = this->mImposedVectorVariables.begin(); vectIt != this->mImposedVectorVariables.end(); vectIt++)
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
      Equations::Tools::sortByType(this->mScalarEquations, this->mScalarPrognosticRange, this->mScalarDiagnosticRange, this->mScalarTrivialRange, this->mScalarWrapperRange);

      // Sort vector equations
      Equations::Tools::sortByType(this->mVectorEquations, this->mVectorPrognosticRange, this->mVectorDiagnosticRange, this->mVectorTrivialRange, this->mVectorWrapperRange);

      // Identifiy the solver indexes by analysing the coupling between the equations
      Equations::Tools::identifySolver(this->mScalarPrognosticRange, this->mVectorPrognosticRange);
      Equations::Tools::identifySolver(this->mScalarDiagnosticRange, this->mVectorDiagnosticRange);
      Equations::Tools::identifySolver(this->mScalarTrivialRange, this->mVectorTrivialRange);
   }
}
