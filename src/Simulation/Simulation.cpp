/** 
 * @file Simulation.cpp
 * @brief Source of the high level simulation
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// First includes
//
#include "Python/PythonHeader.hpp"

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
#include "Simulation/Simulation.hpp"

// Project includes
//
#include "Python/PythonModelWrapper.hpp"
#include "Exceptions/Exception.hpp"
#include "IoTools/Formatter.hpp"
#include "Timers/StageTimer.hpp"
#include "Simulation/SimulationIoTools.hpp"

namespace GeoMHDiSCC {

   Simulation::Simulation()
      : SimulationBase()
   {
   }

   Simulation::~Simulation()
   {
   }

   void Simulation::initAdditionalBase()
   {
      // Debug statement
      DebuggerMacro_enter("initAdditionalBase",1);

      // Get the run configuration
      Array cfgRun = this->mSimIoCtrl.configRun();

      // Set the maximum simulation time
      this->mSimRunCtrl.setMaxSimTime(cfgRun(0));

      // Set the maximum wall time
      this->mSimRunCtrl.setMaxWallTime(cfgRun(1));

      // Debug statement
      DebuggerMacro_leave("initAdditionalBase",1);
   }

   void Simulation::mainRun()
   {
      // Debug statement
      DebuggerMacro_enter("mainRun",1);

      StageTimer::stage("Starting simulation");

      // Start main loop of simulation
      while(this->mSimRunCtrl.status() == RuntimeStatus::GOON)
      {
         // Compute explicit linear terms
         this->explicitEquations();

         // Compute the nonlinear terms
         this->computeNonlinear();

         // Timestep the equations
         this->solveEquations();

         // Write the output
         this->writeOutput();

         // Synchronise computation nodes
         FrameworkMacro::synchronize();
      
         // Update simulation run control
         this->mSimRunCtrl.updateCluster(this->mExecutionTimer.queryTime(ExecutionTimer::TOTAL));
      }

      // Debug statement
      DebuggerMacro_leave("mainRun",1);
   }

   void Simulation::preSolveEquations()
   {  
      StageTimer stage;
      stage.start("initializing fields");

      /// \mhdBug This is not sufficient to recover all fields from previous computation

      // Solve diagnostic equations
      this->explicitDiagnosticEquations(ModelOperator::EXPLICIT_NEXTSTEP);
      this->solveDiagnosticEquations(SolveTiming::AFTER);

      // Solve trivial equations
      this->explicitTrivialEquations(ModelOperator::EXPLICIT_NEXTSTEP);
      this->solveTrivialEquations(SolveTiming::AFTER);

      // Compute physical values
      this->mspBwdGrouper->transform(this->mScalarVariables, this->mVectorVariables, this->mTransformCoordinator);

      // Only compute forward transform for diagnostic and trivial equations
      std::vector<Equations::SharedIScalarEquation> scalEqs;
      std::vector<Equations::SharedIVectorEquation> vectEqs;
      for(ScalarEquation_iterator sIt = this->mScalarDiagnosticRange.first; sIt != this->mScalarDiagnosticRange.second; ++sIt)
      {
         scalEqs.push_back(*sIt);
      }
      for(ScalarEquation_iterator sIt = this->mScalarTrivialRange.first; sIt != this->mScalarTrivialRange.second; ++sIt)
      {
         scalEqs.push_back(*sIt);
      }
      for(VectorEquation_iterator vIt = this->mVectorDiagnosticRange.first; vIt != this->mVectorDiagnosticRange.second; ++vIt)
      {
         vectEqs.push_back(*vIt);
      }
      for(VectorEquation_iterator vIt = this->mVectorTrivialRange.first; vIt != this->mVectorTrivialRange.second; ++vIt)
      {
         vectEqs.push_back(*vIt);
      }
      this->mspFwdGrouper->transform(scalEqs, vectEqs, this->mTransformCoordinator);

      // Solve diagnostic equations
      this->explicitDiagnosticEquations(ModelOperator::EXPLICIT_NONLINEAR);
      this->solveDiagnosticEquations(SolveTiming::BEFORE);

      // Solve trivial equations
      this->explicitTrivialEquations(ModelOperator::EXPLICIT_NONLINEAR);
      this->solveTrivialEquations(SolveTiming::BEFORE);

      // Solve diagnostic equations
      this->explicitDiagnosticEquations(ModelOperator::EXPLICIT_NEXTSTEP);
      this->solveDiagnosticEquations(SolveTiming::AFTER);

      // Solve trivial equations
      this->explicitTrivialEquations(ModelOperator::EXPLICIT_NEXTSTEP);
      this->solveTrivialEquations(SolveTiming::AFTER);

      stage.done();

      // Synchronise all nodes of simulation
      FrameworkMacro::synchronize();
   }

   void Simulation::preRun()
   {
      StageTimer stage;

      // Initialise all values (solve and nonlinear computations except timestep)
      this->preSolveEquations();

      stage.start("Building timestepper");

      // Update CFL condition
      this->mDiagnostics.initialCfl();

      // Synchronise diagnostics
      this->mDiagnostics.synchronize();

      // Init timestepper using clf/100 as starting timestep
      this->mTimestepCoordinator.init(this->mDiagnostics.startTime(), this->mDiagnostics.cfl(), this->mScalarPrognosticRange, this->mVectorPrognosticRange);

      // Finalizing the Python model wrapper
      PythonModelWrapper::finalize();

      stage.done();
      stage.start("write initial ASCII files");

      // Update heavy calculation required for ASCII output
      SimulationIoTools::updateHeavyAscii(this->mSimIoCtrl.beginAscii(), this->mSimIoCtrl.endAscii(), this->mTransformCoordinator);

      // Write initial ASCII output
      this->mSimIoCtrl.writeAscii(this->mTimestepCoordinator.time(), this->mTimestepCoordinator.timestep());

      stage.done();
      stage.start("write initial HDF5 files");

      // Write initial state file
      this->mSimIoCtrl.writeHdf5(this->mTimestepCoordinator.time(), this->mTimestepCoordinator.timestep());

      stage.done();
   }

   void Simulation::explicitEquations()
   {
      // Debug statement
      DebuggerMacro_enter("explicitEquations",2);

      // Explicit trivial equations
      this->explicitTrivialEquations(ModelOperator::EXPLICIT_LINEAR);

      // Explicit diagnostic equations
      this->explicitDiagnosticEquations(ModelOperator::EXPLICIT_LINEAR);

      // Explicit prognostic equations
      this->explicitPrognosticEquations(ModelOperator::EXPLICIT_LINEAR);

      // Debug statement
      DebuggerMacro_leave("explicitEquations",2);
   }

   void Simulation::solveEquations()
   {
      // Debug statement
      DebuggerMacro_enter("solveEquations",2);

      // Solve trivial equations
      this->explicitTrivialEquations(ModelOperator::EXPLICIT_NONLINEAR);
      this->solveTrivialEquations(SolveTiming::BEFORE);

      // Solve diagnostic equations
      this->explicitDiagnosticEquations(ModelOperator::EXPLICIT_NONLINEAR);
      this->solveDiagnosticEquations(SolveTiming::BEFORE);

      // Solve prognostic equations (timestep)
      this->explicitPrognosticEquations(ModelOperator::EXPLICIT_NONLINEAR);
      this->solvePrognosticEquations();

      // Solve diagnostic equations
      this->explicitDiagnosticEquations(ModelOperator::EXPLICIT_NEXTSTEP);
      this->solveDiagnosticEquations(SolveTiming::AFTER);

      // Solve trivial equations
      this->explicitTrivialEquations(ModelOperator::EXPLICIT_NEXTSTEP);
      this->solveTrivialEquations(SolveTiming::AFTER);

      // Update conditions at the end of timestep
      ProfilerMacro_start(ProfilerMacro::CONTROL);
      if(this->mTimestepCoordinator.finishedStep())
      {
         // Update timestepper
         this->mTimestepCoordinator.update();

         // Update CFL condition
         this->mDiagnostics.updateCfl();

         // Synchronise diagnostics
         this->mDiagnostics.synchronize();

         // Adapt timestepper time step
         this->mTimestepCoordinator.adaptTimestep(this->mDiagnostics.cfl(), this->mScalarPrognosticRange, this->mVectorPrognosticRange);
      
         // Update simulation run control
         this->mSimRunCtrl.updateSimulation(this->mTimestepCoordinator.time(), this->mTimestepCoordinator.timestep());
      
         // Update simulation IO control
         this->mSimIoCtrl.update();
      }
      ProfilerMacro_stop(ProfilerMacro::CONTROL);

      // Debug statement
      DebuggerMacro_leave("solveEquations",2);
   }

   void Simulation::explicitPrognosticEquations(const ModelOperator::Id opId)
   {
      // Debug statement
      DebuggerMacro_enter("explicitPrognostic",3);

      DebuggerMacro_start("Explicit prognostic",4);
      ProfilerMacro_start(ProfilerMacro::PROGNOSTICEQUATION);
      this->mTimestepCoordinator.getExplicitInput(opId, this->mScalarPrognosticRange, this->mVectorPrognosticRange, this->mScalarVariables, this->mVectorVariables);
      ProfilerMacro_stop(ProfilerMacro::PROGNOSTICEQUATION);
      DebuggerMacro_stop("Explicit prognostic t = ",4);

      // Debug statement
      DebuggerMacro_leave("explicitPrognostic",3);
   }

   void Simulation::solvePrognosticEquations()
   {
      // Debug statement
      DebuggerMacro_enter("solvePrognostic",3);

      DebuggerMacro_start("Solve prognostic",4);
      ProfilerMacro_start(ProfilerMacro::PROGNOSTICEQUATION);
      this->mTimestepCoordinator.setSolveTime(SolveTiming::PROGNOSTIC);
      this->mTimestepCoordinator.stepForward(this->mScalarPrognosticRange, this->mVectorPrognosticRange, this->mScalarVariables, this->mVectorVariables);
      ProfilerMacro_stop(ProfilerMacro::PROGNOSTICEQUATION);
      DebuggerMacro_stop("Solve prognostic t = ",4);

      // Debug statement
      DebuggerMacro_leave("solvePrognostic",3);
   }

   void Simulation::writeOutput()
   {
      // Debug statement
      DebuggerMacro_enter("writeOutput",2);

      ProfilerMacro_start(ProfilerMacro::IO);
      if(this->mTimestepCoordinator.finishedStep())
      {
         if(this->mSimIoCtrl.isAsciiTime())
         {
            // Update heavy calculation required for ASCII output
            SimulationIoTools::updateHeavyAscii(this->mSimIoCtrl.beginAscii(), this->mSimIoCtrl.endAscii(), this->mTransformCoordinator);
         }

         // Write initial ASCII and HDF5 output files if applicable
         this->mSimIoCtrl.writeFiles(this->mTimestepCoordinator.time(), this->mTimestepCoordinator.timestep());
      }
      ProfilerMacro_stop(ProfilerMacro::IO);

      // Debug statement
      DebuggerMacro_leave("writeOutput",2);
   }

   void Simulation::postRun()
   {
      StageTimer::stage("Post simulation");
      StageTimer  stage;

      // Synchronise all nodes of simulation
      FrameworkMacro::synchronize();

      stage.start("write final ASCII files");

      // Update heavy calculation required for ASCII output
      SimulationIoTools::updateHeavyAscii(this->mSimIoCtrl.beginAscii(), this->mSimIoCtrl.endAscii(), this->mTransformCoordinator);

      // Write final ASCII output
      this->mSimIoCtrl.writeAscii(this->mTimestepCoordinator.time(), this->mTimestepCoordinator.timestep());

      stage.done();
      stage.start("write final HDF5 files");

      // Write final state file
      this->mSimIoCtrl.writeHdf5(this->mTimestepCoordinator.time(), this->mTimestepCoordinator.timestep());

      // Synchronise all nodes of simulation
      FrameworkMacro::synchronize();

      stage.done();
   }

}
