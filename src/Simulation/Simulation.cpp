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
#include "Simulation/Simulation.hpp"

// Project includes
//
#include "Python/PythonModelWrapper.hpp"
#include "Exceptions/Exception.hpp"
#include "IoTools/Formatter.hpp"
#include "Timers/StageTimer.hpp"
#include "Simulation/SimulationIoTools.hpp"

namespace QuICC {

   Simulation::Simulation()
      : SimulationBase()
   {
   }

   Simulation::~Simulation()
   {
   }

   void Simulation::initAdditionalBase()
   {
      // Get the run configuration
      Array cfgRun = this->mSimIoCtrl.configRun();

      // Set the maximum simulation time
      this->mSimRunCtrl.setMaxSimTime(cfgRun(0));

      // Set the maximum wall time
      this->mSimRunCtrl.setMaxWallTime(cfgRun(1));
   }

   void Simulation::mainRun()
   {
      // Reset the profiler if needed
      ProfilerMacro_printInfo();
      ProfilerMacro_reset();
      ProfilerMacro_init();

      StageTimer::stage("Starting simulation");

      // Start main loop of simulation
      while(this->mSimRunCtrl.status() == RuntimeStatus::GOON)
      {
         // Update equation time
         this->updateEquationTime(this->mTimestepCoordinator.time(), this->mTimestepCoordinator.finishedStep());

         // Compute explicit linear terms
         this->explicitEquations();

         // Update pre calculations required for statistcs output
         if(this->mSimIoCtrl.isStatsUpdateTime())
         {
            SimulationIoTools::updateStatsPre(this->mSimIoCtrl.beginStats(), this->mSimIoCtrl.endStats(), this->mTransformCoordinator);
         }

         // Compute the nonlinear terms
         this->computeNonlinear();

         // Update calculations required for statistcs output
         if(this->mSimIoCtrl.isStatsUpdateTime())
         {
            SimulationIoTools::updateStats(this->mSimIoCtrl.beginStats(), this->mSimIoCtrl.endStats(), this->mTransformCoordinator);
         }

         // Timestep the equations
         this->solveEquations();

         // Update calculations required for statistcs output
         if(this->mSimIoCtrl.isStatsUpdateTime())
         {
            SimulationIoTools::updateStatsPost(this->mSimIoCtrl.beginStats(), this->mSimIoCtrl.endStats(), this->mTransformCoordinator);
            this->mSimIoCtrl.disableStats();
         }

         // Write the output
         this->writeOutput();

         // Synchronise computation nodes
         FrameworkMacro::synchronize();

         // Update simulation run control
         this->mSimRunCtrl.updateCluster(this->mExecutionTimer.queryTime(ExecutionTimer::TOTAL));
      }
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

      // Update equation time
      this->updateEquationTime(this->mDiagnostics.startTime(), false);

      // Initialise all values (solve and nonlinear computations except timestep)
      this->preSolveEquations();

      stage.start("Building timestepper");
      // Print timestepper information
      this->mTimestepCoordinator.printInfo(std::cout);

      // Update CFL condition
      this->mDiagnostics.initialCfl();

      // Synchronise diagnostics
      this->mDiagnostics.synchronize();

      // Init timestepper using clf/100 as starting timestep
      this->mTimestepCoordinator.init(this->mDiagnostics.startTime(), this->mDiagnostics.cfl(), this->mDiagnostics.maxError(), this->mScalarPrognosticRange, this->mVectorPrognosticRange);

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
      stage.start("write initial statistics files");

      // Update calculation required for statistics output
      this->mSimIoCtrl.prepareStats(this->mTimestepCoordinator.time(), this->mTimestepCoordinator.timestep());
      SimulationIoTools::updateStatsPre(this->mSimIoCtrl.beginStats(), this->mSimIoCtrl.endStats(), this->mTransformCoordinator);
      SimulationIoTools::updateStats(this->mSimIoCtrl.beginStats(), this->mSimIoCtrl.endStats(), this->mTransformCoordinator);
      SimulationIoTools::updateStatsPost(this->mSimIoCtrl.beginStats(), this->mSimIoCtrl.endStats(), this->mTransformCoordinator);

      // Write initial statistics output
      this->mSimIoCtrl.writeStats();

      stage.done();
   }

   void Simulation::explicitEquations()
   {
      // Explicit trivial equations
      this->explicitTrivialEquations(ModelOperator::EXPLICIT_LINEAR);

      // Explicit diagnostic equations
      this->explicitDiagnosticEquations(ModelOperator::EXPLICIT_LINEAR);

      // Explicit prognostic equations
      this->explicitPrognosticEquations(ModelOperator::EXPLICIT_LINEAR);
   }

   void Simulation::solveEquations()
   {
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
   }

   void Simulation::explicitPrognosticEquations(const ModelOperator::Id opId)
   {
      ProfilerMacro_start(ProfilerMacro::PROGNOSTICEQUATION);
      this->mTimestepCoordinator.getExplicitInput(opId, this->mScalarPrognosticRange, this->mVectorPrognosticRange, this->mScalarVariables, this->mVectorVariables);
      ProfilerMacro_stop(ProfilerMacro::PROGNOSTICEQUATION);
   }

   void Simulation::solvePrognosticEquations()
   {
      ProfilerMacro_start(ProfilerMacro::PROGNOSTICEQUATION);
      this->mTimestepCoordinator.setSolveTime(SolveTiming::PROGNOSTIC);
      this->mTimestepCoordinator.stepForward(this->mScalarPrognosticRange, this->mVectorPrognosticRange, this->mScalarVariables, this->mVectorVariables);
      ProfilerMacro_stop(ProfilerMacro::PROGNOSTICEQUATION);

      this->conserveAngularMomentum();
   }

   void Simulation::conserveAngularMomentum()
   {
      auto spV = this->mVectorVariables.at(PhysicalNames::VELOCITY);
      Transform::TransformCoordinatorType& coord = this->mTransformCoordinator;

      ArrayZ mom;

      // Dealias toroidal variable data
      coord.communicator().dealiasSpectral(spV->dom(0).total().comp(FieldComponents::Spectral::TOR));

      // Recover dealiased BWD data
      Transform::TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVarTor = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

      // Compute energy integral for first dimension
      coord.transform1D().integrate_volume(mom, rInVarTor.data(), Transform::TransformCoordinatorType::Transform1DType::ProjectorType::VOLUME_PROJ, Transform::TransformCoordinatorType::Transform1DType::IntegratorType::VOLUME_R3);

      const auto& res = *spV->dom(0).spRes();
      MHDFloat lfactor = 0.0;
      MHDFloat factor = 1.0;
      int idx = 0;
      bool hasMOrdering = false;
      if(hasMOrdering)
      {
         // Loop over harmonic order m
         for(int k = 0; k < res.cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int m_ = res.cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
            // m = 0, no factor of two
            if(m_ == 0)
            {
               factor = 1.0;
            } else
            {
               factor = std::sqrt(2.0);
            }

            for(int j = 0; j < res.cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
            {
               int l_ = res.cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);
               lfactor = std::sqrt(16.0*Math::PI/(2.0*l_+1.0));

               if(l_ == 1 && m_ == 0)
               {

                  MHDFloat corZ = (factor*lfactor*mom(idx).real())/std::sqrt(64.0/75.0);
                  MHDComplex curZ = spV->dom(0).perturbation().comp(FieldComponents::Spectral::TOR).point(0,j,k);
                  spV->rDom(0).rPerturbation().rComp(FieldComponents::Spectral::TOR).setPoint(curZ - corZ, 0,j,k);
               } else if(l_ == 1 && m_ == 1)
               {
                  MHDFloat corX = (factor*lfactor*mom(idx).real())/std::sqrt(128.0/75.0);
                  MHDFloat corY = (factor*lfactor*mom(idx).imag())/std::sqrt(128.0/75.0);
                  MHDComplex curXY = spV->dom(0).perturbation().comp(FieldComponents::Spectral::TOR).point(0,j,k);
                  spV->rDom(0).rPerturbation().rComp(FieldComponents::Spectral::TOR).setPoint(curXY - MHDComplex(corX, corY), 0, j, k);
               }

               idx += 1;
            }
         }
      } else
      {
         // Loop over harmonic degree l
         for(int k = 0; k < res.cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int l_ = res.cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
            lfactor = std::sqrt(16.0*Math::PI/(2.0*l_+1.0));
            for(int j = 0; j < res.cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
            {
               int m_ = res.cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k);
               if(l_ == 1 && m_ == 0)
               {
                  MHDFloat corZ = (factor*lfactor*mom(idx).real())/std::sqrt(64.0/75.0);
                  MHDComplex curZ = spV->dom(0).perturbation().comp(FieldComponents::Spectral::TOR).point(0,j,k);
                  spV->rDom(0).rPerturbation().rComp(FieldComponents::Spectral::TOR).setPoint(curZ - corZ, 0,j,k);
               } else if(l_ == 1 && m_ == 1)
               {
                  factor = std::sqrt(2.0);
                  MHDFloat corX = (factor*lfactor*mom(idx).real())/std::sqrt(128.0/75.0);
                  MHDFloat corY = (factor*lfactor*mom(idx).imag())/std::sqrt(128.0/75.0);
                  MHDComplex curXY = spV->dom(0).perturbation().comp(FieldComponents::Spectral::TOR).point(0,j,k);
                  spV->rDom(0).rPerturbation().rComp(FieldComponents::Spectral::TOR).setPoint(curXY - MHDComplex(corX, corY), 0, j, k);
               }

               idx += 1;
            }
         }
      }

      // Free BWD storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVarTor);
   }

   void Simulation::writeOutput()
   {
      ProfilerMacro_start(ProfilerMacro::IO);
      this->mSimIoCtrl.writeStats();

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
   }

   void Simulation::updateEquationTime(const MHDFloat time, const bool finished)
   {
      // Create iterators over scalar equations
      std::vector<Equations::SharedIScalarEquation>::iterator scalEqIt;
      // Create iterators over vector equations
      std::vector<Equations::SharedIVectorEquation>::iterator vectEqIt;

      // Loop over all scalar equations
      for(scalEqIt = this->mScalarEquations.begin(); scalEqIt < this->mScalarEquations.end(); ++scalEqIt)
      {
         (*scalEqIt)->setTime(time, finished);
      }

      // Loop over all vector equations
      for(vectEqIt = this->mVectorEquations.begin(); vectEqIt < this->mVectorEquations.end(); ++vectEqIt)
      {
         (*vectEqIt)->setTime(time, finished);
      }
   }

   void Simulation::postRun()
   {
      this->mSimIoCtrl.writeStats();

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

      stage.done();
      stage.start("write final statistics files");

//      // Update calculation required for statistics output
//      this->mSimIoCtrl.prepareStats(this->mTimestepCoordinator.time(), this->mTimestepCoordinator.timestep());
//      SimulationIoTools::updateStatsPre(this->mSimIoCtrl.beginStats(), this->mSimIoCtrl.endStats(), this->mTransformCoordinator);
//      SimulationIoTools::updateStats(this->mSimIoCtrl.beginStats(), this->mSimIoCtrl.endStats(), this->mTransformCoordinator);
//      SimulationIoTools::updateStatsPost(this->mSimIoCtrl.beginStats(), this->mSimIoCtrl.endStats(), this->mTransformCoordinator);
//
//      // Write final statistics output
//      this->mSimIoCtrl.writeStats();

      stage.done();

      // Synchronise all nodes of simulation
      FrameworkMacro::synchronize();
   }

}
