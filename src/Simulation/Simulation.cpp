/** \file Simulation.cpp
 *  \brief Source of the high level simulation
 */

// Configuration includes
//
#include "Profiler/ProfilerMacro.h"
#include "StorageProfiler/StorageProfilerMacro.h"
#include "Framework/FrameworkMacro.h"

// System includes
//

// External includes
//

// Class include
//
#include "Simulation/Simulation.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"
#include "IoTools/Formatter.hpp"

namespace GeoMHDiSCC {

   Simulation::Simulation()
      : mExecutionTimer(), mSimRunCtrl()
   {
   }

   Simulation::~Simulation()
   {
   }

   void Simulation::init()
   {
      // Start timer
      this->mExecutionTimer.start();

      // Make sure to catch raised exception in initialisation steps
      try{
      }
      catch(Exception &e)
      {
         std::cout << e.what() << std::endl;

         throw -1;
      }

      // Stop timer and update initialisation time
      this->mExecutionTimer.stop();
      this->mExecutionTimer.update(ExecutionTimer::INIT);

      // Print message to signal successful completion of initialisation step
      if(FrameworkMacro::allowsIO())
      {
         IoTools::Formatter::printLine(std::cout, '-');
         IoTools::Formatter::printCentered(std::cout, "Simulation initialisation successfull", '*');
         IoTools::Formatter::printLine(std::cout, '-');
         IoTools::Formatter::printNewline(std::cout);
      }

      // Make sure nodes are synchronised after initialisation
      FrameworkMacro::synchronize();
   }

   void Simulation::run()
   {
      // Start timer
      this->mExecutionTimer.start();

      // Execute pre-run steps
      this->preRun();

      // Stop pre-run timing
      this->mExecutionTimer.stop();
      this->mExecutionTimer.update(ExecutionTimer::PRERUN);

      // Start timer
      this->mExecutionTimer.start();

      // Start main loop of simulation
      while(this->mSimRunCtrl.status() == Runtime::Status::GOON)
      {
         // Compute the nonlinear terms
         this->computeNonlinear();

         // Timestep the equations
         this->timestepEquations();

         // Write the output
         this->writeOutput();

         // Synchronise computation nodes
         FrameworkMacro::synchronize();
      }

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
   }

   void Simulation::finalize()
   {
      // Print execution timer infos
      this->mExecutionTimer.printInfo(std::cout);

      // Print profiling infos (if required)
      ProfilerMacro_printInfo();

      // Print storage profiling infos (if required)
      StorageProfilerMacro_printInfo();
   }

   void Simulation::addEquation(int spEq)//SharedIScalarEquation spEq)
   {
   }

   void Simulation::addEquation(double spEq)//SharedIVectorEquation spEq)
   {
   }

   void Simulation::setConfigurationFile(IoConfig::SharedConfigurationReader spCfgFile)
   {
   }

   void Simulation::setInitialStateFile(int spInitFile)//SharedStateFile spInitFile)
   {
   }

   void Simulation::addOutputFile(int spOutFile)//SharedAscii spOutFile)
   {
   }

   void Simulation::addOutputFile(double spOutFile)//SharedHdf5 spOutFile)
   {
   }

   void Simulation::preRun()
   {
      // Print message to signal successful completion of initialisation step
      if(FrameworkMacro::allowsIO())
      {
         IoTools::Formatter::printLine(std::cout, '-');
         IoTools::Formatter::printCentered(std::cout, "... Starting simulation ...", '*');
         IoTools::Formatter::printLine(std::cout, '-');
         IoTools::Formatter::printNewline(std::cout);
      }

      // Write initial ASCII output
      this->mSimIoCtrl.writeAscii();

      // Write initial state file
      this->mSimIoCtrl.writeHdf5();

      // Synchronise all nodes of simulation
      FrameworkMacro::synchronize();
   }

   void Simulation::computeNonlinear()
   {
      // Compute backward transform
      ProfilerMacro_start(ProfilerMacro::BWDTRANSFORM);
      //this->mspBwdGrouper->transform(this->mScalarVariables, this->mVectorVariables, this->transformCoordinator());
      ProfilerMacro_stop(ProfilerMacro::BWDTRANSFORM);

      // compute nonlinear interaction and forward transform
      //this->mspFwdGrouper->transform(this->mScalarEquations, this->mVectorEquations, this->transformCoordinator());
   }

   void Simulation::timestepEquations()
   {
      ProfilerMacro_start(ProfilerMacro::TIMESTEP);
      //this->mTimestepper.stepForward(this->mScalarEquations, this->mVectorEquations);
      ProfilerMacro_stop(ProfilerMacro::TIMESTEP);

      ProfilerMacro_start(ProfilerMacro::CONTROL);
      if(this->mTimestepper.finishedStep())
      {
         // Update timestepper
         this->mTimestepper.update();
      
         // Update simulation control
         this->mSimRunCtrl.update();
      }
      ProfilerMacro_stop(ProfilerMacro::CONTROL);
   }

   void Simulation::writeOutput()
   {
      ProfilerMacro_start(ProfilerMacro::IO);
      if(this->mTimestepper.finishedStep() && this->mSimRunCtrl.doIO())
      {
         // Write initial ASCII output
         this->mSimIoCtrl.writeAscii();
      
         // Write initial state file
         this->mSimIoCtrl.writeHdf5();
      }
      ProfilerMacro_stop(ProfilerMacro::IO);
   }

   void Simulation::postRun()
   {
      // Write final ASCII output
      this->mSimIoCtrl.writeAscii();

      // Write final state file
      this->mSimIoCtrl.writeHdf5();

      // Synchronise all nodes of simulation
      FrameworkMacro::synchronize();
   }

}
