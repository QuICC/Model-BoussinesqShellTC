/** \file Simulation.cpp
 *  \brief Source of the high level simulation
 */

// Configuration includes
//
#include "Framework/FrameworkMacro.h"
#include "Profiler/ProfilerMacro.h"
#include "StorageProfiler/StorageProfilerMacro.h"

// System includes
//

// External includes
//

// Class include
//
#include "Simulation/Simulation.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Exceptions/Exception.hpp"
#include "Base/Enums/RuntimeStatus.hpp"
#include "Timers/ExecutionTimer.hpp"
#include "IoTools/Formatter.hpp"

namespace GeoMHDiSCC {

   Simulation::Simulation()
   {
   }

   void Simulation::init()
   {
      // Start timer
      this->mExecutionTimer.start();

      // Make sure to catch raised exception in initialisation steps
      try{
         // Initialise the simulation
         this->initSimulation();

         // Initialise the equations
         this->initEquations();

         // Initialise the variables
         this->initVariables();

         // Setup the equations
         this->setupEquations();

         // Load information from input files (initial states, etc)
         this->initInput();

         // Initialise the timestepper
         this->initTimestepper();

         // Initialise output files (ASCII diagnostics, state files, etc)
         this->initOutput();

         // Setup output files (ASCII diagnostics, state files, etc)
         this->setupOutput();

         // Cleanup unused memory
         this->cleanupSimulation();
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
      FrameworkMacro::synchronise();
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
      while(this->simCtrl().status() == RuntimeStatus::GOON)
      {
         // Compute the nonlinear terms
         this->computeNonlinear();

         // Timestep the equations
         this->timestepEquations();

         // Write the output
         this->writeOutput();

         // Synchronise computation nodes
         FrameworkMacro::synchronise();
      }

      // Stop main loop timing
      this->mExecutionTimer.stop();
      this->mExecutionTimer.update(ExecutionTimer::RUN);

      // Start timer
      this->mExecutionTimer.start();

      // Execute post-run operations
      this->postRun();

      // Synchronise computation nodes
      FrameworkMacro::synchronise();

      // Stop post-run timing
      this->mExecutionTimer.stop();
      this->mExecutionTimer.update(ExecutionTimer::POSTRUN);

      // Synchronise computation nodes
      FrameworkMacro::synchronise();
   }

   void Simulation::finalise()
   {
      // Print simulation control information
      this->simCtrl().printInfo(std::cout);

      // Print execution timer infos
      this->mExecutionTimer.printInfo(std::cout);

      // Print profiling infos (if required)
      ProfilerMacro_printInfo();

      // Print storage profiling infos (if required)
      StorageProfilerMacro_printInfo(std::cout);

      // Finalise the output files
      this->finaliseOutput();

      // Finalise the simulation base
      this->finaliseSimulation();
   }

}
