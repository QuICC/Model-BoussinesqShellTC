/** 
 * @file StateGenerator.cpp
 * @brief Source of the high level state generator
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

// External includes
//

// Class include
//
#include "Generator/StateGenerator.hpp"

// Project includes
//
#include "Python/PythonModelWrapper.hpp"
#include "Exceptions/Exception.hpp"
#include "IoTools/Formatter.hpp"

namespace GeoMHDiSCC {

   StateGenerator::StateGenerator()
      : SimulationBase()
   {
      this->mForwardIsNonlinear = false;
   }

   StateGenerator::~StateGenerator()
   {
   }

   void StateGenerator::preRun()
   {
      // Finalizing the Python model wrapper
      PythonModelWrapper::finalize();
   }

   void StateGenerator::mainRun()
   {
      // Solve trivial equations
      this->solveTrivialEquations(SolveTiming::BEFORE);

      // Compute nonlinear terms
      this->computeNonlinear();

/// \mhdBug Problem with equations for generating exact states
      // Solve trivial equations
      this->solveTrivialEquations(SolveTiming::AFTER);

      // Solve diagnostic equations
//      this->solveDiagnosticEquations(SolveTiming::AFTER);

      // Synchronise computation nodes
      FrameworkMacro::synchronize();
   }

   void StateGenerator::postRun()
   {
      // Synchronise all nodes of simulation
      FrameworkMacro::synchronize();

      // Write the output
      this->writeOutput();

      // Synchronise all nodes of simulation
      FrameworkMacro::synchronize();
   }

   void StateGenerator::tuneOutput()
   {
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

      // Loop over all files added to the simulation control
      SimulationIoControl::hdf5_iterator  fIt;
      for(fIt = this->mSimIoCtrl.beginHdf5(); fIt != this->mSimIoCtrl.endHdf5(); ++fIt)
      {
         (*fIt)->setSimTime(tstep(0), tstep(1));
      }
   }

   void StateGenerator::writeOutput()
   {
      // Write final state file (using store time and timestep)
      this->mSimIoCtrl.writeHdf5(-1, -1);

      // Synchronise all nodes of simulation
      FrameworkMacro::synchronize();
   }

}
