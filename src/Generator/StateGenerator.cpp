/** 
 * @file StateGenerator.cpp
 * @brief Source of the high level state generator
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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

// External includes
//

// Class include
//
#include "Generator/StateGenerator.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"
#include "IoTools/Formatter.hpp"

namespace GeoMHDiSCC {

   StateGenerator::StateGenerator()
      : SimulationBase()
   {
   }

   StateGenerator::~StateGenerator()
   {
   }

   void StateGenerator::preRun()
   {
      // Debug statement
      DebuggerMacro_enter("preRun",1);

      // Debug statement
      DebuggerMacro_leave("preRun",1);
   }


   void StateGenerator::mainRun()
   {
      // Debug statement
      DebuggerMacro_enter("mainRun",1);

      // Compute nonlinear terms
      this->computeNonlinear();

      // Solve trivial equations
      this->solveTrivialEquations();

      // Solve diagnostic equations
      this->solveDiagnosticEquations();

      // Synchronise computation nodes
      FrameworkMacro::synchronize();

      // Debug statement
      DebuggerMacro_leave("mainRun",1);
   }

   void StateGenerator::postRun()
   {
      // Debug statement
      DebuggerMacro_enter("postRun",1);

      // Write the output
      this->writeOutput();

      // Debug statement
      DebuggerMacro_leave("postRun",1);
   }

   void StateGenerator::tuneOutput()
   {
      // Debug statement
      DebuggerMacro_enter("tuneOutput",2);

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

      // Debug statement
      DebuggerMacro_leave("tuneOutput",2);
   }

   void StateGenerator::writeOutput()
   {
      // Debug statement
      DebuggerMacro_enter("writeOutput",1);

      // Write final state file (using store time and timestep)
      this->mSimIoCtrl.writeHdf5(-1, -1);

      // Synchronise all nodes of simulation
      FrameworkMacro::synchronize();

      // Debug statement
      DebuggerMacro_leave("writeOutput",1);
   }

}
