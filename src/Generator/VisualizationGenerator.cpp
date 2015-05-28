/** 
 * @file VisualizationGenerator.cpp
 * @brief Source of the high level state generator
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

// External includes
//

// Class include
//
#include "Generator/VisualizationGenerator.hpp"

// Project includes
//
#include "Python/PythonModelWrapper.hpp"
#include "Exceptions/Exception.hpp"
#include "IoTools/Formatter.hpp"

namespace GeoMHDiSCC {

   VisualizationGenerator::VisualizationGenerator()
      : SimulationBase()
   {
   }

   VisualizationGenerator::~VisualizationGenerator()
   {
   }

   void VisualizationGenerator::preRun()
   {
      // Finalizing the Python model wrapper
      PythonModelWrapper::finalize();
   }


   void VisualizationGenerator::mainRun()
   {
      // Solve the trivial equations
      this->solveTrivialEquations(SolveTiming::AFTER);

      // Compute nonlinear terms
      this->computeNonlinear();

      // Synchronise computation nodes
      FrameworkMacro::synchronize();
   }

   void VisualizationGenerator::postRun()
   {
      // Synchronise all nodes of simulation
      FrameworkMacro::synchronize();

      // Write the output
      this->writeOutput();

      // Synchronise all nodes of simulation
      FrameworkMacro::synchronize();
   }

   void VisualizationGenerator::tuneOutput()
   {
   }

   void VisualizationGenerator::writeOutput()
   {
      // Write final state file (using store time and timestep)
      this->mSimIoCtrl.writeHdf5(this->mDiagnostics.startTime(), this->mDiagnostics.startTimestep());

      // Synchronise all nodes of simulation
      FrameworkMacro::synchronize();
   }

   void VisualizationGenerator::tuneInitialState(IoVariable::SharedStateFileReader spInitFile)
   {
      // Get time from initial state
      MHDFloat time = spInitFile->time();

      // Get timestep from initial state
      MHDFloat timestep = spInitFile->timestep();

      // Loop over all files added to the simulation control
      SimulationIoControl::hdf5_iterator  fIt;
      for(fIt = this->mSimIoCtrl.beginHdf5(); fIt != this->mSimIoCtrl.endHdf5(); ++fIt)
      {
         (*fIt)->setSimTime(time, timestep);
      }
   }

}
