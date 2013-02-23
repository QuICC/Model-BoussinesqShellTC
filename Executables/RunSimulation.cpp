/** \file RunSimulation.cpp
 *  \brief General executable for the simulation implementations
 */
/// Set the path to the simulation implementation
#define MODELPATH PhysicalModels/GEOMHDISCC_RUNSIM_MODEL.hpp
/// Define small macros allowing to convert to string
#define MAKE_STR_X( _P ) # _P
/// Define small macros allowing to convert to string
#define MAKE_STR( _P ) MAKE_STR_X( _P )
/// Create header include string for the required implementation
#define MODELHEADER MAKE_STR( MODELPATH )
/// Create implementation name macro
#define MODEL mhd::GEOMHDISCC_RUNSIM_MODEL

// Configuration includes
//
#include "Framework/FrameworkMacro.h"

// System includes
//
#include <iostream>

// Project includes
//
#include "Simulation/Simulation.hpp"
#include MODELHEADER

namespace mhd = GeoMHDiSCC;

/**
 * @brief Setup and run the simulation
 */
int run()
{
   // Create simulation
   mhd::SharedSimulation   spSim = MODEL::createSimulation();

   // Exception handling during the initialisation part
   try
   {
      // Initialise the whole simulation
      spSim->init();
   }
   // If exception is thrown, finalise (close files) and return
   catch(int i)
   {
      // Cleanup and close file handles
      spSim->finalize();

      return i;
   }

   // Run the simulation
   spSim->run();

   // Cleanup and close file handles
   spSim->finalize();

   return 0;
}

/**
 * @brief General main, setting up MPI if required
 *
 * The actual program is in run to make sure MPI initialisations
 * are called before anything else and finalisation after destruction
 */
int main(int argc, char* argv[])
{
   // Initilise everything that can't be done inside a class
   mhd::FrameworkMacro::init();

   // Storage for the return code of run
   int code;

   // Compute simulation
   code = run();

   // Finalise everything that can't be done inside a class
   mhd::FrameworkMacro::finalize();

   return code;
}
