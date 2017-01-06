/**
 * @file RunSimulation.cpp
 * @brief General executable for the simulation implementations 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */
//#include <Python.h>

/// Set the path to the simulation implementation
#define MODELPATH PhysicalModels/QUICC_RUNSIM_MODEL.hpp
/// Define small macros allowing to convert to string
#define MAKE_STR_X( _P ) # _P
/// Define small macros allowing to convert to string
#define MAKE_STR( _P ) MAKE_STR_X( _P )
/// Create header include string for the required implementation
#define MODELHEADER MAKE_STR( MODELPATH )

// Configuration includes
//
#include "Framework/FrameworkMacro.h"

// System includes
//
#include <iostream>

// Project includes
//
#include "Exceptions/Exception.hpp"
#include "Simulation/Simulation.hpp"
#include "PhysicalModels/ModelFactory.hpp"
#include MODELHEADER

/**
 * @brief Setup and run the simulation
 */
int run()
{
   int status = 0;

   // Create simulation
   QuICC::SharedSimulation   spSim;

   // Exception handling during the initialisation part
   try
   {
      // Create simulation
      spSim = QuICC::ModelFactory<QuICC::QUICC_RUNSIM_MODEL>::createSimulation();
   }

   // If exception is thrown, finalise (close files) and return
   catch(QuICC::Exception& e)
   {
      std::cerr << e.what() << std::endl;

      status = -1;
   }

   if(status == 0)
   {
      // Run the simulation
      spSim->run();
   }

   // Cleanup and close file handles
   spSim->finalize();

   return status;
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
   QuICC::FrameworkMacro::init();

   // Storage for the return code of run
   int code;

   // Compute simulation
   code = run();

   // Finalise everything that can't be done inside a class
   QuICC::FrameworkMacro::finalize();

   return code;
}
