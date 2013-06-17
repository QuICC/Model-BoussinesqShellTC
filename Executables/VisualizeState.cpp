/** \file VisualizeState.cpp
 *  \brief Simple general executable to visualize a state file for a model
 */
/// Set the path to the simulation implementation
#define MODELPATH PhysicalModels/GEOMHDISCC_RUNSIM_MODEL.hpp
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
#include "Generator/VisualizationGenerator.hpp"
#include "PhysicalModels/VisualizationGeneratorFactory.hpp"
#include MODELHEADER

/**
 * @brief Setup and run the simulation
 */
int run()
{
   int status = 0;

   // Create simulation
   GeoMHDiSCC::SharedVisualizationGenerator   spVis;

   // Exception handling during the initialisation part
   try
   {
      // Create state generator
      spVis = GeoMHDiSCC::VisualizationGeneratorFactory<GeoMHDiSCC::GEOMHDISCC_RUNSIM_MODEL>::createVisualization();
   }

   // If exception is thrown, finalise (close files) and return
   catch(GeoMHDiSCC::Exception& e)
   {
      std::cerr << e.what() << std::endl;

      status = -1;
   }

   if(status == 0)
   {
      // Run the simulation
      spVis->run();
   }

   // Cleanup and close file handles
   spVis->finalize();

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
   GeoMHDiSCC::FrameworkMacro::init();

   // Storage for the return code of run
   int code;

   // Compute simulation
   code = run();

   // Finalise everything that can't be done inside a class
   GeoMHDiSCC::FrameworkMacro::finalize();

   return code;
}
