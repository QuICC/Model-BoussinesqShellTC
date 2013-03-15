/** \file GenerateState.cpp
 *  \brief Simple executable to generate a state file for a model
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
#include "IoVariable/StateFileWriter.hpp"
#include MODELHEADER

typedef GeoMHDiSCC::GEOMHDISCC_RUNSIM_MODEL PModel;

/**
 * @brief Setup and run the simulation
 */
int run()
{
   // Set nCpu for serial run
   int nCpu = 1;

   // Set ID and nCpu in MPI case
   #ifdef GEOMHDISCC_MPI
      // Get MPI size
      int size;
      MPI_Comm_size(MPI_COMM_WORLD, &nCpu);
   #endif //GEOMHDISCC_MPI

   // Setup framework
   GeoMHDiSCC::FrameworkMacro::setup(nCpu);

   // Set type string
   std::string type = PModel::SchemeType::type();

   // Set data regularity
   bool isRegular = PModel::SchemeType::isRegular();

   // Create configuration writer
   GeoMHDiSCC::IoVariable::StateFileWriter state(type, isRegular);

   // Add expected fields
   std::vector<GeoMHDiSCC::PhysicalNames::Id>  ids = PModel::fieldIds();
   std::vector<GeoMHDiSCC::PhysicalNames::Id>::const_iterator  it;
   for(it = ids.begin(); it != ids.end(); ++it)
   {
      state.expect(*it);
   }

   // Set spectral
   GeoMHDiSCC::ArrayI dims(PModel::DIMENSION);
   dims.setConstant(15);
   
   // Create the load splitter
   GeoMHDiSCC::Parallel::LoadSplitter splitter(GeoMHDiSCC::FrameworkMacro::id(), GeoMHDiSCC::FrameworkMacro::nCpu());

   // Initialise the load splitter
   splitter.init<PModel::SchemeType>(dims);

   // Get best splitting resolution object
   std::pair<GeoMHDiSCC::SharedResolution, GeoMHDiSCC::Parallel::SplittingDescription>  best = splitter.bestSplitting();

   // Store the shared resolution object
   GeoMHDiSCC::SharedResolution spRes = best.first;

   // Create scalar variable
   for(it = ids.begin(); it != ids.end(); ++it)
   {
      GeoMHDiSCC::Datatypes::SharedScalarVariableType   spScalar(new GeoMHDiSCC::Datatypes::ScalarVariableType(spRes));
      spScalar->rDom(0).rPerturbation().rData().setRandom();
      spScalar->rDom(0).rPerturbation().rData() *= 1e-2;
      state.addScalar(std::make_pair(*it, spScalar));
   }

   int status = 0;

   // Make sure all information was provided
   if(state.isFull())
   {
      // Initialise state file
      state.init();

      // Write state file
      state.write();

      // Finalize state file
      state.finalize();
   } else
   {
      status = -1;
   }

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
