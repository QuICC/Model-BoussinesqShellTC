/**
 * @file WriteConfig.cpp
 * @brief Simple executable to write a configuration file template for a model 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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
#include "PhysicalModels/PhysicalModelBase.hpp"
#include "Enums/FieldIds.hpp"
#include "Enums/NonDimensional.hpp"
#include "IoConfig/ConfigurationWriter.hpp"
#include "IoConfig/ConfigParts/PhysicalPart.hpp"
#include "IoConfig/ConfigParts/BoundaryPart.hpp"
#include "IoTools/IdToHuman.hpp"
#include MODELHEADER

typedef GeoMHDiSCC::GEOMHDISCC_RUNSIM_MODEL PModel;

/**
 * @brief Setup and run the simulation
 */
int run()
{
   // Set dimension
   int dim = PModel::DIMENSION;

   // Set type string
   std::string type = PModel::SchemeType::type();

   // Box periodicity
   std::vector<bool> isPeriodicBox = GeoMHDiSCC::PhysicalModelBase::isPeriodicBox(PModel::PYNAME);

   // Create configuration writer
   GeoMHDiSCC::IoConfig::ConfigurationWriter writer(dim, isPeriodicBox, type);

   // Create list of field ID strings for boundary conditions
   std::vector<GeoMHDiSCC::PhysicalNames::Id> fields = GeoMHDiSCC::PhysicalModelBase::fieldIds(PModel::PYNAME);
   std::vector<GeoMHDiSCC::PhysicalNames::Id>::iterator fIt;
   std::vector<std::string>   bcNames;
   for(fIt = fields.begin(); fIt != fields.end(); ++fIt)
   {
      bcNames.push_back(GeoMHDiSCC::IoTools::IdToHuman::toTag(*fIt));
   }

   // Create list of nondimensional ID strings for physical parameters
   std::vector<GeoMHDiSCC::NonDimensional::Id> params = GeoMHDiSCC::PhysicalModelBase::paramIds(PModel::PYNAME);
   std::vector<GeoMHDiSCC::NonDimensional::Id>::iterator pIt;
   std::vector<std::string>   ndNames;
   for(pIt = params.begin(); pIt != params.end(); ++pIt)
   {
      ndNames.push_back(GeoMHDiSCC::IoTools::IdToHuman::toTag(*pIt));
   }

   // Add the physical part
   GeoMHDiSCC::IoConfig::SharedPhysicalPart   spPhys(new GeoMHDiSCC::IoConfig::PhysicalPart(ndNames));
   writer.addPart(GeoMHDiSCC::IoConfig::SimulationBlocks::PHYSICAL, spPhys);

   // Add the boundary part
   GeoMHDiSCC::IoConfig::SharedBoundaryPart   spBound(new GeoMHDiSCC::IoConfig::BoundaryPart(bcNames));
   writer.addPart(GeoMHDiSCC::IoConfig::SimulationBlocks::BOUNDARY, spBound);

   // Initialise writer
   writer.init();

   // Write configuration
   writer.write();

   // Finalise writer
   writer.finalize();

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
   GeoMHDiSCC::FrameworkMacro::init();

   // Storage for the return code of run
   int code;

   // Compute simulation
   code = run();

   // Finalise everything that can't be done inside a class
   GeoMHDiSCC::FrameworkMacro::finalize();

   return code;
}
