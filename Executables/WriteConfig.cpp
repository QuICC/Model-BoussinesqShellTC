/**
 * @file WriteConfig.cpp
 * @brief Simple executable to write a configuration file template for a model 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

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
#include "Python/PythonModelWrapper.hpp"
#include "Exceptions/Exception.hpp"
#include "PhysicalModels/PhysicalModelBase.hpp"
#include "Enums/FieldIds.hpp"
#include "Enums/NonDimensional.hpp"
#include "IoConfig/ConfigurationWriter.hpp"
#include "IoConfig/ConfigParts/PhysicalPart.hpp"
#include "IoConfig/ConfigParts/BoundaryPart.hpp"
#include "IoTools/IdToHuman.hpp"
#include MODELHEADER

typedef QuICC::QUICC_RUNSIM_MODEL PModel;

/**
 * @brief Setup and run the simulation
 */
int run()
{
   // Set dimension
   int dim = PModel::DIMENSION;

   // Set type string
   std::string type = PModel::SchemeType::type();

   // Initialize the python model wrapper
   QuICC::PythonModelWrapper::init();
   QuICC::PythonModelWrapper::import(PModel::PYMODULE);
   QuICC::PythonModelWrapper::createModel(PModel::PYCLASS);

   // Box periodicity
   std::vector<bool> isPeriodicBox = QuICC::PhysicalModelBase::isPeriodicBox();

   // Create configuration writer
   QuICC::IoConfig::ConfigurationWriter writer(dim, isPeriodicBox, type);

   // Create list of field ID strings for boundary conditions
   std::vector<QuICC::PhysicalNames::Id> fields = QuICC::PhysicalModelBase::fieldIds();
   std::vector<QuICC::PhysicalNames::Id>::iterator fIt;
   std::vector<std::string>   bcNames;
   for(fIt = fields.begin(); fIt != fields.end(); ++fIt)
   {
      bcNames.push_back(QuICC::IoTools::IdToHuman::toTag(*fIt));
   }

   // Create list of nondimensional ID strings for physical parameters
   std::vector<QuICC::NonDimensional::Id> params = QuICC::PhysicalModelBase::paramIds();
   std::vector<QuICC::NonDimensional::Id>::iterator pIt;
   std::vector<std::string>   ndNames;
   for(pIt = params.begin(); pIt != params.end(); ++pIt)
   {
      ndNames.push_back(QuICC::IoTools::IdToHuman::toTag(*pIt));
   }

   // Add the physical part
   QuICC::IoConfig::SharedPhysicalPart   spPhys(new QuICC::IoConfig::PhysicalPart(ndNames));
   writer.addPart(QuICC::IoConfig::SimulationBlocks::PHYSICAL, spPhys);

   // Add the boundary part
   QuICC::IoConfig::SharedBoundaryPart   spBound(new QuICC::IoConfig::BoundaryPart(bcNames));
   writer.addPart(QuICC::IoConfig::SimulationBlocks::BOUNDARY, spBound);

   // Initialise writer
   writer.init();

   // Write configuration
   writer.write();

   // Finalise writer
   writer.finalize();

   // Finalize the python model wrapper
   QuICC::PythonModelWrapper::finalize();

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
   QuICC::FrameworkMacro::init();

   // Storage for the return code of run
   int code;

   // Compute simulation
   code = run();

   // Finalise everything that can't be done inside a class
   QuICC::FrameworkMacro::finalize();

   return code;
}
