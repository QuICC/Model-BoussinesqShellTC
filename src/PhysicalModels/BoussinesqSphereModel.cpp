/** \file BoussinesqSphereModel.cpp
 *  \brief Source of the Boussinesq sphere physical model
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "PhysicalModels/BoussinesqSphereModel.hpp"

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"

namespace GeoMHDiSCC {

   std::vector<PhysicalNames::Id> BoussinesqSphereModel::fieldIds()
   {
      // Create storage
      std::vector<PhysicalNames::Id> ids;

      // Add temperature
      ids.push_back(PhysicalNames::TEMPERATURE);

      // Add temperature
      ids.push_back(PhysicalNames::VELOCITY);

      return ids;
   }

   std::vector<NonDimensional::Id> BoussinesqSphereModel::paramIds()
   {
      // Create storage
      std::vector<NonDimensional::Id> ids;

      return ids;
   }

   std::vector<bool> BoussinesqSphereModel::isPeriodicBox()
   {
      std::vector<bool> box;

      // X direction is not periodic box
      box.push_back(false);

      // Y direction is periodic box
      box.push_back(false);

      // Z direction is not periodic box
      box.push_back(false);

      return box;
   }

   void BoussinesqSphereModel::addEquations(SharedSimulation spSim)
   {
      throw Exception("Not implemented yet!");
   }

   void BoussinesqSphereModel::addStates(SharedStateGenerator spGen)
   {
      throw Exception("Not implemented yet!");
   }

   void BoussinesqSphereModel::addVisualizers(SharedVisualizationGenerator spVis)
   {
      throw Exception("Not implemented yet!");
   }

   void BoussinesqSphereModel::setVisualizationState(SharedVisualizationGenerator spVis)
   {
      throw Exception("Not implemented yet!");
   }

   void BoussinesqSphereModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
      
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
   }

   void BoussinesqSphereModel::addHdf5OutputFiles(SharedSimulation spSim)
   {
      // Field IDs iterator
      std::vector<GeoMHDiSCC::PhysicalNames::Id>::const_iterator  it;
      std::vector<GeoMHDiSCC::PhysicalNames::Id> ids = BoussinesqSphereModel::fieldIds();

      // Create and add state file to IO
      IoVariable::SharedStateFileWriter spState(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      for(it = ids.begin(); it != ids.end(); ++it)
      {
         spState->expect(*it);
      }
      spSim->addOutputFile(spState);
   }

   SharedSimulationBoundary BoussinesqSphereModel::createBoundary(const std::map<std::string,int>& bcIds)
   {
      // Create shared simulation boundary
      SharedSimulationBoundary  spBcs(new SimulationBoundary());

      // Storage for the dimension ID
      Dimensions::Simulation::Id dimId;

      // Create equation and field keys
      SpectralFieldId eqId;
      SpectralFieldId fieldId;

      throw Exception("Not implemented yet!");

      return spBcs;
   }

   void BoussinesqSphereModel::setInitialState(SharedSimulation spSim)
   {
      // Field IDs iterator
      std::vector<GeoMHDiSCC::PhysicalNames::Id>::const_iterator  it;
      std::vector<GeoMHDiSCC::PhysicalNames::Id> ids = BoussinesqSphereModel::fieldIds();

      // Create and add initial state file to IO
      IoVariable::SharedStateFileReader spInit(new IoVariable::StateFileReader("_initial", SchemeType::type(), SchemeType::isRegular()));

      // Set expected field names
      for(it = ids.begin(); it != ids.end(); ++it)
      {
         spInit->expect(*it);
      }

      // Set simulation state
      spSim->setInitialState(spInit);
   }

}
