/** 
 * @file BoussinesqCylinderModel.cpp
 * @brief Source of the Boussinesq cylinder physical model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "PhysicalModels/BoussinesqCylinderModel.hpp"

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Generator/Visualizers/FieldVisualizer.hpp"

namespace GeoMHDiSCC {

   std::vector<PhysicalNames::Id> BoussinesqCylinderModel::fieldIds()
   {
      // Create storage
      std::vector<PhysicalNames::Id> ids;

      // Add temperature
      ids.push_back(PhysicalNames::TEMPERATURE);

      return ids;
   }

   std::vector<NonDimensional::Id> BoussinesqCylinderModel::paramIds()
   {
      // Create storage
      std::vector<NonDimensional::Id> ids;

      return ids;
   }

   std::vector<bool> BoussinesqCylinderModel::isPeriodicBox()
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

   void BoussinesqCylinderModel::addEquations(SharedSimulation spSim)
   {
      throw Exception("Not implemented yet!");
   }

   void BoussinesqCylinderModel::addStates(SharedStateGenerator spGen)
   {
      // Shared pointer to equation
      //Equations::SharedCylinderExactScalarState spExact;

      // Add initial state generator
      //spExact = spGen->addScalarEquation<Equations::AnnulusExactScalarState>();
      //spExact->setIdentity(PhysicalNames::TEMPERATURE);
      //spExact->setStateType(Equations::AnnulusExactScalarState::CONSTANT);

      // Add output file
      IoVariable::SharedStateFileWriter spOut(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      //spOut->expect(PhysicalNames::TEMPERATURE);
      spGen->addOutputFile(spOut);
   }

   void BoussinesqCylinderModel::addVisualizers(SharedVisualizationGenerator spVis)
   {
      // Shared pointer to basic field visualizer
      Equations::SharedFieldVisualizer spField;

      // Add first field visualization
      //spField = spVis->addScalarEquation<Equations::FieldVisualizer>();
      //spField->setFields(true, false);
      //spField->setIdentity(PhysicalNames::TEMPERATURE);

      // Add output file
      IoVariable::SharedVisualizationFileWriter spOut(new IoVariable::VisualizationFileWriter(SchemeType::type()));
      //spOut->expect(PhysicalNames::TEMPERATURE);
      spVis->addOutputFile(spOut);
   }

   void BoussinesqCylinderModel::setVisualizationState(SharedVisualizationGenerator spVis)
   {
      // Create and add initial state file to IO
      IoVariable::SharedStateFileReader spIn(new IoVariable::StateFileReader("4Visu", SchemeType::type(), SchemeType::isRegular()));

      // Set expected fields
      //spIn->expect(PhysicalNames::TEMPERATURE);

      // Set simulation state
      spVis->setInitialState(spIn);
   }

   void BoussinesqCylinderModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
      
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
   }

   void BoussinesqCylinderModel::addHdf5OutputFiles(SharedSimulation spSim)
   {
      // Field IDs iterator
      std::vector<GeoMHDiSCC::PhysicalNames::Id>::const_iterator  it;
      std::vector<GeoMHDiSCC::PhysicalNames::Id> ids = BoussinesqCylinderModel::fieldIds();

      // Create and add state file to IO
      IoVariable::SharedStateFileWriter spState(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      for(it = ids.begin(); it != ids.end(); ++it)
      {
         spState->expect(*it);
      }
      spSim->addOutputFile(spState);
   }

   SharedSimulationBoundary BoussinesqCylinderModel::createBoundary(const std::map<std::string,int>& bcIds)
   {
      // Create shared simulation boundary
      SharedSimulationBoundary  spBcs(new SimulationBoundary());

      // Storage for the dimension ID
      Dimensions::Simulation::Id dimId;

      // Create equation and field keys
      SpectralFieldId eqId;
      SpectralFieldId fieldId;

      return spBcs;
   }

   void BoussinesqCylinderModel::setInitialState(SharedSimulation spSim)
   {
      // Field IDs iterator
      std::vector<GeoMHDiSCC::PhysicalNames::Id>::const_iterator  it;
      std::vector<GeoMHDiSCC::PhysicalNames::Id> ids = BoussinesqCylinderModel::fieldIds();

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
